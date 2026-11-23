// Test the gapped-plate network against measured retention and conductivity.
//
//   [CORNER_EXP=n] node tools/validate.mjs <DataFiles_v2> [cache.json] [κ] [textures] [film|nofilm|macro]
//
// Data: Hohenbrink et al. (2023), doi:10.5880/fidgeo.2023.012 (CC BY 4.0). Not stored in this repo.
//
// For each sand or loamy-sand sample the network sees only the measured retention curve.
// Its pore-size spread σ and one length scale are fitted to θ(h); conductivity is then
// predicted with no conductivity data and compared with the measured K(h) points.
import fs from "node:fs";
import path from "node:path";
import { DEFAULTS, HEAD_CONST, buildNetwork, computeCurves, filmSaturation, filmK, macroSaturation, macroK } from "../networks/model.js";

const dir = process.argv[2];
if (!dir) { console.error("usage: node tools/validate.mjs <DataFiles_v2> [cache.json]"); process.exit(1); }
const cachePath = process.argv[3];
const KAPPA = process.argv[4] !== undefined && process.argv[4] !== "" ? Number(process.argv[4]) : DEFAULTS.corner;
// USDA texture classes to test; the corner constant was fitted on Sa and LoSa
const TEXTURES = (process.argv[5] || "Sa,LoSa").split(",");
// "film": fit non-capillary film water θ_nc to retention (default); "nofilm": the earlier fit,
// capillary network only, with van Genuchten's θ_r taken as given
// "macro": film water plus a structural macropore set (volume w, median radius r_m), both fitted to retention
const MODE = process.argv[6] || "film";
const FILM = MODE === "film" || MODE === "macro";
const MACRO = MODE === "macro";
// exponent of the corner conductance (r*/r)^n in drained pores; 4 for triangular corners (Tuller & Or 2001)
const CEXP = Number(process.env.CORNER_EXP || 4);
// candidate macropore air-entry heads (mm): r_m = HEAD_CONST / h_m, 0.15 mm … 3 mm
const HM = Array.from({ length: 14 }, (_, i) => 5 * 10 ** (i * 0.1));

/* ───── data ───── */
function readCsv(name) {
  const [head, ...rows] = fs.readFileSync(path.join(dir, name), "utf8").trim().split(/\r?\n/);
  const cols = head.split(",");
  return rows.map((line) => Object.fromEntries(line.split(",").map((v, i) => [cols[i], v])));
}
const basic = new Map(readCsv("BasicProp.csv").map((r) => [r.Sample_ID, r]));
const param = new Map(readCsv("Param.csv").map((r) => [r.Sample_ID, r]));
const group = (rows) => {
  const m = new Map();
  for (const r of rows) { if (!m.has(r.Sample_ID)) m.set(r.Sample_ID, []); m.get(r.Sample_ID).push(r); }
  return m;
};
const ret = group(readCsv("RetMeas.csv")), cond = group(readCsv("CondMeas.csv"));

const CMD = 8.64e6; // m/s → cm/d
const samples = [];
for (const [id, b] of basic) {
  if (!TEXTURES.includes(b.TexClass_USDA)) continue;
  const pr = param.get(id), kRows = cond.get(id) || [], rRows = ret.get(id) || [];
  const thr = +pr.thr_VGM, ths = +pr.ths_VGM;
  const K = kRows.filter((r) => r.MeasType === "HyProp" && +r.k > 0).map((r) => ({ h: 10 ** +r.pF * 10, k: +r.k / CMD }));
  const ks = kRows.filter((r) => r.MeasType === "KSAT" && +r.k > 0).map((r) => +r.k / CMD);
  let R = rRows.filter((r) => (r.MeasType === "HyProp" || r.MeasType === "WP4") && +r.pF > 0.3)
    .map((r) => ({ h: 10 ** +r.pF * 10, th: +r.theta, Se: Math.min(1, Math.max(0, (+r.theta - thr) / (ths - thr))) }));
  if (K.length < 5 || R.length < 10 || !(ths - thr > 0.1)) continue;
  R = R.filter((_, i) => i % Math.ceil(R.length / 60) === 0);
  samples.push({
    id, tex: b.TexClass_USDA, phi: ths - thr, ths, thr, K, R,
    Ksat: ks.length ? Math.exp(ks.reduce((a, v) => a + Math.log(v), 0) / ks.length) : NaN,
    vg: { al: +pr.al_VGM / 10, n: +pr.n_VGM, ks: +pr.ks_VGM / CMD, lam: +pr.lambda_VGM },
  });
}

/* ───── network library ───── */
// Everything in the network scales with one length: h ∝ 1/s and K ∝ s². So one solve per
// (φ, σ) covers every sample; the scale s is fitted afterwards.
// Geometry keeps the calibrated proportions in terms of the volume-median pore radius r50v
// (1.87 r̄ for the calibrated truncated exponential): gap b = 1.39 r50v, plate t = 5.35 r50v / voidScale.
const PHIS = [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55];
const SIGS = [0.25, 0.4, 0.55, 0.7, 0.9, 1.1, 1.4];
function libraryEntry(phi, sig) {
  const rbar = 0.05, r50v = rbar * Math.exp(1.5 * sig * sig);
  const vs = (phi / (1 - phi)) / (0.35 / 0.65);
  const p = {
    ...DEFAULTS, d50: 1, voidScale: false, dist: "lognormal", sigln: sig, bSigma: sig, phi,
    rRatio: rbar, rmaxRatio: Math.exp(3 * sig), bRatio: 1.39 * r50v, tRatio: (5.35 * r50v) / vs,
    N: 300, layers: 14, steps: 40, seed: 3, corner: KAPPA, cornerExp: CEXP,
  };
  const out = computeCurves(buildNetwork(p, "gap"), p, { keepFields: false });
  return { phi, sig, h: out.steps.map((d) => d.h), S: out.steps.map((d) => d.S), K: out.steps.map((d) => d.Kv), aniso: out.summary.aniso };
}
const lib = cachePath && fs.existsSync(cachePath) ? JSON.parse(fs.readFileSync(cachePath, "utf8")) : [];
let added = 0;
for (const phi of PHIS)
  for (const sig of SIGS)
    if (!lib.some((x) => x.phi === phi && x.sig === sig)) { lib.push(libraryEntry(phi, sig)); added++; process.stderr.write("."); }
if (added && cachePath) fs.writeFileSync(cachePath, JSON.stringify(lib));

/** Interpolate a library curve at head h (mm). Returns [S, K] or null beyond the driest step. */
function at(e, h) {
  if (h <= e.h[1]) return [1, e.K[0]];
  const n = e.h.length;
  if (h > e.h[n - 1]) return null;
  let i = 2;
  while (e.h[i] < h) i++;
  const f = Math.log(h / e.h[i - 1]) / Math.log(e.h[i] / e.h[i - 1]);
  return [e.S[i - 1] + f * (e.S[i] - e.S[i - 1]), Math.exp(Math.log(e.K[i - 1]) + f * (Math.log(e.K[i]) - Math.log(e.K[i - 1])))];
}

/* ───── baselines ───── */
function vgKr({ al, n }, h, lam) {
  const m = 1 - 1 / n, Se = (1 + (al * h) ** n) ** -m;
  return Se ** lam * (1 - (1 - Se ** (1 / m)) ** m) ** 2;
}
const rmse = (a) => Math.sqrt(a.reduce((s, v) => s + v * v, 0) / a.length);
const median = (a) => { const s = [...a].sort((x, y) => x - y); return s.length ? (s[(s.length - 1) >> 1] + s[s.length >> 1]) / 2 : NaN; };

/* ───── fit and predict ───── */
const rows = [];
for (const s of samples) {
  let best = null;
  if (!FILM) {
    const phi = PHIS.reduce((a, v) => (Math.abs(v - s.phi) < Math.abs(a - s.phi) ? v : a));
    for (const e of lib.filter((x) => x.phi === phi)) {
      for (let ls = -4; ls <= 4; ls += 0.02) {
        const sc = Math.exp(ls);
        let err = 0;
        for (const pt of s.R) { const v = at(e, pt.h * sc); err += (s.thr + (s.ths - s.thr) * (v ? v[0] : e.S.at(-1)) - pt.th) ** 2; }
        err = Math.sqrt(err / s.R.length);
        if (!best || err < best.err) best = { err, e, sc, tnc: 0 };
      }
    }
  } else {
    // θ = θ_s·S_cap + θ_nc·(S_nc − S_cap) is linear in θ_nc, so it is solved exactly for each
    // (network, scale); the capillary porosity θ_s − θ_nc must match the network's own φ
    const Sc = new Float64Array(s.R.length), Sn = new Float64Array(s.R.length), Sm = new Float64Array(s.R.length);
    for (const e of lib) {
      for (let ls = -4; ls <= 4; ls += MACRO ? 0.04 : 0.02) {
        const sc = Math.exp(ls), ha = e.h[1] / sc;
        let num = 0, den = 0;
        s.R.forEach((pt, k) => {
          const v = at(e, pt.h * sc);
          Sc[k] = v ? v[0] : e.S.at(-1);
          Sn[k] = filmSaturation(pt.h, ha);
          const d = Sn[k] - Sc[k];
          num += (pt.th - s.ths * Sc[k]) * d; den += d * d;
        });
        if (!MACRO) {
          const tnc = Math.min(0.8 * s.ths, Math.max(0, den > 0 ? num / den : 0));
          if (Math.abs(e.phi - (s.ths - tnc)) > 0.026) continue;
          let err = 0;
          s.R.forEach((pt, k) => { err += (s.ths * Sc[k] + tnc * (Sn[k] - Sc[k]) - pt.th) ** 2; });
          err = Math.sqrt(err / s.R.length);
          if (!best || err < best.err) best = { err, e, sc, tnc, w: 0, rm: 0 };
          continue;
        }
        // θ = θ_s·S_cap + θ_nc·(S_nc − S_cap) + w·(S_ma − S_cap): linear in (θ_nc, w)
        for (const hm of HM) {
          if (hm > ha / 2) break; // macropores must be coarser than the grain network
          const rm = HEAD_CONST / hm;
          let a11 = 0, a12 = 0, a22 = 0, b1 = 0, b2 = 0;
          s.R.forEach((pt, k) => {
            const u = Sn[k] - Sc[k], m = macroSaturation(pt.h, rm) - Sc[k], y = pt.th - s.ths * Sc[k];
            Sm[k] = m;
            a11 += u * u; a12 += u * m; a22 += m * m; b1 += u * y; b2 += m * y;
          });
          const cands = [];
          const det = a11 * a22 - a12 * a12;
          if (det > 1e-12) cands.push([(b1 * a22 - b2 * a12) / det, (a11 * b2 - a12 * b1) / det]);
          cands.push([a11 > 0 ? b1 / a11 : 0, 0], [0, a22 > 0 ? b2 / a22 : 0]);
          for (let [tnc, w] of cands) {
            tnc = Math.min(0.8 * s.ths, Math.max(0, tnc)); w = Math.min(0.3, Math.max(0, w));
            if (Math.abs(e.phi - (s.ths - tnc - w)) > 0.026) continue;
            let err = 0;
            s.R.forEach((pt, k) => { err += (s.ths * Sc[k] + tnc * (Sn[k] - Sc[k]) + w * Sm[k] - pt.th) ** 2; });
            err = Math.sqrt(err / s.R.length);
            if (!best || err < best.err) best = { err, e, sc, tnc, w, rm };
          }
        }
      }
    }
    if (!best) continue;
  }
  const { e, sc, tnc, w = 0, rm = 0 } = best;
  const ha = e.h[1] / sc;
  const netAbs = [], netShape = [], muKs = [], muFit = [], seAt = [];
  const v_at = (h) => { const m = 1 - 1 / s.vg.n; return (1 + (s.vg.al * h) ** s.vg.n) ** -m; };
  let shift = 0, used = 0;
  const pts = [];
  for (const pt of s.K) {
    const v = at(e, pt.h * sc);
    if (!v) continue;
    const kf = (FILM ? filmK(pt.h, tnc * filmSaturation(pt.h, ha)) : 0) + (w > 0 ? macroK(pt.h, w, rm) : 0);
    pts.push([Math.log10((sc * sc * v[1] + kf) / pt.k), pt]);
  }
  if (pts.length < 5) continue;
  for (const [d] of pts) { netAbs.push(d); shift += d; used++; }
  shift /= used;
  for (const [d, pt] of pts) {
    netShape.push(d - shift); // curve shape only: best vertical shift removed
    seAt.push(v_at(pt.h));
    if (s.Ksat > 0) muKs.push(Math.log10((s.Ksat * vgKr(s.vg, pt.h, 0.5)) / pt.k));
    muFit.push(Math.log10((s.vg.ks * vgKr(s.vg, pt.h, s.vg.lam)) / pt.k));
  }
  // slope of the shape residual against log10 h: > 0 means the network falls too slowly with drying
  const slope = (() => { const X = pts.map(([, pt]) => Math.log10(pt.h)), mx = X.reduce((a, v) => a + v, 0) / X.length; let sxy = 0, sxx = 0; pts.forEach(([d], k) => { sxy += (X[k] - mx) * (d - shift); sxx += (X[k] - mx) ** 2; }); return sxx > 0 ? sxy / sxx : NaN; })();
  const hMed = median(pts.map(([, pt]) => Math.log10(pt.h / 10)));
  const muShape = (() => { const d = pts.map(([, pt]) => Math.log10(vgKr(s.vg, pt.h, 0.5) / pt.k)); const m = d.reduce((a, v) => a + v, 0) / d.length; return rmse(d.map((v) => v - m)); })();
  rows.push({
    id: s.id, tex: s.tex, phi: s.phi, sig: e.sig, retErr: best.err, n: used, tnc, w, rm,
    r50v_um: 1000 * sc * 0.05 * Math.exp(1.5 * e.sig * e.sig),
    netAbs: rmse(netAbs), bias: shift, netShape: rmse(netShape), muShape, slope, hMed,
    muKs: muKs.length ? rmse(muKs) : NaN, muFit: rmse(muFit),
    muBias: muKs.length ? muKs.reduce((a, v) => a + v, 0) / muKs.length : NaN, se: median(seAt),
    ksatRatio: s.Ksat > 0 ? Math.log10((sc * sc * e.K[0] + (w > 0 ? macroK(0, w, rm) : 0)) / s.Ksat) : NaN,
  });
}

const f = (v, d = 2) => (Number.isFinite(v) ? v.toFixed(d) : "–");
console.log("sample      tex   φ    σ    fit(S) pts  net|abs  bias   net|shape  Mualem|shape  Mualem+Ks  VGM fitted  log Ksat net/meas");
for (const r of rows)
  console.log(`${r.id.padEnd(11)} ${r.tex.padEnd(5)} ${f(r.phi)} ${f(r.sig)} ${f(r.retErr, 3)}  ${String(r.n).padStart(3)}  ${f(r.netAbs).padStart(6)} ${f(r.bias).padStart(6)}  ${f(r.netShape).padStart(8)}  ${f(r.muShape).padStart(11)}  ${f(r.muKs).padStart(9)}  ${f(r.muFit).padStart(9)}  ${f(r.ksatRatio).padStart(8)}`);
const col = (k) => rows.map((r) => r[k]).filter(Number.isFinite);
console.log(`\n${rows.length} samples (of ${samples.length} sand / loamy sand with ≥5 K points). Medians, in decades of log10 K:`);
console.log(`  film water: ${FILM ? "fitted" : "off"}`);
console.log(`  retention fit RMSE (θ, volumetric)           ${f(median(col("retErr")), 3)}`);
console.log(`  network, absolute, no K data used         ${f(median(col("netAbs")))}   (median bias ${f(median(col("bias")))})`);
console.log(`  Mualem ℓ=0.5 scaled to measured K_sat     ${f(median(col("muKs")))}   (${col("muKs").length} samples with K_sat)`);
console.log(`  network, shape only (best shift removed)  ${f(median(col("netShape")))}`);
console.log(`  Mualem ℓ=0.5, shape only                  ${f(median(col("muShape")))}`);
console.log(`  VGM with K_s and ℓ fitted to the K data   ${f(median(col("muFit")))}   (lower bound)`);
console.log(`  Mualem+K_sat median bias ${f(median(col("muBias")))}; median saturation at the K points ${f(median(col("se")))}; κ = ${KAPPA}, corner exponent ${CEXP}`);
for (const tex of TEXTURES) {
  const sub = rows.filter((r) => r.tex === tex);
  if (!sub.length) continue;
  const b = median(sub.map((r) => r.bias));
  // bias moves ~0.95 decades per decade of κ (measured on the sands), so this κ would remove it
  console.log(`  ${tex.padEnd(6)} ${String(sub.length).padStart(3)} samples, θ fit ${f(median(sub.map((r) => r.retErr)), 3)}, θ_nc ${f(median(sub.map((r) => r.tnc)), 3)}${MACRO ? `, w ${f(median(sub.map((r) => r.w)), 3)} (${sub.filter((r) => r.w > 0).length} >0), r_m ${f(median(sub.filter((r) => r.w > 0).map((r) => r.rm)), 2)} mm` : ""}, bias ${f(b)}, absolute ${f(median(sub.map((r) => r.netAbs)))}, shape ${f(median(sub.map((r) => r.netShape)))} (slope ${f(median(sub.map((r) => r.slope)))}/decade, K points at pF ${f(median(sub.map((r) => r.hMed)), 1)}), Mualem+K_sat ${f(median(sub.map((r) => r.muKs).filter(Number.isFinite)))}, σ fit ${f(median(sub.map((r) => r.sig)))}, κ that removes bias ≈ ${(KAPPA * 10 ** (-b / 0.95)).toExponential(1)}`);
}
for (const lab of [...new Set(rows.map((r) => r.id.split("_")[0]))]) {
  const sub = rows.filter((r) => r.id.startsWith(lab));
  console.log(`  ${lab}: ${sub.length} samples, median bias ${f(median(sub.map((r) => r.bias)))}, absolute ${f(median(sub.map((r) => r.netAbs)))}, shape ${f(median(sub.map((r) => r.netShape)))}`);
}
console.log(`  log10 K_sat network / measured            ${f(median(col("ksatRatio")))}   (${col("ksatRatio").length} samples)`);
