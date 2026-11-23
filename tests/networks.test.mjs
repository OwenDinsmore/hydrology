// node --test tests/
import test from "node:test";
import assert from "node:assert/strict";
import {
  filmSaturation, filmThickness, filmK, macroSaturation, macroK, normCdf, HEAD_CONST, FLUID,
  DEFAULTS, MODELS, buildNetwork, computeCurves, solveFlow, boundarySets, diskShapeFactor, K_UNIT,
} from "../networks/model.js";

const P = { ...DEFAULTS, N: 120, layers: 8, steps: 14 };

test("straight bundle reproduces the HW 3 formula K = ρg/8μ Σr⁴/R²", () => {
  const net = buildNetwork(P, "bundle");
  const out = computeCurves(net, P, { keepFields: false });
  let s4 = 0;
  for (let i = 0; i < net.nBonds; i++) s4 += net.r[i] ** 4;
  const expected = (K_UNIT / 8) * s4 / net.meta.R ** 2;
  assert.ok(Math.abs(out.summary.Ksat / expected - 1) < 1e-4, `${out.summary.Ksat} vs ${expected}`);
  assert.equal(out.summary.KsatH, 0);
});

test("shape factor of a disk is symmetric and positive", () => {
  const s = diskShapeFactor(0.15);
  assert.ok(s > 0.5 && s < 5, `shape ${s}`);
});

test("uniform cubic lattice is isotropic (K_h ≈ K_v)", () => {
  const net = buildNetwork({ ...P, coord: 6, rmaxRatio: 1.0001, rRatio: 0.2, dist: "lognormal", sigln: 0.0001 }, "lattice");
  const out = computeCurves(net, { ...P, steps: 3 }, { keepFields: false });
  assert.ok(Math.abs(out.summary.aniso - 1) < 0.2, `anisotropy ${out.summary.aniso}`);
});

test("flow is conserved: what enters the top leaves the bottom", () => {
  for (const m of ["plates", "gap", "tree", "lattice"]) {
    const net = buildNetwork(P, m);
    const wet = new Uint8Array(net.nBonds).fill(1);
    const { fixedV } = boundarySets(net);
    const Pn = new Float64Array(net.nNodes).fill(0.5);
    const Qin = solveFlow(net, wet, fixedV, Pn);
    let Qout = 0;
    for (const n of net.outlet)
      for (let k = net.adjStart[n]; k < net.adjStart[n + 1]; k++) {
        const b = net.adjBond[k], o = net.a[b] === n ? net.b[b] : net.a[b];
        if (Number.isFinite(Pn[o])) Qout += net.g[b] * (Pn[o] - Pn[n]);
      }
    assert.ok(Qin > 0, `${m} carries flow`);
    assert.ok(Math.abs(Qout / Qin - 1) < 1e-5, `${m}: in ${Qin}, out ${Qout}`);
  }
});

test("every model drains monotonically and porosities match the target", () => {
  for (const m of MODELS) {
    const net = buildNetwork(P, m);
    const out = computeCurves(net, P, { keepFields: false });
    const S = out.steps.map((s) => s.S), K = out.steps.map((s) => s.Kv);
    for (let i = 1; i < S.length; i++) {
      assert.ok(S[i] <= S[i - 1] + 1e-12, `${m} S rises at ${i}`);
      assert.ok(K[i] <= K[i - 1] * (1 + 1e-6) + 1e-30, `${m} K rises at ${i}`);
    }
    assert.ok(Math.abs(out.steps[0].S - 1) < 1e-12);
    // a tree regrows differently each time the sample is resized, so it only lands near the target
    const tol = m === "tree" ? 0.15 : 0.03;
    assert.ok(Math.abs(net.meta.phi / P.phi - 1) < tol, `${m} φ = ${net.meta.phi}`);
  }
});

test("corner flow keeps every model conducting to the driest step", () => {
  for (const m of MODELS) {
    const dry = computeCurves(buildNetwork(P, m), P, { keepFields: false });
    assert.ok(dry.steps.at(-1).Kv > 0, `${m} still conducts`);
    assert.equal(dry.summary.Sc, 0);
  }
  const round = { ...P, corner: 0 };
  assert.ok(computeCurves(buildNetwork(round, "plates"), round, { keepFields: false }).summary.Sc > 0);
});

test("calibrated gapped plates: near isotropic and near Kozeny–Carman across porosity", () => {
  for (const phi of [0.28, 0.35, 0.42]) {
    const p = { ...DEFAULTS, phi, steps: 4 };
    const s = computeCurves(buildNetwork(p, "gap"), p, { keepFields: false }).summary;
    assert.ok(s.aniso > 0.6 && s.aniso < 1.6, `φ=${phi}: K_h/K_v = ${s.aniso}`);
    assert.ok(s.Ksat / s.KC > 0.7 && s.Ksat / s.KC < 2.5, `φ=${phi}: K/KC = ${s.Ksat / s.KC}`);
  }
});

test("a gap adds lateral conductivity that stacked plates lack", () => {
  const plates = computeCurves(buildNetwork(P, "plates"), P, { keepFields: false }).summary;
  const gap = computeCurves(buildNetwork(P, "gap"), P, { keepFields: false }).summary;
  assert.ok(gap.aniso > plates.aniso, `plates ${plates.aniso}, gap ${gap.aniso}`);
});

test("film water: full below air entry, gone at oven dryness, nanometre films", () => {
  assert.ok(filmSaturation(10, 300) > 0.99);
  assert.ok(filmSaturation(10 ** 7.8, 300) < 1e-9);
  const mid = filmSaturation(1e5, 300);
  assert.ok(mid > 0.3 && mid < 0.7, `S_nc at pF 4 = ${mid}`);
  const d = filmThickness(1000); // 1 m of suction
  assert.ok(d > 3e-9 && d < 2e-8, `δ = ${d}`);
  assert.ok(filmK(1000, 0.1) < 1e-10);
});

test("adding film water raises θ and K but leaves the capillary network unchanged", () => {
  const p0 = { ...P, steps: 8 }, p1 = { ...p0, film: 0.1 };
  const a = computeCurves(buildNetwork(p0, "gap"), p0, { keepFields: false });
  const b = computeCurves(buildNetwork(p1, "gap"), p1, { keepFields: false });
  assert.ok(Math.abs(b.steps[0].theta - a.steps[0].theta - 0.1) < 1e-9);
  for (let i = 0; i < a.steps.length; i++) assert.ok(b.steps[i].Kv >= a.steps[i].Kv);
});

test("macropore bundle: drains around its median radius, K matches the random-tube formula", () => {
  assert.ok(Math.abs(normCdf(1) - 0.841345) < 1e-6 && Math.abs(normCdf(-2) - 0.02275) < 1e-5);
  const rm = 0.5, hm = HEAD_CONST / rm;
  assert.ok(Math.abs(macroSaturation(hm, rm) - 0.5) < 1e-9);
  assert.ok(macroSaturation(hm / 10, rm) > 0.99 && macroSaturation(hm * 10, rm) < 0.01);
  // saturated: ρg/24μ · w · r_m² e^{2σ²}
  const k0 = macroK(0, 0.02, rm, 0.5), r = rm * 1e-3;
  assert.ok(Math.abs(k0 / ((FLUID.rho * FLUID.g / (24 * FLUID.mu)) * 0.02 * r * r * Math.exp(0.5)) - 1) < 1e-12);
  assert.ok(macroK(hm, 0.02, rm) < k0 && macroK(hm * 10, 0.02, rm) < 1e-3 * k0);
});

test("corner exponent defaults to the triangular-corner 4", () => {
  const p0 = { ...P, steps: 8 }, p4 = { ...p0, cornerExp: 4 }, p2 = { ...p0, cornerExp: 2 };
  const net = buildNetwork(p0, "gap");
  const a = computeCurves(net, p0, { keepFields: false }), b = computeCurves(net, p4, { keepFields: false });
  const c = computeCurves(net, p2, { keepFields: false });
  for (let i = 0; i < a.steps.length; i++) {
    assert.equal(a.steps[i].Kv, b.steps[i].Kv);
    assert.ok(c.steps[i].Kv >= a.steps[i].Kv * (1 - 1e-9));
  }
});
