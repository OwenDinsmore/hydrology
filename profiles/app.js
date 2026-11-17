import { DEFAULTS, TEXTURES, CASES, soil, qRef, qMax, profile, thetaAt, thetaE } from "./model.js";

const MULTIPLES = ["sand", "loam", "clay"];
const SEC_PER_DAY = 86400;

/* ───────────────────────── model run ───────────────────────── */
// Solve the five cases. ref = "KL" (fractions of hydrostatic K at the surface, as in the notes) or "qmax".
function run(input, ref) {
  const p = soil(input);
  const kl = qRef(p);
  const qm = qMax(p);
  const base = ref === "qmax" ? qm : kl;
  const profiles = CASES.map((cs) => {
    const a = cs.ofRef ? cs.frac * base : cs.frac;
    const label = cs.ofRef ? `${cs.frac === 0.5 ? "½" : ".97"} ${ref === "qmax" ? "q_max" : "K(L)"}` : cs.label;
    return { ...cs, label, short: shortLabel(cs), ...profile(p, a) };
  });
  return { p, input, ref, kl, qm, profiles };
}
function shortLabel(cs) {
  return { "perc-half": "−½", "perc-quarter": "−¼", hydro: "0", "evap-half": "+½", "evap-97": "+.97" }[cs.key];
}

/* ───────────────────────── formatting ───────────────────────── */
const SUP = { "-": "⁻", 0: "⁰", 1: "¹", 2: "²", 3: "³", 4: "⁴", 5: "⁵", 6: "⁶", 7: "⁷", 8: "⁸", 9: "⁹" };
const fmt = {
  sci(x, n) {
    const [m, e] = x.toExponential(n - 1).split("e");
    const pow = "10" + [...String(Number(e))].map((ch) => SUP[ch]).join("");
    return Number(m) === 1 ? pow : `${m} × ${pow}`;
  },
  sig(x, n = 3) {
    if (!isFinite(x)) return "–";
    if (x === 0) return "0";
    const a = Math.abs(x);
    if (a >= 1e5 || a < 1e-2) return fmt.sci(x, n);
    return Number(x.toPrecision(n)).toString();
  },
  f(x, d = 2) { return isFinite(x) ? x.toFixed(d) : "–"; },
  minus(s) { return String(s).replace(/^-/, "−"); },
};
const sub = (s) => s.replace("K_s", "K<sub>s</sub>").replace("q_max", "q<sub>max</sub>");

/* ───────────────────────── chart ───────────────────────── */
const svgNS = "http://www.w3.org/2000/svg";
function el(name, attrs = {}, parent) {
  const e = document.createElementNS(svgNS, name);
  for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v);
  if (parent) parent.appendChild(e);
  return e;
}
function niceStep(range, target) {
  const raw = range / target;
  const mag = Math.pow(10, Math.floor(Math.log10(raw)));
  for (const m of [1, 2, 2.5, 5, 10]) if (m * mag >= raw) return m * mag;
  return 10 * mag;
}

class ProfileChart {
  constructor(svg, tooltip, { mini = false } = {}) {
    this.svg = svg;
    this.tooltip = tooltip;
    this.mini = mini;
    this.m = mini ? { t: 22, r: 10, b: 30, l: 40 } : { t: 30, r: 16, b: 40, l: 56 };
    this.gBand = el("g", {}, svg);
    this.gGrid = el("g", { class: "grid" }, svg);
    this.gRef = el("g", {}, svg);
    this.gSeries = el("g", {}, svg);
    this.gAxis = el("g", { class: "axis" }, svg);
    this.gTicks = el("g", { class: "tick" }, svg);
    this.gLabels = el("g", {}, svg);
    this.gHover = el("g", { visibility: "hidden" }, svg);
    this.hit = el("rect", { class: "hit" }, svg);
    this.xText = el("text", { class: "axis-label", "text-anchor": "middle" }, svg);
    this.xText.innerHTML = `Θ<tspan dy="3" font-size="10">E</tspan>`;
    this.yText = el("text", { class: "axis-label", "text-anchor": "middle" }, svg);
    this.yText.textContent = mini ? "z (cm)" : "height above water table z (cm)";
    if (tooltip) {
      this.hit.addEventListener("pointermove", (e) => this.onHover(e));
      this.hit.addEventListener("pointerleave", () => this.hideHover());
    }
  }
  x(v) { return this.m.l + v * (this.W - this.m.l - this.m.r); }
  y(v) { return this.H - this.m.b - (v / this.L) * (this.H - this.m.t - this.m.b); }
  draw(res) {
    this.res = res;
    this.L = res.p.L;
    this.W = this.svg.clientWidth || 600;
    this.H = this.svg.clientHeight || 400;
    const { m, W, H } = this;
    const x0 = m.l, x1 = W - m.r, y0 = m.t, y1 = H - m.b;
    this.svg.setAttribute("viewBox", `0 0 ${W} ${H}`);
    for (const [k, v] of Object.entries({ x: x0, y: y0, width: x1 - x0, height: y1 - y0 })) this.hit.setAttribute(k, v);
    this.xText.setAttribute("x", (x0 + x1) / 2);
    this.xText.setAttribute("y", H - 5);
    this.yText.setAttribute("transform", `translate(12 ${(y0 + y1) / 2}) rotate(-90)`);
    for (const g of [this.gBand, this.gGrid, this.gRef, this.gSeries, this.gAxis, this.gTicks, this.gLabels]) g.innerHTML = "";

    // grid and ticks
    const xStep = this.mini ? 0.5 : 0.2;
    for (let v = 0; v <= 1 + 1e-9; v += xStep) {
      const xx = this.x(v);
      el("line", { x1: xx, x2: xx, y1: y0, y2: y1 }, this.gGrid);
      el("text", { x: xx, y: y1 + 16, "text-anchor": "middle" }, this.gTicks).textContent = Number(v.toFixed(1));
    }
    const zStep = niceStep(this.L, this.mini ? 4 : 8);
    for (let v = 0; v <= this.L + 1e-9; v += zStep) {
      const yy = this.y(v);
      el("line", { x1: x0, x2: x1, y1: yy, y2: yy }, this.gGrid);
      el("text", { x: x0 - 7, y: yy + 4, "text-anchor": "end" }, this.gTicks).textContent = Number(v.toFixed(1));
    }
    el("path", { d: `M${x0} ${y0} V${y1} H${x1}` }, this.gAxis);

    // surface and air-entry height
    el("line", { class: "surface", x1: x0, x2: x1, y1: this.y(this.L), y2: this.y(this.L) }, this.gRef);
    if (!this.mini) {
      el("text", { class: "ref-label", x: x0 + 5, y: this.y(this.L) + 13 }, this.gRef).textContent = "surface";
      el("text", { class: "ref-label", x: x1 - 4, y: y1 - 5, "text-anchor": "end" }, this.gRef).textContent = "water table";
    }
    const sA = res.p.sA;
    if (sA < this.L) {
      const yy = this.y(sA);
      el("line", { class: "ref-line", x1: x0, x2: x1, y1: yy, y2: yy }, this.gRef);
      if (!this.mini) el("text", { class: "ref-label", x: x0 + 5, y: yy - 4 }, this.gRef).textContent = `z = |h_A| = ${fmt.sig(sA, 3)}`;
    }

    // dry band where an upward flux has no steady solution
    const dry = res.profiles.filter((pr) => pr.dryAt != null).map((pr) => pr.dryAt);
    if (dry.length) {
      const zd = Math.min(...dry);
      el("rect", { class: "dry-band", x: x0, y: y0, width: x1 - x0, height: this.y(zd) - y0 }, this.gBand);
      el("text", { class: "ref-label", x: x1 - 4, y: this.y(zd) - 4, "text-anchor": "end" }, this.gBand).textContent = "q > q_max: dries out";
    }

    // series
    for (const pr of res.profiles) {
      const d = pr.pts.filter((q) => q.z <= this.L + 1e-9).map((q) => `${this.x(q.theta).toFixed(2)} ${this.y(q.z).toFixed(2)}`);
      el("path", { class: `series s-${pr.key}`, d: "M" + d.join(" L") }, this.gSeries);
      const top = pr.pts.at(-1);
      el("circle", { class: `end-dot dot-${pr.key}`, cx: this.x(top.theta), cy: this.y(top.z), r: this.mini ? 3 : 4 }, this.gSeries);
    }

    // direct labels above the surface, pushed apart so they never overlap
    if (!this.mini) {
      const gap = 30;
      const items = res.profiles
        .map((pr) => ({ pr, x: this.x(pr.pts.at(-1).theta), top: pr.pts.at(-1) }))
        .sort((a, b) => a.x - b.x);
      const pos = items.map((it) => it.x);
      for (let iter = 0; iter < 50; iter++) {
        for (let i = 1; i < pos.length; i++) {
          const over = gap - (pos[i] - pos[i - 1]);
          if (over > 0) { pos[i - 1] -= over / 2; pos[i] += over / 2; }
        }
        for (let i = 0; i < pos.length; i++) pos[i] = Math.min(x1 - 12, Math.max(x0 + 12, pos[i]));
      }
      items.forEach((it, i) => {
        const yTop = this.y(it.top.z);
        el("line", { class: "end-tick", x1: it.x, y1: yTop - 4, x2: pos[i], y2: y0 - 12 }, this.gLabels);
        el("text", { class: "end-label", x: pos[i], y: y0 - 15, "text-anchor": "middle" }, this.gLabels).textContent = it.pr.short;
      });
    }
    this.hideHover();
  }
  onHover(e) {
    const rect = this.svg.getBoundingClientRect();
    const sy = ((e.clientY - rect.top) / rect.height) * this.H;
    const z = Math.min(this.L, Math.max(0, ((this.H - this.m.b - sy) / (this.H - this.m.t - this.m.b)) * this.L));
    const yy = this.y(z);
    this.gHover.innerHTML = "";
    this.gHover.setAttribute("visibility", "visible");
    el("line", { class: "hover-row", x1: this.m.l, x2: this.W - this.m.r, y1: yy, y2: yy }, this.gHover);
    const rows = this.res.profiles.map((pr) => {
      const t = thetaAt(pr, z);
      el("circle", { class: `hover-pt dot-${pr.key}`, cx: this.x(t), cy: yy, r: 4.5 }, this.gHover);
      return `<tr><td><i class="sw" style="border-color:var(--${COLOR_VAR[pr.key]})"></i>${sub(pr.label)}</td><td>${t.toFixed(3)}</td></tr>`;
    });
    this.tooltip.hidden = false;
    this.tooltip.innerHTML = `z = ${z.toFixed(1)} cm<table>${rows.join("")}</table>`;
    const tx = (this.x(0.5) / this.W) * rect.width;
    this.tooltip.style.left = `${tx}px`;
    this.tooltip.style.top = `${(yy / this.H) * rect.height}px`;
  }
  hideHover() {
    this.gHover.setAttribute("visibility", "hidden");
    if (this.tooltip) this.tooltip.hidden = true;
  }
}
const COLOR_VAR = { "perc-half": "perc-1", "perc-quarter": "perc-2", hydro: "hydro", "evap-half": "evap-1", "evap-97": "evap-2" };

/* ───────────────────────── UI ───────────────────────── */
const $ = (sel) => document.querySelector(sel);
const form = $("#params");
const mainSvg = $("#main-svg");
const mainChart = new ProfileChart(mainSvg, mainSvg.parentElement.querySelector(".tooltip"));
let res = null;

// presets
const presetSel = $("#preset");
presetSel.innerHTML =
  `<option value="hw">HW soil (h_A = −30, λ = 1.2)</option>` +
  TEXTURES.map((t) => `<option value="${t.key}">${t.name}</option>`).join("") +
  `<option value="custom" hidden>Custom</option>`;

// small multiples
const multiples = MULTIPLES.map((key) => {
  const t = TEXTURES.find((x) => x.key === key);
  const wrap = document.createElement("div");
  wrap.className = "mini";
  wrap.innerHTML = `<p class="mini__title">${t.name} <small>|h<sub>A</sub>| ${t.hA * -1}, λ ${t.lam.toFixed(2)}</small></p>
    <div class="chart"><svg role="img" aria-label="Profiles for ${t.name}"></svg><div class="tooltip" hidden></div></div>`;
  $("#multiples").appendChild(wrap);
  const svg = wrap.querySelector("svg");
  return { t, chart: new ProfileChart(svg, wrap.querySelector(".tooltip"), { mini: true }) };
});

// legend
$("#key").innerHTML = CASES.map((cs) => `<li><i class="ln s-${cs.key}"></i><span data-key="${cs.key}"></span></li>`).join("");

function readInput() {
  const f = form.elements;
  const num = (name, lo, hi, dflt) => {
    const v = parseFloat(f[name].value);
    return Number.isFinite(v) ? Math.min(hi, Math.max(lo, v)) : dflt;
  };
  return {
    hA: num("hA", -500, -1, DEFAULTS.hA),
    lam: num("lam", 0.05, 5, DEFAULTS.lam),
    Ks: num("Ks", 1e-9, 10, DEFAULTS.Ks),
    L: num("L", 5, 1000, DEFAULTS.L),
  };
}
function writeInput(v) {
  const f = form.elements;
  f.hA.value = v.hA; f.lam.value = v.lam; f.Ks.value = v.Ks; f.L.value = v.L;
  syncOutputs();
}
function syncOutputs() {
  const f = form.elements;
  f.hA_val.value = fmt.minus(Number(f.hA.value).toFixed(1).replace(/\.0$/, ""));
  f.lam_val.value = Number(f.lam.value).toFixed(2);
  f.L_val.value = f.L.value;
}

function update() {
  const input = readInput();
  const ref = form.elements.ref.value;
  res = run(input, ref);
  const { p, kl, qm, profiles } = res;
  const by = Object.fromEntries(profiles.map((pr) => [pr.key, pr]));
  const Ks = input.Ks;
  const mmDay = (a) => a * Ks * 10 * SEC_PER_DAY;

  // readouts
  $("#c-out").textContent = p.c.toFixed(3);
  $("#eta-out").textContent = p.eta.toFixed(2);
  $("#qref-out").textContent = fmt.sig(kl * Ks, 3);
  $("#qref-note").textContent = `${fmt.sig(mmDay(kl), 3)} mm/day`;
  $("#qmax-out").textContent = fmt.sig(qm * Ks, 3);
  $("#qmax-note").textContent = `${(qm / kl).toFixed(2)} × K(L), ${fmt.sig(mmDay(qm), 3)} mm/day`;
  $("#tstar-out").textContent = Math.pow(0.5, 1 / p.c).toFixed(3);
  $("#tsurf-out").textContent = by["evap-97"].thetaSurf.toFixed(3);
  $("#tsurf-note").textContent = `for ${by["evap-97"].label}`;

  // legend
  for (const pr of profiles) {
    $(`#key [data-key="${pr.key}"]`).innerHTML = pr.key === "hydro" ? "0 (rest)" : sub(pr.label);
  }

  mainChart.draw(res);
  for (const mtp of multiples) {
    mtp.chart.draw(run({ hA: mtp.t.hA, lam: mtp.t.lam, Ks: mtp.t.Ks, L: input.L }, ref));
  }

  writeShapes(res, by, mmDay);
  writeForces(res, by);
  writeTexture(res, input, ref);
  writeTable(res);

  const pDry = profiles.find((pr) => pr.dryAt != null);
  $("#status").innerHTML = pDry
    ? `${sub(pDry.label)} is above q<sub>max</sub>: no steady profile reaches the surface. Θ<sub>E</sub> = 0 above z = ${pDry.dryAt.toFixed(1)} cm.`
    : `K<sub>s</sub> only sets the absolute fluxes. Every flux here is a multiple of K<sub>s</sub>, so the profiles do not depend on it.`;
}

function writeShapes({ p, kl, qm, ref }, by, mmDay) {
  const zs = ["perc-half", "perc-quarter", "hydro", "evap-half"].map((k) => by[k].zSat);
  $("#shape-fringe").innerHTML =
    `<strong>Every curve starts saturated.</strong> Suction stays below |h<sub>A</sub>| up to z = |h<sub>A</sub>| / (1 + q/K<sub>s</sub>):
    ${zs.slice(0, 3).map((z) => z.toFixed(0)).join(", ")} cm for −½ K<sub>s</sub>, −¼ K<sub>s</sub> and rest, and ${zs[3].toFixed(1)} cm for the upward cases.
    Below that the Brooks–Corey formula would give Θ<sub>E</sub> &gt; 1, so it is capped at 1.`;
  $("#shape-hydro").innerHTML =
    `<strong>Rest (q = 0).</strong> Suction equals height, so Θ<sub>E</sub> = (z/|h<sub>A</sub>|)<sup>−λ</sup>, falling to
    ${by.hydro.thetaSurf.toFixed(3)} at the surface. The other four curves sit on either side of this one.`;
  const pq = by["perc-quarter"], ph = by["perc-half"];
  $("#shape-perc").innerHTML = ph.saturated && pq.saturated
    ? `<strong>Percolation.</strong> The tension-saturated zone reaches the surface, so both downward cases stay at Θ<sub>E</sub> = 1.`
    : `<strong>Percolation is wetter than rest everywhere.</strong> Each curve turns vertical at the Θ<sub>E</sub> where K = |q|:
      Θ* = (|q|/K<sub>s</sub>)<sup>1/c</sup> = ${thetaE(p, p.sA * Math.pow(2, 1 / p.eta)).toFixed(3)} for −½ K<sub>s</sub> and
      ${thetaE(p, p.sA * Math.pow(4, 1 / p.eta)).toFixed(3)} for −¼ K<sub>s</sub>. Above that the column drains at unit gradient and
      Θ<sub>E</sub> barely changes (${ph.thetaSurf.toFixed(3)} and ${pq.thetaSurf.toFixed(3)} at the surface${ph.saturated ? "; −½ K<sub>s</sub> is still saturated there" : ""}).`;
  const e1 = by["evap-half"], e2 = by["evap-97"];
  $("#shape-evap").innerHTML =
    `<strong>Capillary rise is drier than rest, and the gap opens toward the surface.</strong>
    Θ<sub>E</sub> at the surface is ${e1.thetaSurf.toFixed(3)} for ${sub(e1.label)} and ${e2.thetaSurf.toFixed(3)} for ${sub(e2.label)},
    against ${by.hydro.thetaSurf.toFixed(3)} at rest. Low in the column K is still large and the curves follow the rest profile;
    near the surface K is small and the curves bend left.`;
  const expo = p.lam / (p.eta - 1);
  const at97qmax = profile(p, 0.97 * qm).thetaSurf;
  $("#shape-check").innerHTML = ref === "KL"
    ? `<strong>Check on the note "at .97, Θ<sub>E</sub> at the surface ≈ 0".</strong> Not with these numbers: it is ${e2.thetaSurf.toFixed(2)}.
      K(L) is not the most the soil can lift. The largest steady rise is q<sub>max</sub> = ${(qm / kl).toFixed(2)} K(L)
      (${fmt.sig(mmDay(qm), 3)} mm/day, Gardner 1958), where suction reaches ∞ exactly at the surface.
      Near q<sub>max</sub> the surface Θ<sub>E</sub> falls as (1 − q/q<sub>max</sub>)<sup>λ/(λc−1)</sup> = (1 − q/q<sub>max</sub>)<sup>${expo.toFixed(2)}</sup>,
      so even .97 q<sub>max</sub> leaves ${at97qmax.toFixed(2)}. Switch "Upward fluxes are fractions of" to q<sub>max</sub> to see it.`
    : `<strong>Measured against q<sub>max</sub> = ${(qm / kl).toFixed(2)} K(L).</strong> Θ<sub>E</sub> at the surface only reaches 0 at q<sub>max</sub> itself.
      Near it the surface value falls as (1 − q/q<sub>max</sub>)<sup>${expo.toFixed(2)}</sup>, so .97 q<sub>max</sub> still leaves ${at97qmax.toFixed(2)}.`;
}

function writeForces({ p }, by) {
  const ph = by["perc-half"];
  $("#force-perc").innerHTML =
    `<strong>Percolation (q &lt; 0).</strong> Water moves down, so drag points up and carries part of the weight.
    Suction only has to rise at 1 − |q|/K per cm, so it builds slowly and the column stays wet. That is why the saturated zone
    stretches to ${ph.zSat.toFixed(0)} cm for −½ K<sub>s</sub>. Higher up, K falls until it equals |q|. There drag carries all of the
    weight, the suction gradient is zero, and the water falls at a steady rate set by gravity and friction alone, like a
    raindrop at terminal velocity. Θ<sub>E</sub> stops changing, so the curve goes vertical.`;
  $("#force-evap").innerHTML =
    `<strong>Capillary rise (q &gt; 0).</strong> Water moves up, so drag points down and adds to gravity. Suction has to rise
    at 1 + q/K per cm to pull the water up against both. K falls as s<sup>−${p.eta.toFixed(1)}</sup>, so drag grows very fast as the soil dries:
    a little more suction means a lot less K, which needs more suction still. The runaway is what bends the curves toward
    Θ<sub>E</sub> = 0 near the surface. q<sub>max</sub> is the flux where the runaway reaches s = ∞ right at the surface. A larger
    evaporative demand cannot be met, so the surface dries and evaporation becomes soil-limited.`;
}

function writeTexture({ p }, input, ref) {
  const soils = [{ key: "hw", name: "HW soil", ...DEFAULTS }, ...TEXTURES];
  const rows = soils.map((t) => {
    const q = soil({ hA: t.hA, lam: t.lam, L: input.L });
    const qm = qMax(q);
    return { t, q, thetaStar: Math.pow(0.5, 1 / q.c), surf: Math.pow(Math.max(1, input.L / q.sA), -q.lam), qmax: qm * t.Ks * 10 * SEC_PER_DAY };
  });
  $("#tex-rows").innerHTML = rows.map((r) => `<tr data-key="${r.t.key}">
    <td class="l">${r.t.name}</td><td>${fmt.sig(r.q.sA, 3)}</td><td>${r.q.lam.toFixed(2)}</td><td>${r.q.c.toFixed(1)}</td>
    <td>${r.thetaStar.toFixed(3)}</td><td>${r.surf.toFixed(2)}</td><td>${fmt.sig(r.qmax, 2)}</td></tr>`).join("");

  const hw = rows[0];
  const ch = rows.slice(1);
  const lo = (f) => ch.reduce((a, b) => (f(b) < f(a) ? b : a));
  const hi = (f) => ch.reduce((a, b) => (f(b) > f(a) ? b : a));
  $("#tex-intro").innerHTML =
    `Three soil numbers set the curves. With every flux scaled by K<sub>s</sub>, K<sub>s</sub> itself drops out of Θ<sub>E</sub>(z).
    Water table at L = ${input.L} cm.`;
  $("#tex-points").innerHTML =
    `<strong>|h<sub>A</sub>| sets the saturated zone.</strong> Finer soils have smaller pores, a larger air-entry suction and a taller
    tension-saturated zone (${lo((r) => r.q.sA).t.name.toLowerCase()} ${fmt.sig(lo((r) => r.q.sA).q.sA, 2)} cm, ${hi((r) => r.q.sA).t.name.toLowerCase()} ${fmt.sig(hi((r) => r.q.sA).q.sA, 2)} cm).
    <br><strong>λ sets how fast Θ<sub>E</sub> falls above it.</strong> The HW soil (λ = ${hw.q.lam.toFixed(1)}) has nearly uniform pores,
    so it drains sharply. The texture classes have λ = ${lo((r) => r.q.lam).q.lam.toFixed(2)}–${hi((r) => r.q.lam).q.lam.toFixed(2)} and stay wet far
    above the water table.
    <br><strong>c sets percolation.</strong> Θ* = (|q|/K<sub>s</sub>)<sup>1/c</sup> depends on nothing else. A large c (fine soil)
    means K drops fast with Θ<sub>E</sub>, so a small drop in Θ<sub>E</sub> already brings K down to |q|: Θ* is
    ${hi((r) => r.thetaStar).thetaStar.toFixed(2)} in ${hi((r) => r.thetaStar).t.name.toLowerCase()} against ${hw.thetaStar.toFixed(2)} in the HW soil.
    <br><strong>Upward supply is a trade-off.</strong> Coarse soils have a large K<sub>s</sub> but lose it quickly as they dry. Fine soils start
    low but keep K and have a large |h<sub>A</sub>|. q<sub>max</sub> ranges from ${fmt.sig(lo((r) => r.qmax).qmax, 2)} mm/day
    (${lo((r) => r.qmax).t.name.toLowerCase()}) to ${fmt.sig(hi((r) => r.qmax).qmax, 2)} mm/day (${hi((r) => r.qmax).t.name.toLowerCase()}),
    and ${fmt.sig(hw.qmax, 2)} mm/day in the HW soil, whose K falls as s<sup>−${hw.q.eta.toFixed(1)}</sup>.`;
}

function writeTable({ p, profiles }) {
  $("#prof-head").innerHTML = `<tr><th scope="col">z <em>cm</em></th>${profiles.map((pr) => `<th scope="col">${sub(pr.label)}</th>`).join("")}</tr>`;
  const step = niceStep(p.L, 16);
  const zs = [];
  for (let z = 0; z <= p.L + 1e-9; z += step) zs.push(z);
  if (Math.abs(zs.at(-1) - p.L) > 1e-6) zs.push(p.L);
  $("#prof-rows").innerHTML = zs.reverse().map((z) =>
    `<tr><td>${Number(z.toFixed(1))}</td>${profiles.map((pr) => `<td>${thetaAt(pr, z).toFixed(3)}</td>`).join("")}</tr>`).join("");
}

/* ───────────────────────── wiring ───────────────────────── */
form.addEventListener("input", (e) => {
  if (e.target.name === "preset") return;
  if (["hA", "lam", "Ks"].includes(e.target.name)) presetSel.value = "custom";
  syncOutputs();
  update();
});
form.addEventListener("submit", (e) => e.preventDefault());

function loadPreset(key) {
  const cur = readInput();
  if (key === "hw") writeInput(DEFAULTS);
  else {
    const t = TEXTURES.find((x) => x.key === key);
    if (!t) return;
    writeInput({ hA: t.hA, lam: Number(t.lam.toFixed(2)), Ks: t.Ks, L: cur.L });
  }
  presetSel.value = key;
  update();
}
presetSel.addEventListener("change", () => loadPreset(presetSel.value));
$("#tex-rows").addEventListener("click", (e) => {
  const tr = e.target.closest("tr[data-key]");
  if (tr) loadPreset(tr.dataset.key);
});

document.querySelectorAll(".tab[data-tab]").forEach((tab) => {
  tab.addEventListener("click", () => {
    document.querySelectorAll(".tab[data-tab]").forEach((t) => {
      const on = t === tab;
      t.setAttribute("aria-pressed", on);
      document.getElementById(t.dataset.tab).hidden = !on;
    });
  });
});

$("#reset").addEventListener("click", () => {
  form.elements.ref.value = "KL";
  loadPreset("hw");
});
$("#csv").addEventListener("click", () => {
  const { p, profiles, input } = res;
  const rows = [["z_cm", ...profiles.map((pr) => `thetaE_${pr.key}`)]];
  for (let z = 0; z <= p.L + 1e-9; z += 1) rows.push([z, ...profiles.map((pr) => thetaAt(pr, z).toFixed(5))]);
  const meta = `# h_A_cm=${input.hA}, lambda=${input.lam}, Ks_cm_s=${input.Ks}, L_cm=${p.L}, c=${p.c}, ` +
    profiles.map((pr) => `q_${pr.key}_cm_s=${(pr.a * input.Ks).toPrecision(4)}`).join(", ") + "\n";
  const blob = new Blob([meta + rows.map((r) => r.join(",")).join("\n")], { type: "text/csv" });
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob);
  a.download = `profiles_L${p.L}.csv`;
  a.click();
  setTimeout(() => URL.revokeObjectURL(a.href), 1000);
});

/* ───────────────────────── boot ───────────────────────── */
writeInput(DEFAULTS);
update();
new ResizeObserver(() => { if (res) update(); }).observe(mainSvg);
