/* ───────────────────────── defaults ───────────────────────── */
const DEFAULTS = {
  n: 1_000_000,
  rmean: 0.01,     // mm
  cut: 3,          // largest pore, in multiples of r̄
  Rs: 20,          // mm, sample radius
  seed: 42,
  sigma: 0.072,    // N/m
  rho: 1000,       // kg/m3
  g: 9.8,          // m/s2
  mu: 0.001,       // kg/(m s)
  thetaMin: 0.01,
};
const GRID_DECADES = 3;
const GRID_PER_DECADE = 100;
const LOCAL_HALF_WIDTH = 5;
const THETA_FLOOR = 1e-4;
const RANGE_ROWS = [0.1, 0.01, 0.001, 0.0001];
const TABLE_ROWS = [1, 0.8, 0.6, 0.4, 0.2, 0.1, 0.05, 0.01, 0.001, 0.0001];

/* ───────────────────────── RNG ───────────────────────── */
function mulberry32(seed) {
  let a = seed >>> 0;
  return function () {
    a = (a + 0x6d2b79f5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

/* ───────────────────────── model ───────────────────────── */
function runModel(p) {
  const rng = mulberry32(p.seed);
  const rCap = p.cut * p.rmean;
  const drawn = new Float64Array(p.n);
  let N = 0;
  for (let i = 0; i < p.n; i++) {
    const x = -p.rmean * Math.log(1 - rng());
    if (x <= rCap) drawn[N++] = x;
  }
  const r = drawn.slice(0, N).sort();

  // cumulative r² and r⁴ (ascending)
  const cum2 = new Float64Array(N + 1);
  const cum4 = new Float64Array(N + 1);
  for (let i = 0; i < N; i++) {
    const r2 = r[i] * r[i];
    cum2[i + 1] = cum2[i] + r2;
    cum4[i + 1] = cum4[i] + r2 * r2;
  }

  const R2 = p.Rs * p.Rs;
  const phi = cum2[N] / R2;
  // h[mm] = -C / r[mm]  with  C = 2 σ · 1e6 / (ρ g)  [mm²]
  const C = (2 * p.sigma * 1e6) / (p.rho * p.g);
  const hA = C / r[N - 1];
  // K = (ρ g / 8 μ) Σ r⁴ / R², with mm⁴ / mm² → m²
  const kOf = (sum4) => ((p.rho * p.g) / (8 * p.mu)) * (sum4 / R2) * 1e-6;
  const Ksat = kOf(cum4[N]);

  function countBelow(rstar) {
    let lo = 0, hi = N;
    while (lo < hi) {
      const mid = (lo + hi) >> 1;
      if (r[mid] <= rstar) lo = mid + 1; else hi = mid;
    }
    return lo;
  }

  const grid = [];
  for (let i = 0; i <= GRID_DECADES * GRID_PER_DECADE; i++) {
    const h = hA * Math.pow(10, i / GRID_PER_DECADE);
    const rstar = C / h;
    const k = countBelow(rstar);
    grid.push({
      h, rstar, drained: N - k,
      thetaE: cum2[k] / cum2[N],
      theta: cum2[k] / R2,
      K: kOf(cum4[k]),
    });
  }
  grid[0].thetaE = 1;
  grid[0].theta = phi;
  grid[0].K = Ksat;
  const curve = grid.filter((g) => g.thetaE > 0);

  // exponents from a centered slope along the log-spaced curve
  const local = [];
  for (let i = 0; i < curve.length; i++) {
    const a = curve[Math.max(0, i - LOCAL_HALF_WIDTH)];
    const b = curve[Math.min(curve.length - 1, i + LOCAL_HALF_WIDTH)];
    const dT = Math.log(a.thetaE / b.thetaE);
    if (!(dT > 0) || curve[i].thetaE < THETA_FLOOR) continue;
    const lam = dT / Math.log(b.h / a.h);
    local.push({
      thetaE: curve[i].thetaE,
      lam,
      c: Math.log(a.K / b.K) / dT,
      predicted: 1 + 2 / lam,
    });
  }

  // closed forms for the same truncated exponential, as a check on the sampling
  const t = p.cut, e = Math.exp(-t);
  const p3 = 1 - e * (1 + t + t * t / 2);
  const p5 = 1 - e * (1 + t + t * t / 2 + t ** 3 / 6 + t ** 4 / 24);
  const phiExact = (p.n * p.rmean ** 2 * 2 * p3) / R2;
  const KsatExact = kOf(p.n * p.rmean ** 4 * 24 * p5);

  return { p, N, r, phi, C, hA, Ksat, curve, local, phiExact, KsatExact };
}

// least squares through (Θ_E = 1, h_A) and (Θ_E = 1, K_sat), both fixed by the model form
function fit(model, thetaMin) {
  const pts = model.curve.filter((g) => g.thetaE >= thetaMin);
  let tt = 0, th = 0, tk = 0;
  for (const g of pts) {
    const t = Math.log(g.thetaE);
    tt += t * t;
    th += t * Math.log(g.h / model.hA);
    tk += t * Math.log(g.K / model.Ksat);
  }
  return { lam: -tt / th, c: tk / tt, pts };
}

// RMS misfit in decades over the fit range
function misfit(model, pts, lam, c) {
  let eh = 0, ek = 0;
  for (const g of pts) {
    const t = Math.log(g.thetaE);
    eh += (Math.log(g.h / model.hA) + t / lam) ** 2;
    ek += (Math.log(g.K / model.Ksat) - c * t) ** 2;
  }
  return { h: Math.sqrt(eh / pts.length) / Math.LN10, K: Math.sqrt(ek / pts.length) / Math.LN10 };
}

/* ───────────────────────── formatting ───────────────────────── */
const SUP = { "-": "⁻", 0: "⁰", 1: "¹", 2: "²", 3: "³", 4: "⁴", 5: "⁵", 6: "⁶", 7: "⁷", 8: "⁸", 9: "⁹" };
const THETA_E = `Θ<tspan dy="3" font-size="10">E</tspan><tspan dy="-3"> </tspan>`;
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
    if (a >= 1e5 || a < 1e-3) return fmt.sci(x, n);
    return Number(x.toPrecision(n)).toString();
  },
  pow(d) { return d >= 0 && d <= 4 ? Math.pow(10, d).toLocaleString() : d < 0 && d >= -3 ? Math.pow(10, d).toFixed(-d) : fmt.sci(Math.pow(10, d), 1); },
};

/* ───────────────────────── chart (SVG) ───────────────────────── */
const svgNS = "http://www.w3.org/2000/svg";
function el(name, attrs = {}, parent) {
  const e = document.createElementNS(svgNS, name);
  for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v);
  if (parent) parent.appendChild(e);
  return e;
}

function logTicks(lo, hi) {
  const dLo = Math.floor(Math.log10(lo)), dHi = Math.ceil(Math.log10(hi));
  const steps = Math.log10(hi / lo) < 2.5 ? [1, 2, 5] : [1];
  const out = [];
  for (let d = dLo; d <= dHi; d++) {
    for (const s of steps) {
      const v = s * Math.pow(10, d);
      if (v >= lo * 0.999 && v <= hi * 1.001) out.push(v);
    }
  }
  return out;
}

class Chart {
  constructor(svg, tooltip, { xLabel, yLabel, yLog = true, tip }) {
    this.svg = svg;
    this.tooltip = tooltip;
    this.yLog = yLog;
    this.tip = tip;
    this.m = { t: 8, r: 14, b: 40, l: 66 };
    this.gBand = el("g", {}, svg);
    this.gGrid = el("g", { class: "grid" }, svg);
    this.gSeries = el("g", {}, svg);
    this.gAxis = el("g", { class: "axis" }, svg);
    this.gTicks = el("g", { class: "tick" }, svg);
    this.cross = el("line", { class: "crosshair", visibility: "hidden" }, svg);
    this.dot = el("circle", { class: "hover-dot", r: 5, visibility: "hidden" }, svg);
    this.hit = el("rect", { class: "hit" }, svg);
    this.xText = el("text", { class: "axis-label", "text-anchor": "middle" }, svg);
    this.xText.innerHTML = xLabel;
    this.yText = el("text", { class: "axis-label", "text-anchor": "middle" }, svg);
    this.yText.textContent = yLabel;
    this.hit.addEventListener("pointermove", (e) => this.onHover(e));
    this.hit.addEventListener("pointerleave", () => this.hideHover());
  }
  x(v) {
    const { m, W } = this;
    return m.l + (Math.log10(v / this.xMin) / Math.log10(this.xMax / this.xMin)) * (W - m.l - m.r);
  }
  y(v) {
    const { m, H } = this;
    const hgt = H - m.t - m.b;
    const frac = this.yLog
      ? Math.log10(v / this.yMin) / Math.log10(this.yMax / this.yMin)
      : (v - this.yMin) / (this.yMax - this.yMin);
    return m.t + hgt - frac * hgt;
  }
  // series: [{ cls, pts: [[x, y], …] }], hover: [{ x, y, row }]
  draw({ xMin, xMax, yMin, yMax, series, hover, excludeBelow }) {
    Object.assign(this, { xMin, xMax, yMin, yMax, hover });
    // the svg is sized by its panel, so the drawing space follows it
    this.W = this.svg.clientWidth || 640;
    this.H = this.svg.clientHeight || 420;
    const { m, W, H } = this;
    const x0 = m.l, x1 = W - m.r, y0 = m.t, y1 = H - m.b;
    this.svg.setAttribute("viewBox", `0 0 ${W} ${H}`);
    for (const [k, v] of Object.entries({ x: x0, y: y0, width: x1 - x0, height: y1 - y0 })) this.hit.setAttribute(k, v);
    this.xText.setAttribute("x", (x0 + x1) / 2);
    this.xText.setAttribute("y", H - 6);
    this.yText.setAttribute("transform", `translate(13 ${(y0 + y1) / 2}) rotate(-90)`);
    for (const g of [this.gGrid, this.gTicks, this.gAxis, this.gSeries, this.gBand]) g.innerHTML = "";

    if (excludeBelow && excludeBelow > xMin) {
      const xx = this.x(excludeBelow);
      el("rect", { class: "excluded", x: x0, y: y0, width: xx - x0, height: y1 - y0 }, this.gBand);
      el("text", { class: "excluded-label", x: xx - 6, y: y0 + 16, "text-anchor": "end" }, this.gBand).textContent = "outside fit";
    }
    for (const v of logTicks(xMin, xMax)) {
      const xx = this.x(v);
      el("line", { x1: xx, x2: xx, y1: y0, y2: y1 }, this.gGrid);
      el("text", { x: xx, y: y1 + 18, "text-anchor": "middle" }, this.gTicks).textContent = fmt.pow(Math.round(Math.log10(v)));
    }
    const yTicks = this.yLog ? logTicks(yMin, yMax) : linTicks(yMin, yMax);
    for (const v of yTicks) {
      const yy = this.y(v);
      el("line", { x1: x0, x2: x1, y1: yy, y2: yy }, this.gGrid);
      el("text", { x: x0 - 8, y: yy + 4, "text-anchor": "end" }, this.gTicks).textContent = this.yLog ? fmt.sig(v, 1) : v.toFixed(1);
    }
    el("path", { d: `M${x0} ${y0} V${y1} H${x1}` }, this.gAxis);

    const inside = (a) => a[0] >= xMin && a[0] <= xMax * 1.0001 && a[1] >= yMin && a[1] <= yMax;
    for (const s of series) {
      const d = s.pts.filter(inside).map((a) => `${this.x(a[0]).toFixed(2)} ${this.y(a[1]).toFixed(2)}`);
      if (d.length > 1) el("path", { class: `series ${s.cls || ""}`, d: "M" + d.join(" L") }, this.gSeries);
    }
    this.hideHover();
  }
  onHover(e) {
    if (!this.hover?.length) return;
    const rect = this.svg.getBoundingClientRect();
    const sx = ((e.clientX - rect.left) / rect.width) * this.W;
    let best = this.hover[0], bestD = Infinity;
    for (const a of this.hover) {
      const d = Math.abs(this.x(a.x) - sx);
      if (d < bestD) { bestD = d; best = a; }
    }
    const xx = this.x(best.x), yy = this.y(best.y);
    for (const c of [this.cross, this.dot]) c.setAttribute("visibility", "visible");
    this.cross.setAttribute("x1", xx); this.cross.setAttribute("x2", xx);
    this.cross.setAttribute("y1", this.m.t); this.cross.setAttribute("y2", this.H - this.m.b);
    this.dot.setAttribute("cx", xx); this.dot.setAttribute("cy", yy);
    this.tooltip.hidden = false;
    this.tooltip.innerHTML = this.tip(best.row);
    this.tooltip.style.left = `${(xx / this.W) * rect.width}px`;
    this.tooltip.style.top = `${(yy / this.H) * rect.height}px`;
  }
  hideHover() {
    for (const c of [this.cross, this.dot]) c.setAttribute("visibility", "hidden");
    this.tooltip.hidden = true;
  }
}
function linTicks(lo, hi) {
  const out = [];
  for (let v = Math.ceil(lo * 2) / 2; v <= hi + 1e-9; v += 0.5) out.push(v);
  return out;
}
function makeChart(id, opts) {
  const svg = document.getElementById(id);
  return new Chart(svg, svg.parentElement.querySelector(".tooltip"), opts);
}

/* ───────────────────────── UI wiring ───────────────────────── */
const $ = (sel) => document.querySelector(sel);
const form = $("#params");
let model = null;
let best = null;

const swrcChart = makeChart("swrc-svg", {
  xLabel: `saturation ${THETA_E}`,
  yLabel: "suction |h| (mm)",
  tip: (g) => `Θ<sub>E</sub> = ${fmt.sig(g.thetaE, 3)}<br>|h| = ${fmt.sig(g.h, 4)} mm<br>r* = ${fmt.sig(g.rstar, 3)} mm`,
});
const kChart = makeChart("k-svg", {
  xLabel: `saturation ${THETA_E}`,
  yLabel: "K (m/s)",
  tip: (g) => `Θ<sub>E</sub> = ${fmt.sig(g.thetaE, 3)}<br>K = ${fmt.sig(g.K, 3)} m/s<br>K/K<sub>sat</sub> = ${fmt.sig(g.K / model.Ksat, 3)}`,
});
const localChart = makeChart("local-svg", {
  xLabel: `saturation ${THETA_E}`,
  yLabel: "exponent",
  yLog: false,
  tip: (a) => `Θ<sub>E</sub> = ${fmt.sig(a.thetaE, 3)}<br>local λ = ${fmt.sig(a.lam, 3)}<br>local c = ${fmt.sig(a.c, 3)}<br>1 + 2/λ = ${fmt.sig(a.predicted, 3)}`,
});

function readParams() {
  const f = form.elements;
  const num = (name, lo, hi) => {
    const v = parseFloat(f[name].value);
    return Number.isFinite(v) ? Math.min(hi, Math.max(lo, v)) : DEFAULTS[name];
  };
  return {
    n: Math.round(num("n", 1000, 2_000_000)),
    rmean: num("rmean", 1e-4, 1),
    cut: num("cut", 1, 20),
    Rs: num("Rs", 1, 1000),
    seed: Math.round(num("seed", 0, 999999)),
    sigma: num("sigma", 1e-3, 1),
    rho: num("rho", 1, 20000),
    g: num("g", 0.01, 100),
    mu: num("mu", 1e-5, 1),
    thetaMin: parseFloat(f.thetaMin.value),
  };
}

function writeParams(p) {
  const f = form.elements;
  for (const name of ["n", "rmean", "cut", "Rs", "seed", "sigma", "rho", "g", "mu"]) f[name].value = p[name];
  f.thetaMin.value = String(p.thetaMin);
}

function setExponents(lam, c) {
  const f = form.elements;
  f.lam.value = lam; f.lam_val.value = lam.toFixed(2);
  f.c.value = c; f.c_val.value = c.toFixed(2);
}

function currentExponents() {
  const f = form.elements;
  return { lam: f.lam.valueAsNumber, c: f.c.valueAsNumber };
}

function drawCurves() {
  const { lam, c } = currentExponents();
  const thetaMin = readParams().thetaMin;
  const { curve, hA, Ksat, local } = model;
  const shown = curve.filter((g) => g.thetaE >= THETA_FLOOR);
  const hTop = shown.at(-1).h * 1.05;
  const kBottom = shown.at(-1).K * 0.8;

  swrcChart.draw({
    xMin: THETA_FLOOR, xMax: 1, yMin: hA * 0.8, yMax: hTop,
    excludeBelow: thetaMin,
    series: [
      { pts: curve.map((g) => [g.thetaE, g.h]) },
      { cls: "series--fit", pts: shown.map((g) => [g.thetaE, hA * Math.pow(g.thetaE, -1 / lam)]) },
    ],
    hover: shown.map((g) => ({ x: g.thetaE, y: g.h, row: g })),
  });
  kChart.draw({
    xMin: THETA_FLOOR, xMax: 1, yMin: kBottom, yMax: Ksat * 1.6,
    excludeBelow: thetaMin,
    series: [
      { pts: curve.map((g) => [g.thetaE, g.K]) },
      { cls: "series--fit", pts: shown.map((g) => [g.thetaE, Ksat * Math.pow(g.thetaE, c)]) },
      { cls: "series--pred", pts: shown.map((g) => [g.thetaE, Ksat * Math.pow(g.thetaE, 1 + 2 / lam)]) },
    ],
    hover: shown.map((g) => ({ x: g.thetaE, y: g.K, row: g })),
  });
  localChart.draw({
    xMin: THETA_FLOOR, xMax: 1, yMin: 1.5,
    yMax: Math.max(3, Math.ceil(2 * Math.max(...local.map((a) => a.predicted))) / 2),
    series: [
      { pts: local.map((a) => [a.thetaE, a.c]) },
      { cls: "series--pred", pts: local.map((a) => [a.thetaE, a.predicted]) },
    ],
    hover: local.map((a) => ({ x: a.thetaE, y: a.c, row: a })),
  });

  const pred = 1 + 2 / lam;
  const m = misfit(model, best.pts, lam, c);
  const gap = (100 * (pred - c)) / c;
  $("#lam-out").textContent = lam.toFixed(2);
  $("#c-out").textContent = c.toFixed(2);
  $("#c-bc").textContent = pred.toFixed(2);
  $("#lam-rms").textContent = `off by ±${m.h.toFixed(3)} decades`;
  $("#c-rms").textContent = `off by ±${m.K.toFixed(3)} decades`;
  $("#c-gap").textContent = `${gap >= 0 ? "+" : ""}${gap.toFixed(1)}% vs c`;

  const bestPred = 1 + 2 / best.lam;
  const bestGap = (100 * (bestPred - best.c)) / best.c;
  $("#answer").innerHTML =
    `<strong>Close. Fitting down to Θ<sub>E</sub> = ${thetaMin} gives λ = ${best.lam.toFixed(2)} and c = ${best.c.toFixed(2)}.
    Theory says ${bestPred.toFixed(2)}, which is ${Math.abs(bestGap).toFixed(0)}% ${bestGap >= 0 ? "high" : "low"}.</strong>`;
}

function drawTables() {
  const { curve, Ksat } = model;
  const nearest = (t) => curve.reduce((a, b) =>
    Math.abs(Math.log(b.thetaE / t)) < Math.abs(Math.log(a.thetaE / t)) ? b : a);
  $("#rows").innerHTML = TABLE_ROWS.map(nearest).map((g) => `<tr>
    <td>${fmt.sig(g.h, 4)}</td><td>${fmt.sig(g.rstar, 3)}</td><td>${fmt.sig(g.thetaE, 3)}</td>
    <td>${fmt.sig(g.theta, 3)}</td><td>${fmt.sig(g.K, 3)}</td><td>${fmt.sig(g.K / Ksat, 3)}</td>
  </tr>`).join("");

  const active = readParams().thetaMin;
  $("#ranges").innerHTML = RANGE_ROWS.map((t) => {
    const f = fit(model, t);
    const pred = 1 + 2 / f.lam;
    const gap = (100 * (pred - f.c)) / f.c;
    return `<tr class="${t === active ? "is-active" : ""}"><td>${t}</td><td>${f.lam.toFixed(2)}</td><td>${f.c.toFixed(2)}</td>
      <td>${pred.toFixed(2)}</td><td>${gap >= 0 ? "+" : ""}${gap.toFixed(1)}%</td></tr>`;
  }).join("");
}

function refit() {
  best = fit(model, readParams().thetaMin);
  setExponents(best.lam, best.c);
  drawCurves();
  drawTables();
}

function run() {
  const t0 = performance.now();
  model = runModel(readParams());
  const ms = performance.now() - t0;
  $("#phi").textContent = fmt.sig(model.phi, 3);
  $("#kept").textContent = `${model.N.toLocaleString()} pores kept`;
  $("#ha").innerHTML = `${fmt.sig(model.hA, 3)} <small>mm</small>`;
  $("#ksat").textContent = fmt.sig(model.Ksat, 3);
  $("#status").innerHTML =
    `Sampled / closed form: φ ${fmt.sig(model.phi, 4)} / ${fmt.sig(model.phiExact, 4)}, K<sub>sat</sub> ${fmt.sig(model.Ksat, 4)} / ${fmt.sig(model.KsatExact, 4)}. ${ms.toFixed(0)} ms.`;
  refit();
}

let timer = null;
form.addEventListener("input", (e) => {
  const name = e.target.name;
  if (name === "lam" || name === "c") {
    const f = form.elements;
    f.lam_val.value = f.lam.valueAsNumber.toFixed(2);
    f.c_val.value = f.c.valueAsNumber.toFixed(2);
    drawCurves();
  } else if (name === "thetaMin") {
    refit();
  } else {
    clearTimeout(timer);
    timer = setTimeout(run, 200);
  }
});
form.addEventListener("submit", (e) => e.preventDefault());

document.querySelectorAll(".tab[data-tab]").forEach((tab) => {
  tab.addEventListener("click", () => {
    document.querySelectorAll(".tab[data-tab]").forEach((t) => {
      const on = t === tab;
      t.setAttribute("aria-pressed", on);
      document.getElementById(t.dataset.tab).hidden = !on;
    });
  });
});

$("#best").addEventListener("click", refit);
$("#reseed").addEventListener("click", () => {
  form.elements.seed.value = Math.floor(Math.random() * 1e6);
  run();
});
$("#reset").addEventListener("click", () => {
  writeParams(DEFAULTS);
  run();
});
$("#csv").addEventListener("click", () => {
  const rows = [["h_mm", "r_star_mm", "theta", "theta_E", "K_m_per_s", "K_over_Ksat"]];
  for (const g of model.curve) {
    rows.push([-g.h, g.rstar, g.theta, g.thetaE, g.K, g.K / model.Ksat].map((v) => Number(v.toPrecision(6))));
  }
  const p = model.p;
  const meta = `# porosity=${model.phi}, Ksat_m_per_s=${model.Ksat}, h_A_mm=${model.hA}, N_drilled=${p.n}, N_kept=${model.N}, rmean_mm=${p.rmean}, cutoff_x_rmean=${p.cut}, R_mm=${p.Rs}, seed=${p.seed}\n`;
  const blob = new Blob([meta + rows.map((r) => r.join(",")).join("\n")], { type: "text/csv" });
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob);
  a.download = `conductivity_N${p.n}_seed${p.seed}.csv`;
  a.click();
  setTimeout(() => URL.revokeObjectURL(a.href), 1000);
});

/* ───────────────────────── boot ───────────────────────── */
writeParams(DEFAULTS);
run();
new ResizeObserver(() => { if (model) drawCurves(); }).observe(document.getElementById("swrc-svg"));
