import * as THREE from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";

/* ───────────────────────── defaults ───────────────────────── */
const DEFAULTS = {
  dist: "exponential",
  rmean: 0.01,     // mm
  shape: 2,        // gamma k
  sigln: 0.8,      // lognormal sigma
  cutoff: true,
  rmax: 0.2,       // mm
  n: 1_000_000,
  seed: 42,
  D: 20,           // mm
  L: 20,           // mm
  sigma: 0.072,    // N/m
  gamma: 0,        // deg
  rho: 1000,       // kg/m3
  g: 9.8,          // m/s2
};
const MAX_DRAWN_PORES = 700;
const CURVE_POINTS = 240;

/* ───────────────────────── RNG & samplers ───────────────────────── */
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
function makeNormal(rng) {
  let spare = null;
  return function () {
    if (spare !== null) { const s = spare; spare = null; return s; }
    let u, v, s;
    do { u = rng() * 2 - 1; v = rng() * 2 - 1; s = u * u + v * v; } while (s >= 1 || s === 0);
    const m = Math.sqrt((-2 * Math.log(s)) / s);
    spare = v * m;
    return u * m;
  };
}
function makeGamma(rng, normal, k) {
  // Marsaglia & Tsang; boost for k < 1
  const boost = k < 1;
  const kk = boost ? k + 1 : k;
  const d = kk - 1 / 3;
  const c = 1 / Math.sqrt(9 * d);
  return function () {
    let x, v, u;
    for (;;) {
      do { x = normal(); v = 1 + c * x; } while (v <= 0);
      v = v * v * v;
      u = rng();
      if (u < 1 - 0.0331 * x * x * x * x) break;
      if (Math.log(u) < 0.5 * x * x + d * (1 - v + Math.log(v))) break;
    }
    const gkk = d * v;
    return boost ? gkk * Math.pow(rng(), 1 / k) : gkk;
  };
}
function makeSampler(p, rng) {
  const normal = makeNormal(rng);
  if (p.dist === "gamma") {
    const g = makeGamma(rng, normal, p.shape);
    const scale = p.rmean / p.shape;
    return () => g() * scale;
  }
  if (p.dist === "lognormal") {
    const mu = Math.log(p.rmean) - 0.5 * p.sigln * p.sigln;
    return () => Math.exp(mu + p.sigln * normal());
  }
  return () => -p.rmean * Math.log(1 - rng());
}

/* ───────────────────────── model ───────────────────────── */
function capillaryConstant(p) {
  // h[mm] = -C / r[mm]  with  C = 2 σ cosγ · 1e6 / (ρ g)  [mm²]
  return (2 * p.sigma * Math.cos((p.gamma * Math.PI) / 180) * 1e6) / (p.rho * p.g);
}

function runModel(p) {
  const rng = mulberry32(p.seed);
  const sample = makeSampler(p, rng);
  const N = p.n;
  const r = new Float64Array(N);
  let clamped = 0;
  for (let i = 0; i < N; i++) {
    let x = sample();
    if (p.cutoff) {
      let tries = 0;
      while (x > p.rmax && tries < 60) { x = sample(); tries++; }
      if (x > p.rmax) { x = p.rmax; clamped++; }
    }
    r[i] = x;
  }
  r.sort();

  // cumulative r² (ascending)
  const cum = new Float64Array(N + 1);
  let s = 0;
  for (let i = 0; i < N; i++) { s += r[i] * r[i]; cum[i + 1] = s; }

  const R = p.D / 2;
  const phi = s / (R * R);
  const C = capillaryConstant(p);

  const rMax = r[N - 1];
  const rMin = r[0];
  const rMedian = r[N >> 1];
  const realizedMean = r.reduce((a, b) => a + b, 0) / N;

  // head range: air entry of the largest pore → head at which only 0.2% of the pore volume remains
  const hAirEntry = C / rMax;
  let iq = 0;
  while (iq < N && cum[iq + 1] < 0.002 * s) iq++;
  const rQ = r[Math.min(N - 1, iq)];
  const hLo = Math.max(1e-3, hAirEntry / 4);
  const hHi = Math.min(1e12, (C / Math.max(rQ, 1e-9)) * 1.3);

  // SWRC curve on a log grid in |h|
  const curve = [];
  const lgLo = Math.log10(hLo), lgHi = Math.log10(hHi);
  for (let i = 0; i < CURVE_POINTS; i++) {
    const h = Math.pow(10, lgLo + ((lgHi - lgLo) * i) / (CURVE_POINTS - 1));
    curve.push({ h, theta: thetaAt(h) });
  }

  function countBelow(rstar) {
    // number of pores with r <= rstar (binary search on ascending r)
    let lo = 0, hi = N;
    while (lo < hi) {
      const mid = (lo + hi) >> 1;
      if (r[mid] <= rstar) lo = mid + 1; else hi = mid;
    }
    return lo;
  }
  function thetaAt(h) {
    if (h <= 0) return phi;
    const rstar = C / h;
    return cum[countBelow(rstar)] / (R * R);
  }
  function stateAt(h) {
    const rstar = h <= 0 ? Infinity : C / h;
    const kept = countBelow(rstar);
    const theta = cum[kept] / (R * R);
    return { h, rstar, theta, sat: phi > 0 ? theta / phi : 0, drained: N - kept };
  }

  return {
    p, r, cum, N, R, phi, C, rMax, rMin, rMedian, realizedMean, clamped,
    hAirEntry, hLo, hHi, curve, thetaAt, stateAt,
  };
}

/* ───────────────────────── formatting ───────────────────────── */
const fmt = {
  sig(x, n = 3) {
    if (!isFinite(x)) return "–";
    if (x === 0) return "0";
    const a = Math.abs(x);
    if (a >= 1e5 || a < 1e-3) return x.toExponential(n - 1).replace("e+", "e");
    return Number(x.toPrecision(n)).toString();
  },
  head(h) {
    if (h === 0) return "0";
    if (h >= 1e5) return "−" + h.toExponential(2).replace("e+", "e");
    if (h >= 100) return "−" + Math.round(h).toLocaleString();
    return "−" + fmt.sig(h, 3);
  },
  pct(x) { return (100 * x).toFixed(1) + "%"; },
  int(x) { return Math.round(x).toLocaleString(); },
};

/* ───────────────────────── chart (SVG) ───────────────────────── */
const svgNS = "http://www.w3.org/2000/svg";
function el(name, attrs = {}, parent) {
  const e = document.createElementNS(svgNS, name);
  for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v);
  if (parent) parent.appendChild(e);
  return e;
}

class Chart {
  constructor(svg, tooltip) {
    this.svg = svg;
    this.tooltip = tooltip;
    this.W = 640; this.H = 420;
    this.m = { t: 18, r: 22, b: 48, l: 64 };
    this.model = null;
    this.build();
  }
  build() {
    const s = this.svg;
    s.innerHTML = "";
    this.gGrid = el("g", { class: "grid" }, s);
    this.gArea = el("path", { class: "area" }, s);
    this.gSeries = el("path", { class: "series" }, s);
    this.gAxis = el("g", { class: "axis" }, s);
    this.gTicks = el("g", { class: "tick" }, s);
    this.aeLabel = el("text", { class: "ae-label", "text-anchor": "end" }, s);
    this.markerH = el("line", { class: "marker-line" }, s);
    this.markerV = el("line", { class: "marker-line" }, s);
    this.marker = el("circle", { class: "marker", r: 6 }, s);
    this.crossX = el("line", { class: "crosshair", visibility: "hidden" }, s);
    this.crossY = el("line", { class: "crosshair", visibility: "hidden" }, s);
    this.hoverDot = el("circle", { class: "hover-dot", r: 5, visibility: "hidden" }, s);
    this.hit = el("rect", { class: "hit" }, s);
    el("text", { class: "axis-label", x: this.m.l + (this.W - this.m.l - this.m.r) / 2, y: this.H - 10, "text-anchor": "middle" }, s)
      .textContent = "volumetric water content θ";
    const yl = el("text", { class: "axis-label", "text-anchor": "middle",
      transform: `translate(16 ${this.m.t + (this.H - this.m.t - this.m.b) / 2}) rotate(-90)` }, s);
    yl.textContent = "matric suction |h|  (mm, log)";

    this.hit.addEventListener("pointermove", (e) => this.onHover(e));
    this.hit.addEventListener("pointerleave", () => this.hideHover());
  }
  x(theta) { const w = this.W - this.m.l - this.m.r; return this.m.l + (theta / this.xMax) * w; }
  y(h) {
    const hgt = this.H - this.m.t - this.m.b;
    const lg = Math.log10(h), lo = Math.log10(this.hLo), hi = Math.log10(this.hHi);
    return this.m.t + hgt - ((lg - lo) / (hi - lo)) * hgt;
  }
  setModel(model) {
    this.model = model;
    this.hLo = model.hLo; this.hHi = model.hHi;
    this.xMax = niceMax(model.phi);
    const { m, W, H } = this;
    const x0 = m.l, x1 = W - m.r, y0 = m.t, y1 = H - m.b;

    // grid + ticks
    this.gGrid.innerHTML = ""; this.gTicks.innerHTML = ""; this.gAxis.innerHTML = "";
    const dLo = Math.ceil(Math.log10(this.hLo)), dHi = Math.floor(Math.log10(this.hHi));
    for (let d = dLo; d <= dHi; d++) {
      const yy = this.y(Math.pow(10, d));
      el("line", { x1: x0, x2: x1, y1: yy, y2: yy }, this.gGrid);
      const t = el("text", { x: x0 - 8, y: yy + 4, "text-anchor": "end" }, this.gTicks);
      t.textContent = d >= 0 && d <= 4 ? Math.pow(10, d).toLocaleString() : `1e${d}`;
    }
    const xStep = niceStep(this.xMax);
    for (let v = 0; v <= this.xMax + 1e-9; v += xStep) {
      const xx = this.x(v);
      el("line", { x1: xx, x2: xx, y1: y0, y2: y1 }, this.gGrid);
      const t = el("text", { x: xx, y: y1 + 18, "text-anchor": "middle" }, this.gTicks);
      t.textContent = Number(v.toFixed(3)).toString();
    }
    el("path", { d: `M${x0} ${y0} V${y1} H${x1}` }, this.gAxis);
    this.hit.setAttribute("x", x0); this.hit.setAttribute("y", y0);
    this.hit.setAttribute("width", x1 - x0); this.hit.setAttribute("height", y1 - y0);

    // series
    const pts = model.curve.map((c) => `${this.x(c.theta).toFixed(2)} ${this.y(c.h).toFixed(2)}`);
    this.gSeries.setAttribute("d", "M" + pts.join(" L"));
    this.gArea.setAttribute("d", `M${x0} ${this.y(model.curve[0].h).toFixed(2)} L` + pts.join(" L") + ` L${x0} ${this.y(model.curve.at(-1).h).toFixed(2)} Z`);

    // air-entry annotation
    const yae = this.y(model.hAirEntry);
    this.aeLabel.setAttribute("x", x1 - 4);
    this.aeLabel.setAttribute("y", yae - 6);
    this.aeLabel.textContent = `air entry ${fmt.sig(model.hAirEntry)} mm`;
  }
  setMarker(state) {
    if (!this.model) return;
    const h = Math.max(state.h, this.hLo);
    const xx = this.x(state.theta), yy = this.y(h);
    this.marker.setAttribute("cx", xx); this.marker.setAttribute("cy", yy);
    this.markerH.setAttribute("x1", this.m.l); this.markerH.setAttribute("x2", xx);
    this.markerH.setAttribute("y1", yy); this.markerH.setAttribute("y2", yy);
    this.markerV.setAttribute("x1", xx); this.markerV.setAttribute("x2", xx);
    this.markerV.setAttribute("y1", yy); this.markerV.setAttribute("y2", this.H - this.m.b);
  }
  onHover(e) {
    if (!this.model) return;
    const rect = this.svg.getBoundingClientRect();
    const sy = ((e.clientY - rect.top) / rect.height) * this.H;
    const hgt = this.H - this.m.t - this.m.b;
    const frac = 1 - (sy - this.m.t) / hgt;
    const lg = Math.log10(this.hLo) + frac * (Math.log10(this.hHi) - Math.log10(this.hLo));
    const h = Math.pow(10, lg);
    const st = this.model.stateAt(h);
    const xx = this.x(st.theta), yy = this.y(h);
    for (const c of [this.crossX, this.crossY, this.hoverDot]) c.setAttribute("visibility", "visible");
    this.crossX.setAttribute("x1", this.m.l); this.crossX.setAttribute("x2", this.W - this.m.r);
    this.crossX.setAttribute("y1", yy); this.crossX.setAttribute("y2", yy);
    this.crossY.setAttribute("x1", xx); this.crossY.setAttribute("x2", xx);
    this.crossY.setAttribute("y1", this.m.t); this.crossY.setAttribute("y2", this.H - this.m.b);
    this.hoverDot.setAttribute("cx", xx); this.hoverDot.setAttribute("cy", yy);
    const tt = this.tooltip;
    tt.hidden = false;
    tt.innerHTML = `h = ${fmt.head(h)} mm<br>θ = ${fmt.sig(st.theta, 4)}<br>S = ${fmt.pct(st.sat)} · r* = ${fmt.sig(st.rstar)} mm`;
    tt.style.left = `${(xx / this.W) * rect.width}px`;
    tt.style.top = `${(yy / this.H) * rect.height}px`;
  }
  hideHover() {
    for (const c of [this.crossX, this.crossY, this.hoverDot]) c.setAttribute("visibility", "hidden");
    this.tooltip.hidden = true;
  }
}
function niceStep(max) {
  const raw = max / 5;
  const p = Math.pow(10, Math.floor(Math.log10(raw)));
  const m = raw / p;
  return (m < 1.5 ? 1 : m < 3.5 ? 2 : m < 7.5 ? 5 : 10) * p;
}
function niceMax(v) {
  if (v <= 0) return 1;
  const step = niceStep(v);
  return Math.ceil((v * 1.04) / step) * step;
}

/* ───────────────────────── 3D cylinder ───────────────────────── */
class CylinderView {
  constructor(canvas, container) {
    this.canvas = canvas;
    this.container = container;
    this.renderer = new THREE.WebGLRenderer({ canvas, antialias: true, alpha: true });
    this.renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    this.scene = new THREE.Scene();
    this.camera = new THREE.PerspectiveCamera(34, 1, 0.1, 2000);
    this.camera.position.set(38, 26, 46);
    this.controls = new OrbitControls(this.camera, canvas);
    this.controls.enableDamping = true;
    this.controls.dampingFactor = 0.08;
    this.controls.enablePan = false;
    this.controls.autoRotate = true;
    this.controls.autoRotateSpeed = 0.6;
    canvas.addEventListener("pointerdown", () => { this.controls.autoRotate = false; }, { once: true });

    this.scene.add(new THREE.HemisphereLight(0xfff8ea, 0xcbbfa8, 1.4));
    const key = new THREE.DirectionalLight(0xffffff, 1.6);
    key.position.set(30, 50, 20);
    this.scene.add(key);
    const rim = new THREE.DirectionalLight(0xbcd0ff, 0.6);
    rim.position.set(-30, 10, -40);
    this.scene.add(rim);

    this.group = new THREE.Group();
    this.scene.add(this.group);
    this.waterColor = new THREE.Color("#1d5bd6");
    this.airColor = new THREE.Color("#e6d4b4");

    this.resize();
    new ResizeObserver(() => this.resize()).observe(container);
    this.animate = this.animate.bind(this);
    requestAnimationFrame(this.animate);
  }
  resize() {
    const w = this.container.clientWidth, h = this.container.clientHeight;
    if (!w || !h) return;
    this.renderer.setSize(w, h, false);
    this.camera.aspect = w / h;
    this.camera.updateProjectionMatrix();
  }
  animate() {
    this.controls.update();
    this.renderer.render(this.scene, this.camera);
    requestAnimationFrame(this.animate);
  }
  setModel(model) {
    const { R, p } = model;
    const L = p.L;
    this.group.clear();

    // sample pores to draw
    const rng = mulberry32(p.seed ^ 0x9e3779b9);
    const count = Math.min(MAX_DRAWN_PORES, model.N);
    this.poreR = new Float64Array(count);
    this.poreX = new Float64Array(count);
    this.poreZ = new Float64Array(count);
    this.poreV = new Float64Array(count);
    this.L = L;
    const geo = new THREE.CylinderGeometry(1, 1, 1, 10, 1, false);
    const mat = new THREE.MeshStandardMaterial({ roughness: 0.55, metalness: 0.05 });
    const mesh = new THREE.InstancedMesh(geo, mat, count);
    const dummy = new THREE.Object3D();
    const rRef = model.rMedian > 0 ? model.rMedian : p.rmean;
    for (let i = 0; i < count; i++) {
      const idx = Math.floor(rng() * model.N);
      const rr = model.r[idx];
      this.poreR[i] = rr;
      const rv = R * Math.min(0.08, Math.max(0.003, 0.011 * Math.sqrt(rr / rRef)));
      // uniform in disk, kept inside the shell
      const a = rng() * Math.PI * 2;
      const d = Math.sqrt(rng()) * (R - rv - 0.05 * R);
      this.poreX[i] = Math.cos(a) * d; this.poreZ[i] = Math.sin(a) * d; this.poreV[i] = rv;
      dummy.position.set(this.poreX[i], 0, this.poreZ[i]);
      dummy.scale.set(rv, L * 1.004, rv);
      dummy.updateMatrix();
      mesh.setMatrixAt(i, dummy.matrix);
      mesh.setColorAt(i, this.waterColor);
    }
    mesh.instanceMatrix.needsUpdate = true;
    this.pores = mesh;
    this.dummy = dummy;
    this.group.add(mesh);

    // translucent shell
    const shellGeo = new THREE.CylinderGeometry(R, R, L, 72, 1, false);
    const shellMat = new THREE.MeshPhysicalMaterial({
      color: 0xd8cbb0, transparent: true, opacity: 0.22, roughness: 0.35,
      transmission: 0, depthWrite: false, side: THREE.DoubleSide,
    });
    const shell = new THREE.Mesh(shellGeo, shellMat);
    shell.renderOrder = 2;
    this.group.add(shell);

    // ink outline on rims + a few generatrix lines
    const edges = new THREE.EdgesGeometry(new THREE.CylinderGeometry(R, R, L, 72, 1, true), 25);
    const outline = new THREE.LineSegments(edges, new THREE.LineBasicMaterial({ color: 0x221b12, transparent: true, opacity: 0.55 }));
    this.group.add(outline);
    const gen = [];
    for (let k = 0; k < 6; k++) {
      const a = (k / 6) * Math.PI * 2;
      gen.push(Math.cos(a) * R, -L / 2, Math.sin(a) * R, Math.cos(a) * R, L / 2, Math.sin(a) * R);
    }
    const genGeo = new THREE.BufferGeometry();
    genGeo.setAttribute("position", new THREE.Float32BufferAttribute(gen, 3));
    this.group.add(new THREE.LineSegments(genGeo, new THREE.LineBasicMaterial({ color: 0x221b12, transparent: true, opacity: 0.18 })));

    // frame the object
    const span = Math.max(L, 2 * R);
    const dist = span * 3.1;
    const dir = this.camera.position.clone().normalize();
    this.camera.position.copy(dir.multiplyScalar(dist));
    this.controls.minDistance = span * 1.2;
    this.controls.maxDistance = span * 6;
    this.controls.target.set(0, 0, 0);
    this.controls.update();
  }
  setHead(rstar) {
    if (!this.pores) return;
    const n = this.poreR.length;
    const d = this.dummy;
    for (let i = 0; i < n; i++) {
      const wet = this.poreR[i] <= rstar;
      this.pores.setColorAt(i, wet ? this.waterColor : this.airColor);
      // a drained pore collapses to a thin sand-coloured thread so the column visibly empties
      const rv = this.poreV[i] * (wet ? 1 : 0.35);
      d.position.set(this.poreX[i], 0, this.poreZ[i]);
      d.scale.set(rv, this.L * 1.004, rv);
      d.updateMatrix();
      this.pores.setMatrixAt(i, d.matrix);
    }
    this.pores.instanceColor.needsUpdate = true;
    this.pores.instanceMatrix.needsUpdate = true;
  }
}

/* ───────────────────────── UI wiring ───────────────────────── */
const $ = (sel) => document.querySelector(sel);
const form = $("#params");
const headSlider = $("#head");
const chart = new Chart($("#swrc-svg"), $("#tooltip"));
const view = new CylinderView($("#gl"), $("#viewport"));
let model = null;
let sweeping = null;

function readParams() {
  const f = form.elements;
  const num = (name, fallback) => {
    const v = parseFloat(f[name].value);
    return Number.isFinite(v) ? v : fallback;
  };
  const p = {
    dist: f.dist.value,
    rmean: Math.max(1e-5, num("rmean", DEFAULTS.rmean)),
    shape: num("shape", DEFAULTS.shape),
    sigln: num("sigln", DEFAULTS.sigln),
    cutoff: f.cutoff.checked,
    rmax: Math.max(1e-4, num("rmax", DEFAULTS.rmax)),
    n: Math.round(Math.min(2_000_000, Math.max(1000, num("n", DEFAULTS.n)))),
    seed: Math.round(num("seed", DEFAULTS.seed)),
    D: Math.max(0.5, num("D_num", DEFAULTS.D)),
    L: Math.max(0.5, num("L_num", DEFAULTS.L)),
    sigma: num("sigma", DEFAULTS.sigma),
    gamma: num("gamma", DEFAULTS.gamma),
    rho: num("rho", DEFAULTS.rho),
    g: num("g", DEFAULTS.g),
  };
  return p;
}

function writeParams(p) {
  const f = form.elements;
  f.dist.value = p.dist;
  f.rmean.value = p.rmean; f.rmean_log.value = Math.log10(p.rmean);
  f.shape.value = p.shape; f.shape_out.value = p.shape.toFixed(1);
  f.sigln.value = p.sigln; f.sigln_out.value = p.sigln.toFixed(2);
  f.cutoff.checked = p.cutoff;
  f.rmax.value = p.rmax; f.rmax_log.value = Math.log10(p.rmax);
  f.n.value = p.n; f.n_log.value = Math.log10(p.n);
  f.seed.value = p.seed;
  f.D.value = p.D; f.D_num.value = p.D;
  f.L.value = p.L; f.L_num.value = p.L;
  f.sigma.value = p.sigma; f.gamma.value = p.gamma; f.rho.value = p.rho; f.g.value = p.g;
  syncDistFields();
}

function syncDistFields() {
  const dist = form.elements.dist.value;
  form.querySelectorAll("[data-only]").forEach((n) => { n.hidden = n.dataset.only !== dist; });
  const cut = form.elements.cutoff.checked;
  form.elements.rmax.disabled = !cut; form.elements.rmax_log.disabled = !cut;
}

// log sliders ↔ number inputs
form.addEventListener("input", (e) => {
  const t = e.target;
  const f = form.elements;
  if (t.dataset.logFor) {
    const v = Math.pow(10, parseFloat(t.value));
    const target = f[t.dataset.logFor];
    target.value = t.dataset.logFor === "n" ? Math.round(v) : Number(v.toPrecision(3));
  } else if (t.name === "rmean" || t.name === "rmax" || t.name === "n") {
    const v = parseFloat(t.value);
    if (v > 0) f[t.name + "_log"].value = Math.log10(v);
  } else if (t.name === "D" || t.name === "L") {
    f[t.name + "_num"].value = t.value;
  } else if (t.dataset.mirror) {
    f[t.dataset.mirror].value = t.value;
  } else if (t.name === "shape") {
    f.shape_out.value = parseFloat(t.value).toFixed(1);
  } else if (t.name === "sigln") {
    f.sigln_out.value = parseFloat(t.value).toFixed(2);
  }
  if (t.name === "dist" || t.name === "cutoff") syncDistFields();
  // cylinder geometry & fluid constants are cheap: re-run live; the ensemble needs the button
  if (["D", "D_num", "L", "L_num", "sigma", "gamma", "rho", "g"].includes(t.name)) {
    scheduleLiveRun();
  }
  markDirty(true);
});
form.addEventListener("change", (e) => {
  // selects/checkboxes fire change only in some browsers; mirror the input handler's side effects
  if (e.target.name === "dist" || e.target.name === "cutoff") syncDistFields();
});
form.addEventListener("submit", (e) => { e.preventDefault(); run(); });

let liveTimer = null;
function scheduleLiveRun() {
  clearTimeout(liveTimer);
  liveTimer = setTimeout(() => run({ quiet: true }), 160);
}
function markDirty(d) {
  $("#run").classList.toggle("is-dirty", d);
  $("#run").textContent = d ? "Re-run model ↻" : "Re-run model";
}

// head slider ↔ |h| on a log scale over the model's range (0 at the saturated end)
function headFromSlider() {
  if (!model) return 0;
  const t = headSlider.valueAsNumber / 1000;
  if (t <= 0) return 0;
  const lo = Math.log10(model.hLo), hi = Math.log10(model.hHi);
  return Math.pow(10, lo + t * (hi - lo));
}
function sliderFromHead(h) {
  if (!model || h <= 0) return 0;
  const lo = Math.log10(model.hLo), hi = Math.log10(model.hHi);
  return Math.round(1000 * Math.min(1, Math.max(0, (Math.log10(h) - lo) / (hi - lo))));
}

function updateHead() {
  if (!model) return;
  const h = headFromSlider();
  const st = model.stateAt(h);
  $("#head-val").textContent = fmt.head(h);
  $("#theta").textContent = fmt.sig(st.theta, 4);
  $("#sat").textContent = fmt.pct(st.sat);
  $("#rstar").textContent = isFinite(st.rstar) ? fmt.sig(st.rstar) : "∞";
  $("#drained").textContent = fmt.int(st.drained);
  $("#drained-pct").textContent = `${fmt.pct(st.drained / model.N)} of N`;
  chart.setMarker(st);
  view.setHead(st.rstar);
}
for (const ev of ["input", "change"]) headSlider.addEventListener(ev, () => { stopSweep(); updateHead(); });

function run({ quiet = false } = {}) {
  stopSweep();
  const p = readParams();
  const runBtn = $("#run");
  runBtn.disabled = true;
  if (!quiet) $("#status").textContent = `sampling ${p.n.toLocaleString()} pores…`;
  // let the UI paint before the heavy loop
  setTimeout(() => {
    const t0 = performance.now();
    const prevH = model ? headFromSlider() : 0;
    model = runModel(p);
    const ms = performance.now() - t0;
    chart.setModel(model);
    view.setModel(model);
    $("#phi").textContent = fmt.sig(model.phi, 4);
    $("#hae").textContent = fmt.sig(model.hAirEntry);
    $("#legend-note").textContent =
      `${Math.min(MAX_DRAWN_PORES, model.N).toLocaleString()} of ${model.N.toLocaleString()} pores drawn · r̄ realized ${fmt.sig(model.realizedMean)} mm · r_max ${fmt.sig(model.rMax)} mm`;

    const warn = $("#warn");
    const msgs = [];
    if (model.phi > 1) {
      msgs.push(`Porosity ${fmt.sig(model.phi, 3)} exceeds 1: the summed pore area (${fmt.sig(Math.PI * model.phi * model.R * model.R, 3)} mm²) is larger than the cylinder cross-section (${fmt.sig(Math.PI * model.R * model.R, 3)} mm²), so the pores would overlap. Lower N or r̄ for a physical sample.`);
    }
    if (model.clamped > 0) {
      msgs.push(`${model.clamped.toLocaleString()} samples were clamped to r_max because the cutoff sits well inside the distribution; the realized mean is ${fmt.sig(model.realizedMean)} mm.`);
    }
    warn.hidden = msgs.length === 0;
    warn.textContent = msgs.join(" ");

    headSlider.value = sliderFromHead(prevH);
    updateHead();
    runBtn.disabled = false;
    markDirty(false);
    $("#status").textContent = `C = 2σcosγ/(ρg) = ${fmt.sig(model.C, 4)} mm² → h ≈ −${fmt.sig(model.C, 3)}/r · ${model.N.toLocaleString()} pores in ${ms.toFixed(0)} ms`;
  }, 30);
}

/* ───────────────────────── actions ───────────────────────── */
$("#run").addEventListener("click", () => run());
$("#reseed").addEventListener("click", () => {
  form.elements.seed.value = Math.floor(Math.random() * 1e6);
  run();
});
$("#reset").addEventListener("click", () => {
  writeParams({ ...DEFAULTS });
  run();
});
$("#csv").addEventListener("click", () => {
  if (!model) return;
  const rows = [["h_mm", "abs_h_mm", "theta", "saturation", "r_star_mm"]];
  for (const c of model.curve) {
    const st = model.stateAt(c.h);
    rows.push([-c.h, c.h, st.theta, st.sat, st.rstar].map((v) => Number(v.toPrecision(6))));
  }
  const meta = `# porosity=${model.phi}, N=${model.N}, dist=${model.p.dist}, rmean_mm=${model.p.rmean}, rmax_mm=${model.p.cutoff ? model.p.rmax : "none"}, D_mm=${model.p.D}, L_mm=${model.p.L}, sigma=${model.p.sigma}, gamma_deg=${model.p.gamma}, rho=${model.p.rho}, g=${model.p.g}, seed=${model.p.seed}\n`;
  const blob = new Blob([meta + rows.map((r) => r.join(",")).join("\n")], { type: "text/csv" });
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob);
  a.download = `swrc_${model.p.dist}_N${model.N}_seed${model.p.seed}.csv`;
  a.click();
  setTimeout(() => URL.revokeObjectURL(a.href), 1000);
});

function stopSweep() {
  if (sweeping) { cancelAnimationFrame(sweeping); sweeping = null; }
  $("#sweep").setAttribute("aria-pressed", "false");
  $("#sweep").textContent = "Drain sweep";
}
$("#sweep").addEventListener("click", () => {
  if (sweeping) { stopSweep(); return; }
  if (!model) return;
  $("#sweep").setAttribute("aria-pressed", "true");
  $("#sweep").textContent = "Stop sweep";
  const start = performance.now();
  const duration = 9000;
  const step = (now) => {
    const t = ((now - start) % duration) / duration;
    headSlider.value = Math.round(t * 1000);
    updateHead();
    sweeping = requestAnimationFrame(step);
  };
  sweeping = requestAnimationFrame(step);
});

/* ───────────────────────── boot ───────────────────────── */
writeParams({ ...DEFAULTS });
run();
