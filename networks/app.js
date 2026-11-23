import * as THREE from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";
import { DEFAULTS, MODELS, MODEL_LABEL, KIND, mualemFromRetention } from "./model.js";

const FILM_PRESET = { loamy: 0.15, silt: 0.33 };
const PRESETS = {
  // film water: median θ_nc fitted to the Hohenbrink et al. (2023) samples of each class
  coarse: { d50: 0.8, dist: "exponential", film: 0, phi: 0.35 },
  fine: { d50: 0.2, dist: "exponential", film: 0, phi: 0.35 },
  loamy: { d50: 0.1, dist: "lognormal", sigln: 0.9, film: FILM_PRESET.loamy, phi: 0.25 },
  silt: { d50: 0.03, dist: "lognormal", sigln: 1.2, film: FILM_PRESET.silt, phi: 0.12 },
};

const NOTES = {
  bundle:
    "<strong>Baseline (HW 2 and 3).</strong> Straight capillaries run the full length and never touch. Nothing moves sideways, so K<sub>h</sub> = 0. K<sub>sat</sub> matches ρg/8μ · Σr⁴/R² exactly.",
  plates:
    "<strong>Idea 1.</strong> Glass plates one grain thick, each drilled with random holes. Holes in touching plates connect only where their circles overlap, through a Sampson orifice (Q = ΔP r³/3μ). This builds Childs and Collis-George's (1950) cut-and-random-rejoin explicitly, the picture behind Mualem (1976).",
  gap:
    "<strong>Idea 2, Guido's version.</strong> The same plates, separated by a rough gap of aperture b. Every hole drains into the gap, which carries water sideways at b³/12μ per unit width and drains once suction passes 2σ/ρgb. b → 0 recovers idea 1. The defaults (r̄ = 0.10 d<sub>50</sub>, b = 0.26 d<sub>50</sub>) are calibrated so K<sub>h</sub> ≈ K<sub>v</sub> and K<sub>sat</sub> is within 1.5× of Kozeny–Carman.",
  tree:
    "<strong>Idea 3.</strong> Vertical trunks sprout side branches that keep bifurcating by Murray's law (r<sub>child</sub> = 2<sup>−1/3</sup> r). A pure tree has no loops, the Bethe-lattice limit, so a drained branch cuts off everything behind it. Tip fusion joins tips to nearby branches and closes loops.",
  lattice:
    "<strong>Reference.</strong> A classic pore-throat network (Fatt 1956) on a cubic lattice, thinned to coordination z. Guido noted that randomly rotated tubes converge to this. It is the model the new ones have to beat.",
};

const MODEL_COLOR = { bundle: "#6e7178", plates: "#c2550c", gap: "#00789a", tree: "#3d8b37", lattice: "#6d4fc0" };
const WATER = new THREE.Color("#00789a");
const SAND = new THREE.Color("#ddd0b6");
const STILL = new THREE.Color("#b4c3cb");
const RAMP = [new THREE.Color("#cfe6ee"), new THREE.Color("#2b9bbb"), new THREE.Color("#0a2f5c")];
const FLOW_DECADES = 4;

/* ───────────────────────── 3D view ───────────────────────── */
const Y_UP = new THREE.Vector3(0, 1, 0);
const Z_AX = new THREE.Vector3(0, 0, 1);

class NetworkView {
  constructor(canvas, container) {
    this.container = container;
    this.renderer = new THREE.WebGLRenderer({ canvas, antialias: true, alpha: true });
    this.renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    this.scene = new THREE.Scene();
    this.camera = new THREE.PerspectiveCamera(32, 1, 0.001, 1000);
    this.camera.position.set(2.2, 1.3, 2.6);
    this.controls = new OrbitControls(this.camera, canvas);
    this.controls.enableDamping = true;
    this.controls.dampingFactor = 0.08;
    this.controls.enablePan = false;
    this.controls.autoRotate = true;
    this.controls.autoRotateSpeed = 0.5;
    canvas.addEventListener("pointerdown", () => { this.controls.autoRotate = false; }, { once: true });

    this.scene.add(new THREE.HemisphereLight(0xffffff, 0xc9ccd0, 1.5));
    const key = new THREE.DirectionalLight(0xffffff, 1.5);
    key.position.set(3, 5, 2);
    this.scene.add(key);
    const rim = new THREE.DirectionalLight(0xcfe3ff, 0.6);
    rim.position.set(-3, 1, -4);
    this.scene.add(rim);

    this.group = new THREE.Group();
    this.scene.add(this.group);
    this.dummy = new THREE.Object3D();
    this.color = new THREE.Color();
    this.resize();
    new ResizeObserver(() => this.resize()).observe(container);
    const tick = () => { this.controls.update(); this.renderer.render(this.scene, this.camera); requestAnimationFrame(tick); };
    requestAnimationFrame(tick);
  }
  resize() {
    const w = this.container.clientWidth, h = this.container.clientHeight;
    if (!w || !h) return;
    this.renderer.setSize(w, h, false);
    this.camera.aspect = w / h;
    this.camera.updateProjectionMatrix();
  }
  clear() {
    this.group.traverse((o) => { o.geometry?.dispose(); o.material?.dispose?.(); });
    this.group.clear();
  }
  setNetwork(net, sameFrame, cut = this.cut) {
    this.clear();
    this.cut = cut;
    const { meta } = net;
    const R = meta.model === "lattice" ? meta.Rh : meta.R;
    const L = meta.L;
    const minR = R * 0.0035;
    // model (x, y, z) with z up → three (x, z, y)
    const P = (n, v) => v.set(net.x[n], net.z[n], net.y[n]);

    const tubes = [], gaps = [];
    for (let i = 0; i < net.nBonds; i++) {
      const k = net.kind[i];
      if (k === KIND.ORIFICE || k === KIND.JUNCTION) continue;
      (k === KIND.GAP ? gaps : tubes).push(i);
    }
    this.tubes = this.makeInstances(net, tubes, new THREE.CylinderGeometry(1, 1, 1, 8, 1, true), (i, pa, pb, d) => ({
      quat: new THREE.Quaternion().setFromUnitVectors(Y_UP, d), sx: Math.max(net.r[i], minR), sy: pa.distanceTo(pb),
    }), { roughness: 0.5, metalness: 0.05 });
    this.gaps = this.makeInstances(net, gaps, new THREE.BoxGeometry(1, 1, 1), (i, pa, pb, d) => ({
      quat: new THREE.Quaternion().setFromUnitVectors(Z_AX, d), sx: net.w[i], sy: Math.max(net.r[i], minR * 0.8), sz: pa.distanceTo(pb),
    }), { roughness: 0.3, transparent: true, opacity: 0.55, depthWrite: false });

    // glass plates
    if (meta.model === "plates" || meta.model === "gap") {
      const t = meta.t, b0 = meta.b0 || 0;
      const n = Math.round((L + b0) / (t + b0));
      const plateMat = new THREE.MeshPhysicalMaterial({
        color: 0xcdd6db, transparent: true, opacity: 0.07, roughness: 0.2, depthWrite: false, side: THREE.DoubleSide,
      });
      const ringMat = new THREE.LineBasicMaterial({ color: 0x1b1d21, transparent: true, opacity: 0.22 });
      // with the cutaway on, the +x/+z quadrant of every plate is removed
      const a0 = cut ? Math.PI / 2 : 0, span = cut ? 1.5 * Math.PI : 2 * Math.PI;
      const plateGeo = new THREE.CylinderGeometry(R, R, t, 64, 1, false, a0, span);
      const ringGeo = new THREE.BufferGeometry().setFromPoints(
        Array.from({ length: 97 }, (_, k) => {
          const th = a0 + (k / 96) * span;
          return new THREE.Vector3(R * Math.sin(th), 0, R * Math.cos(th));
        }),
      );
      for (let k = 0; k < n; k++) {
        const yc = L / 2 - k * (t + b0) - t / 2;
        const m = new THREE.Mesh(plateGeo, plateMat);
        m.position.y = yc;
        m.renderOrder = 3;
        this.group.add(m);
        for (const dy of [t / 2, -t / 2]) {
          const ring = new THREE.Line(ringGeo, ringMat);
          ring.position.y = yc + dy;
          this.group.add(ring);
        }
      }
    }

    // glass cylinder, as in HW 2
    const pad = R * 1.04;
    const shell = new THREE.Mesh(
      new THREE.CylinderGeometry(pad, pad, L * 1.02, 72, 1, true),
      new THREE.MeshPhysicalMaterial({ color: 0xd6dde1, transparent: true, opacity: 0.16, roughness: 0.3, depthWrite: false, side: THREE.DoubleSide }),
    );
    shell.renderOrder = 4;
    this.group.add(shell);
    const edges = new THREE.EdgesGeometry(new THREE.CylinderGeometry(pad, pad, L * 1.02, 72, 1, true), 25);
    this.group.add(new THREE.LineSegments(edges, new THREE.LineBasicMaterial({ color: 0x1b1d21, transparent: true, opacity: 0.5 })));
    const gen = [];
    for (let k = 0; k < 6; k++) {
      const a = (k / 6) * Math.PI * 2;
      gen.push(Math.cos(a) * pad, -L * 0.51, Math.sin(a) * pad, Math.cos(a) * pad, L * 0.51, Math.sin(a) * pad);
    }
    const genGeo = new THREE.BufferGeometry();
    genGeo.setAttribute("position", new THREE.Float32BufferAttribute(gen, 3));
    this.group.add(new THREE.LineSegments(genGeo, new THREE.LineBasicMaterial({ color: 0x1b1d21, transparent: true, opacity: 0.14 })));

    // electrodes for the horizontal solve: −x rim at high head, +x rim at low head
    this.electrodes = new THREE.Group();
    const th = (40 * Math.PI) / 180;
    for (const [c0, col] of [[Math.PI, 0xc2550c], [0, 0x1b1d21]]) {
      const arc = new THREE.Mesh(
        new THREE.CylinderGeometry(pad * 1.01, pad * 1.01, L * 1.02, 24, 1, true, Math.PI / 2 - c0 - th, 2 * th),
        new THREE.MeshBasicMaterial({ color: col, transparent: true, opacity: 0.12, side: THREE.DoubleSide, depthWrite: false }),
      );
      this.electrodes.add(arc);
    }
    this.electrodes.visible = false;
    this.group.add(this.electrodes);

    if (!sameFrame) {
      const span = Math.max(L, 2 * R);
      const dir = this.camera.position.clone().normalize();
      this.camera.position.copy(dir.multiplyScalar(span * 2.5));
      this.camera.near = span / 100;
      this.camera.far = span * 50;
      this.camera.updateProjectionMatrix();
      this.controls.minDistance = span * 0.8;
      this.controls.maxDistance = span * 8;
      this.controls.target.set(0, 0, 0);
    }
    this.controls.update();
  }
  makeInstances(net, ids, geo, place, matOpts) {
    const mesh = new THREE.InstancedMesh(geo, new THREE.MeshStandardMaterial(matOpts), Math.max(ids.length, 1));
    mesh.count = ids.length;
    const pos = new Float32Array(ids.length * 3), quat = new Float32Array(ids.length * 4), scl = new Float32Array(ids.length * 3);
    const pa = new THREE.Vector3(), pb = new THREE.Vector3(), d = new THREE.Vector3();
    ids.forEach((i, n) => {
      pa.set(net.x[net.a[i]], net.z[net.a[i]], net.y[net.a[i]]);
      pb.set(net.x[net.b[i]], net.z[net.b[i]], net.y[net.b[i]]);
      d.subVectors(pb, pa).normalize();
      const s = place(i, pa, pb, d);
      pos[3 * n] = (pa.x + pb.x) / 2; pos[3 * n + 1] = (pa.y + pb.y) / 2; pos[3 * n + 2] = (pa.z + pb.z) / 2;
      quat[4 * n] = s.quat.x; quat[4 * n + 1] = s.quat.y; quat[4 * n + 2] = s.quat.z; quat[4 * n + 3] = s.quat.w;
      scl[3 * n] = s.sx; scl[3 * n + 1] = s.sy; scl[3 * n + 2] = s.sz ?? s.sx;
    });
    this.group.add(mesh);
    return { mesh, ids, pos, quat, scl, isGap: geo.type === "BoxGeometry" };
  }
  update(wet, flux, mode) {
    this.electrodes.visible = mode === "h";
    let qmax = 0;
    if (flux) for (let i = 0; i < flux.length; i++) { const q = Math.abs(flux[i]); if (q > qmax) qmax = q; }
    for (const set of [this.tubes, this.gaps]) {
      const { mesh, ids, pos, quat, scl, isGap } = set;
      const d = this.dummy;
      ids.forEach((i, n) => {
        let thin = 1;
        if (this.cut && pos[3 * n] > 0 && pos[3 * n + 2] > 0) thin = 0;
        else if (!wet[i]) { this.color.copy(SAND); thin = isGap ? 0 : 0.35; }
        else if (mode === "water") this.color.copy(WATER);
        else {
          const q = Math.abs(flux[i]);
          if (!(qmax > 0) || q < qmax * 10 ** -FLOW_DECADES) { this.color.copy(STILL); thin = isGap ? 0.25 : 0.6; }
          else {
            const u = 1 + Math.log10(q / qmax) / FLOW_DECADES; // 0 … 1
            if (u < 0.5) this.color.lerpColors(RAMP[0], RAMP[1], u * 2);
            else this.color.lerpColors(RAMP[1], RAMP[2], (u - 0.5) * 2);
          }
        }
        d.position.set(pos[3 * n], pos[3 * n + 1], pos[3 * n + 2]);
        d.quaternion.set(quat[4 * n], quat[4 * n + 1], quat[4 * n + 2], quat[4 * n + 3]);
        if (isGap) d.scale.set(scl[3 * n] * (thin ? 1 : 0), scl[3 * n + 1] * thin, scl[3 * n + 2] * (thin ? 1 : 0));
        else d.scale.set(scl[3 * n] * thin, scl[3 * n + 1] * (thin ? 1 : 0), scl[3 * n + 2] * thin);
        d.updateMatrix();
        mesh.setMatrixAt(n, d.matrix);
        mesh.setColorAt(n, this.color);
      });
      mesh.instanceMatrix.needsUpdate = true;
      if (mesh.instanceColor) mesh.instanceColor.needsUpdate = true;
    }
  }
}

/* ───────────────────────── charts ───────────────────────── */
const SVGNS = "http://www.w3.org/2000/svg";
function el(name, attrs = {}, parent) {
  const e = document.createElementNS(SVGNS, name);
  for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v);
  if (parent) parent.appendChild(e);
  return e;
}
function axisScale(ax, a, b) {
  if (ax.log) {
    const l0 = Math.log10(ax.min), l1 = Math.log10(ax.max);
    return (v) => a + ((Math.log10(Math.max(v, ax.min)) - l0) / (l1 - l0)) * (b - a);
  }
  return (v) => a + ((v - ax.min) / (ax.max - ax.min)) * (b - a);
}
function ticks(ax) {
  if (ax.log) {
    const out = [];
    for (let e = Math.ceil(Math.log10(ax.min) - 1e-9); e <= Math.floor(Math.log10(ax.max) + 1e-9); e++) out.push(10 ** e);
    return out;
  }
  const span = ax.max - ax.min, raw = span / 5, mag = 10 ** Math.floor(Math.log10(raw));
  const step = [1, 2, 2.5, 5, 10].map((m) => m * mag).find((s) => s >= raw);
  const out = [];
  for (let v = Math.ceil(ax.min / step) * step; v <= ax.max + 1e-12; v += step) out.push(+v.toPrecision(6));
  return out;
}
function fmtTick(v, ax) {
  if (ax.log) {
    const e = Math.round(Math.log10(v));
    if (e >= 0 && e <= 4) return String(10 ** e);
    return `10${String(e).replace("-", "⁻").replace(/\d/g, (d) => "⁰¹²³⁴⁵⁶⁷⁸⁹"[d])}`;
  }
  return ax.fmt ? ax.fmt(v) : String(v);
}
function plot(svg, cfg) {
  const W = svg.clientWidth || 400, H = svg.clientHeight || 200;
  svg.setAttribute("viewBox", `0 0 ${W} ${H}`);
  svg.replaceChildren();
  const m = { l: 50, r: 10, t: 8, b: 30 };
  const sx = axisScale(cfg.x, m.l, W - m.r), sy = axisScale(cfg.y, H - m.b, m.t);
  const grid = el("g", { class: "grid" }, svg), axis = el("g", { class: "axis" }, svg);
  for (const v of ticks(cfg.x)) {
    el("line", { x1: sx(v), x2: sx(v), y1: m.t, y2: H - m.b }, grid);
    el("text", { x: sx(v), y: H - m.b + 14, "text-anchor": "middle", class: "tick-label" }, axis).textContent = fmtTick(v, cfg.x);
  }
  for (const v of ticks(cfg.y)) {
    el("line", { x1: m.l, x2: W - m.r, y1: sy(v), y2: sy(v) }, grid);
    el("text", { x: m.l - 5, y: sy(v) + 3.5, "text-anchor": "end", class: "tick-label" }, axis).textContent = fmtTick(v, cfg.y);
  }
  el("path", { d: `M${m.l},${m.t} V${H - m.b} H${W - m.r}` }, axis);
  el("text", { x: (m.l + W - m.r) / 2, y: H - 3, "text-anchor": "middle", class: "axis-label" }, svg).textContent = cfg.x.label;
  el("text", { x: 11, y: (m.t + H - m.b) / 2, "text-anchor": "middle", class: "axis-label", transform: `rotate(-90 11 ${(m.t + H - m.b) / 2})` }, svg).textContent = cfg.y.label;
  const clip = `c${Math.random().toString(36).slice(2)}`;
  el("rect", { x: m.l, y: m.t, width: W - m.l - m.r, height: H - m.t - m.b }, el("clipPath", { id: clip }, el("defs", {}, svg)));
  const body = el("g", { "clip-path": `url(#${clip})` }, svg);
  for (const s of cfg.series) {
    let d = "", pen = false;
    for (const [x, y] of s.pts) {
      const ok = Number.isFinite(x) && Number.isFinite(y) && (!cfg.y.log || y > 0) && (!cfg.x.log || x > 0);
      if (!ok) {
        // a curve that reaches zero on a log axis falls off the bottom
        if (pen && cfg.y.log && y <= 0) d += `V${H - m.b + 4}`;
        pen = false;
        continue;
      }
      d += `${pen ? "L" : "M"}${sx(x).toFixed(1)},${sy(Math.max(y, cfg.y.min * 0.5)).toFixed(1)}`;
      pen = true;
    }
    const path = el("path", { d, class: `series ${s.cls || ""}` }, body);
    if (s.color) path.style.stroke = s.color;
  }
  if (cfg.marker) {
    const [x, y] = cfg.marker;
    if (Number.isFinite(x) && Number.isFinite(y) && (!cfg.y.log || y > 0)) {
      el("line", { x1: sx(x), x2: sx(x), y1: m.t, y2: H - m.b, class: "marker-line" }, body);
      el("circle", { cx: sx(x), cy: sy(y), r: 4.5, class: "marker" }, body);
    }
  }
}

/* ───────────────────────── formatting ───────────────────────── */
const sci = (v, d = 2) => {
  if (!Number.isFinite(v)) return "–";
  if (v === 0) return "0";
  const e = Math.floor(Math.log10(Math.abs(v)));
  const m = v / 10 ** e;
  return `${m.toFixed(d - 1)}×10${String(e).replace("-", "⁻").replace(/\d/g, (c) => "⁰¹²³⁴⁵⁶⁷⁸⁹"[c])}`;
};
const fix = (v, d) => (Number.isFinite(v) ? v.toFixed(d) : "–");

/* ───────────────────────── state & wiring ───────────────────────── */
const $ = (s) => document.querySelector(s);
const form = $("#params");
const view = new NetworkView($("#gl"), $("#viewport"));
const worker = new Worker(new URL("./worker.js", import.meta.url), { type: "module" });
let params = { ...DEFAULTS, preset: "fine" };
let result = null;
let compare = null;
let mode = "water";
let cutaway = true;
view.cut = cutaway;
let reqId = 0, cmpId = 0;
let sweep = null;
let lastModel = null;

const modelBox = $("#models");
for (const m of MODELS) {
  const b = document.createElement("button");
  b.type = "button";
  b.className = "model-btn";
  b.dataset.model = m;
  b.innerHTML = `<i style="background:${MODEL_COLOR[m]}"></i>${MODEL_LABEL[m]}`;
  b.addEventListener("click", () => { params.model = m; syncForm(); run(); });
  modelBox.appendChild(b);
}

function syncForm() {
  for (const elx of form.elements) {
    if (!elx.name || !(elx.name in params)) continue;
    if (elx.type === "checkbox") elx.checked = !!params[elx.name];
    else elx.value = params[elx.name];
    const out = elx.parentElement.querySelector("output");
    if (out) out.textContent = formatParam(elx.name, params[elx.name]);
  }
  for (const b of modelBox.children) b.setAttribute("aria-pressed", String(b.dataset.model === params.model));
  for (const n of form.querySelectorAll("[data-models]")) n.hidden = !n.dataset.models.split(" ").includes(params.model);
  for (const n of form.querySelectorAll("[data-only]")) n.hidden = params.dist !== n.dataset.only;
  $("#note").innerHTML = NOTES[params.model];
}
function formatParam(name, v) {
  if (["N", "layers", "generations", "segsPerGen", "branchAngle"].includes(name)) return String(v);
  if (name === "bRatio" || name === "film") return Number(v).toFixed(3);
  return Number(v).toFixed(2);
}
function readField(elx) {
  if (elx.type === "checkbox") return elx.checked;
  if (elx.tagName === "SELECT") return elx.value;
  return Number(elx.value);
}

let timer = null;
form.addEventListener("input", (ev) => {
  const t = ev.target;
  if (!t.name) return;
  if (t.name === "preset") {
    params.preset = t.value;
    Object.assign(params, PRESETS[t.value] || {});
  } else {
    params[t.name] = readField(t);
    if (["d50", "dist", "sigln"].includes(t.name)) params.preset = "custom";
  }
  syncForm();
  clearTimeout(timer);
  timer = setTimeout(run, t.type === "range" ? 280 : 120);
});
$("#reset").addEventListener("click", () => { params = { ...DEFAULTS, preset: "fine", model: params.model }; syncForm(); run(); });
$("#reseed").addEventListener("click", () => { params.seed = (params.seed * 7919 + 13) % 99991; syncForm(); run(); });

function run() {
  stopSweep();
  const id = ++reqId;
  $("#busy").hidden = false;
  $("#status").textContent = `Building ${MODEL_LABEL[params.model]} and solving ${params.steps + 1} drainage steps in two directions…`;
  worker.postMessage({ id, type: "run", params: { ...params } });
  if (compare) $("#cmp-status").textContent = "Parameters changed since the last comparison. Run again to update.";
}

worker.onmessage = (ev) => {
  const msg = ev.data;
  if (msg.type === "error") { $("#busy").hidden = true; $("#status").textContent = `Error: ${msg.message}`; return; }
  if (msg.type === "progress") { $("#cmp-status").textContent = `Solving ${MODEL_LABEL[msg.model]}…`; return; }
  if (msg.type === "compare") { if (msg.id === cmpId) { compare = msg.results; drawCompare(); } return; }
  if (msg.type !== "run" || msg.id !== reqId) return;
  $("#busy").hidden = true;
  result = msg;
  const sameModel = lastModel === msg.net.meta.model;
  lastModel = msg.net.meta.model;
  view.setNetwork(msg.net, sameModel);
  const slider = $("#step");
  slider.max = String(msg.steps.length - 1);
  if (Number(slider.value) > msg.steps.length - 1) slider.value = slider.max;
  const s = msg.summary;
  $("#net-title").textContent = MODEL_LABEL[s.model];
  $("#phi").textContent = fix(s.phi, 3);
  $("#size").textContent = `${s.nodes.toLocaleString()} nodes · ${s.bonds.toLocaleString()} bonds`;
  $("#ksat").textContent = sci(s.Ksat);
  $("#kc").textContent = `m/s · ${fix(s.Ksat / s.KC, 2)} × Kozeny–Carman`;
  $("#muerr").textContent = `Mualem assumes 0.5 · off by ${fix(s.mualemErr, 2)} decades`;
  $("#sc-note").textContent = params.corner > 0 ? "corners keep every pore connected" : "saturation where K_v → 0";
  $("#aniso").textContent = s.KsatH > 0 ? fix(s.aniso, s.aniso < 0.01 ? 4 : 3) : "0";
  $("#ksath").textContent = s.KsatH > 0 ? `K_h = ${sci(s.KsatH)} m/s` : "no sideways path";
  $("#lfit").textContent = fix(s.lFit, 2);
  $("#hae").textContent = fix(s.hAE, 0);
  $("#sc").textContent = s.Sc > 0 ? `S = ${fix(s.Sc, 2)}` : "never";
  $("#status").textContent = `${s.nodes.toLocaleString()} nodes, ${s.bonds.toLocaleString()} bonds, ${msg.steps.length} drainage steps × 2 directions in ${(msg.ms / 1000).toFixed(1)} s.`;
  updateStep();
};

function updateStep() {
  if (!result) return;
  const s = Number($("#step").value);
  const B = result.net.nBonds;
  const st = result.steps[s];
  const wet = result.masks.subarray(s * B, (s + 1) * B);
  const flux = mode === "v" ? result.fluxV.subarray(s * B, (s + 1) * B) : mode === "h" ? result.fluxH.subarray(s * B, (s + 1) * B) : null;
  view.update(wet, flux, mode);
  $("#h-val").textContent = fix(st.h, 0);
  $("#s-val").textContent = fix(st.S, 3);
  $("#kr-val").textContent = sci(st.Kv / result.steps[0].Kv);
  drawLegend();
  drawModelCharts(s);
}

function drawLegend() {
  const lg = $("#legend");
  if (mode === "water") {
    lg.innerHTML = `<span><i class="sw sw--water"></i>water-filled</span><span><i class="sw sw--air"></i>drained (air)</span>`;
  } else {
    lg.innerHTML = `<span><i class="sw sw--ramp"></i>flux, 10⁻⁴ → 1 of max</span><span><i class="sw sw--still"></i>water, no flow</span><span><i class="sw sw--air"></i>drained</span>` +
      (mode === "h" ? `<span class="legend__note">Water enters at the orange rim and leaves at the dark rim.</span>` : `<span class="legend__note">Water enters at the top face and leaves at the bottom.</span>`);
  }
}

function drawModelCharts(s) {
  const st = result.steps;
  const hs = st.filter((d) => d.h > 0).map((d) => d.h);
  const hMin = Math.min(...hs) * 0.7, hMax = Math.max(...hs) * 1.3;
  const phi = result.summary.phi;
  plot($("#ret-svg"), {
    x: { log: true, min: 10 ** Math.floor(Math.log10(hMin)), max: 10 ** Math.ceil(Math.log10(hMax)), label: "suction |h| (mm)" },
    y: { min: 0, max: Math.ceil(phi * 20) / 20, label: "θ" },
    series: [{ pts: st.map((d) => [Math.max(d.h, hMin), d.theta]) }],
    marker: [Math.max(st[s].h, hMin), st[s].theta],
  });
  const K0 = st[0].Kv, H0 = st[0].Kh;
  const l = result.summary.lFit;
  const fit = Number.isFinite(l) ? mualemFromRetention(st.map((d) => d.h), st.map((d) => d.S), l).Kr : [];
  plot($("#k-svg"), {
    x: { min: 0, max: 1, label: "saturation S" },
    y: { log: true, min: 1e-6, max: 1, label: "K / K_sat" },
    series: [
      { pts: st.map((d, i) => [d.S, result.mualem[i]]), cls: "series--pred" },
      { pts: st.map((d, i) => [d.S, fit[i]]), cls: "series--fit" },
      ...(H0 > 0 ? [{ pts: st.map((d) => [d.S, d.Kh / H0]), cls: "series--h" }] : []),
      { pts: st.map((d) => [d.S, d.Kv / K0]) },
    ],
    marker: [st[s].S, st[s].Kv / K0],
  });
}

function drawCompare() {
  if (!compare) return;
  $("#cmp-status").textContent = "Same soil, same cylinder, every model. Vertical K_r against saturation.";
  $("#cmp-key").innerHTML = MODELS.map((m) => `<li><i class="ln" style="border-color:${MODEL_COLOR[m]}"></i>${MODEL_LABEL[m]}</li>`).join("");
  plot($("#cmp-svg"), {
    x: { min: 0, max: 1, label: "saturation S" },
    y: { log: true, min: 1e-6, max: 1, label: "K / K_sat" },
    series: MODELS.map((m) => {
      const st = compare[m].steps;
      return { pts: st.map((d) => [d.S, d.Kv / st[0].Kv]), color: MODEL_COLOR[m] };
    }),
  });
  $("#cmp-table tbody").innerHTML = MODELS.map((m) => {
    const s = compare[m].summary;
    return `<tr class="${m === params.model ? "is-active" : ""}"><td>${MODEL_LABEL[m]}</td><td>${fix(s.phi, 2)}</td><td>${sci(s.Ksat)}</td><td>${s.KsatH > 0 ? fix(s.aniso, 3) : "0"}</td><td>${fix(s.lFit, 2)}</td><td>${s.Sc > 0 ? fix(s.Sc, 2) : "–"}</td></tr>`;
  }).join("");
}

$("#compare").addEventListener("click", () => {
  cmpId++;
  $("#cmp-status").textContent = "Starting…";
  worker.postMessage({ id: cmpId, type: "compare", params: { ...params } });
});

$("#step").addEventListener("input", () => { stopSweep(); updateStep(); });
for (const b of $("#colour").children)
  b.addEventListener("click", () => {
    mode = b.dataset.mode;
    for (const o of $("#colour").children) o.setAttribute("aria-pressed", String(o === b));
    updateStep();
  });

function stopSweep() {
  if (!sweep) return;
  clearInterval(sweep);
  sweep = null;
  $("#sweep").setAttribute("aria-pressed", "false");
}
$("#sweep").addEventListener("click", () => {
  if (sweep) return stopSweep();
  const slider = $("#step");
  if (Number(slider.value) >= Number(slider.max)) slider.value = "0";
  $("#sweep").setAttribute("aria-pressed", "true");
  sweep = setInterval(() => {
    const v = Number(slider.value) + 1;
    if (v > Number(slider.max)) return stopSweep();
    slider.value = String(v);
    updateStep();
  }, 160);
});

for (const t of document.querySelectorAll(".tab"))
  t.addEventListener("click", () => {
    for (const o of document.querySelectorAll(".tab")) o.setAttribute("aria-pressed", String(o === t));
    for (const body of document.querySelectorAll(".tab-body")) body.hidden = body.dataset.body !== t.dataset.tab;
    if (t.dataset.tab === "compare") drawCompare();
    else if (result) drawModelCharts(Number($("#step").value));
  });

$("#cut").addEventListener("click", () => {
  cutaway = !cutaway;
  $("#cut").setAttribute("aria-pressed", String(cutaway));
  if (result) { view.setNetwork(result.net, true, cutaway); updateStep(); }
});

$("#csv").addEventListener("click", () => {
  if (!result) return;
  const st = result.steps;
  const rows = ["h_mm,theta,S,Kv_m_s,Kh_m_s,Kr_v,Kr_h,Kr_mualem_l0.5"];
  st.forEach((d, i) => rows.push([d.h, d.theta, d.S, d.Kv, d.Kh, d.Kv / st[0].Kv, st[0].Kh > 0 ? d.Kh / st[0].Kh : 0, result.mualem[i]].map((v) => +v.toPrecision(6)).join(",")));
  const a = document.createElement("a");
  a.href = URL.createObjectURL(new Blob([rows.join("\n")], { type: "text/csv" }));
  a.download = `${params.model}-network-curves.csv`;
  a.click();
  URL.revokeObjectURL(a.href);
});

new ResizeObserver(() => {
  if (!result) return;
  if (!$("[data-body=model]").hidden) drawModelCharts(Number($("#step").value));
  else drawCompare();
}).observe($(".panel--curves"));

syncForm();
run();
