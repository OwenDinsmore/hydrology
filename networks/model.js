// Pore-network models for the three layered/branching ideas, plus two references.
//
// Every model is a graph of nodes joined by bonds. Bonds hold all the water and all
// the capillary barriers; nodes are zero-volume junctions that never block flow.
//
// Units: lengths in mm, heads in mm, conductances in mm³ with μ factored out
// (a tube is π r⁴ / 8ℓ). K_UNIT turns G·(length/area) into m/s.

export const FLUID = { sigma: 0.072, rho: 1000, g: 9.8, mu: 1.0e-3 };
export const K_UNIT = (FLUID.rho * FLUID.g / FLUID.mu) * 1e-6;
export const HEAD_CONST = (2 * FLUID.sigma * 1e6) / (FLUID.rho * FLUID.g); // h[mm] = HEAD_CONST / r[mm]

export const KIND = { TUBE: 0, ORIFICE: 1, GAP: 2, LATERAL: 3, TRUNK: 4, BRANCH: 5, THROAT: 6, JUNCTION: 7 };
export const MODELS = ["bundle", "plates", "gap", "tree", "lattice"];
export const MODEL_LABEL = {
  bundle: "Straight bundle",
  plates: "1 · Stacked plates",
  gap: "2 · Gapped plates",
  tree: "3 · Branching tree",
  lattice: "Pore-throat lattice",
};

export const CORNER_AREA = 1 - Math.PI / (3 * Math.sqrt(3)); // corner water of a triangle at r* = r
export const DEFAULTS = {
  model: "gap",
  d50: 0.2,        // median grain diameter, mm
  rRatio: 0.1,     // mean pore radius / d50 at φ = 0.35
  rmaxRatio: 3,    // largest pore / mean pore
  dist: "exponential",
  sigln: 0.6,
  phi: 0.35,       // porosity of the drilled holes (plates) or of the whole sample
  N: 220,          // pores per layer; sets the sample radius
  layers: 14,
  tRatio: 1,       // plate thickness / d50
  seed: 7,
  access: true,    // air must reach a pore through drained pores before it drains
  lateral: 0,      // drilled horizontal holes per vertical hole (idea 2a)
  bRatio: 0.26,    // gap aperture / d50 at φ = 0.35 (Guido's gap)
  voidScale: true, // pore radius and gap aperture grow with the void ratio φ/(1−φ)
  bSigma: 0.4,     // log-normal roughness of the gap aperture
  trunkFactor: 1.6,
  pBranch: 0.45,
  generations: 3,
  segsPerGen: 2,
  branchAngle: 35,
  sideRatio: 0.6,
  anastomosis: 0,
  coord: 4,
  film: 0,         // non-capillary water on grain surfaces at air entry, θ_nc (volumetric)
  corner: 0.00015, // conductance kept by a drained pore, × (r*/r)⁴; fitted to measured sands; 0 = round tubes
  steps: 26,
};

/* ───────────────────────── random numbers ───────────────────────── */
export function mulberry32(seed) {
  let a = seed >>> 0;
  return function () {
    a = (a + 0x6d2b79f5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}
function normal(rng) {
  let u = 0;
  while (u === 0) u = rng();
  return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * rng());
}
function poisson(rng, lam) {
  if (lam <= 0) return 0;
  const L = Math.exp(-lam);
  let k = 0, p = 1;
  do { k++; p *= rng(); } while (p > L);
  return k - 1;
}

/** Hydraulic radius of a grain pack scales with the void ratio; 1 at φ = 0.35. */
export function voidScale(p) {
  return p.voidScale ? (p.phi / (1 - p.phi)) / (0.35 / 0.65) : 1;
}

export function makeRadius(p, rng) {
  const rbar = p.rRatio * p.d50 * voidScale(p);
  const rmax = p.rmaxRatio * rbar;
  const rmin = 0.02 * rbar;
  if (p.dist === "lognormal") {
    const s = p.sigln, mu = Math.log(rbar) - 0.5 * s * s;
    return () => {
      for (;;) {
        const r = Math.exp(mu + s * normal(rng));
        if (r <= rmax && r >= rmin) return r;
      }
    };
  }
  return () => {
    for (;;) {
      const r = -rbar * Math.log(1 - rng());
      if (r <= rmax && r >= rmin) return r;
    }
  };
}

/** Shared sample size so every model sees the same cylinder. */
export function geometry(p, scale = 1) {
  const draw = makeRadius(p, mulberry32(0x51ed ^ p.seed));
  let s2 = 0;
  const n = 20000;
  for (let i = 0; i < n; i++) { const r = draw(); s2 += r * r; }
  const Er2 = s2 / n;
  const R = scale * Math.sqrt((p.N * Er2) / p.phi);
  const t = p.tRatio * p.d50;
  return { R, t, L: p.layers * t, Er2, scale };
}

/* ───────────────────────── spatial hashes ───────────────────────── */
class Hash2 {
  constructor(cell) { this.c = cell; this.m = new Map(); }
  key(i, j) { return (i + 50000) * 100003 + (j + 50000); }
  add(x, y, id) {
    const k = this.key(Math.floor(x / this.c), Math.floor(y / this.c));
    let a = this.m.get(k);
    if (!a) this.m.set(k, (a = []));
    a.push(id);
  }
  each(x, y, rad, fn) {
    const c = this.c;
    const i0 = Math.floor((x - rad) / c), i1 = Math.floor((x + rad) / c);
    const j0 = Math.floor((y - rad) / c), j1 = Math.floor((y + rad) / c);
    for (let i = i0; i <= i1; i++)
      for (let j = j0; j <= j1; j++) {
        const a = this.m.get(this.key(i, j));
        if (a) for (const id of a) fn(id);
      }
  }
}
class Hash3 {
  constructor(cell) { this.c = cell; this.m = new Map(); }
  key(i, j, k) { return ((i + 1000) * 2003 + (j + 1000)) * 2003 + (k + 1000); }
  add(x, y, z, id) {
    const c = this.c, k = this.key(Math.floor(x / c), Math.floor(y / c), Math.floor(z / c));
    let a = this.m.get(k);
    if (!a) this.m.set(k, (a = []));
    a.push(id);
  }
  each(x, y, z, rad, fn) {
    const c = this.c;
    for (let i = Math.floor((x - rad) / c); i <= Math.floor((x + rad) / c); i++)
      for (let j = Math.floor((y - rad) / c); j <= Math.floor((y + rad) / c); j++)
        for (let k = Math.floor((z - rad) / c); k <= Math.floor((z + rad) / c); k++) {
          const a = this.m.get(this.key(i, j, k));
          if (a) for (const id of a) fn(id);
        }
  }
}

/** Random sequential placement of non-overlapping disks inside a circle of radius R. */
function placeDisks(radii, R, rng, pad = 0) {
  const order = radii.map((r, i) => i).sort((a, b) => radii[b] - radii[a]);
  const rmax = radii[order[0]];
  const h = new Hash2(2 * rmax + pad);
  const x = new Float64Array(radii.length), y = new Float64Array(radii.length);
  for (const i of order) {
    const r = radii[i];
    let px = 0, py = 0;
    for (let tries = 0; tries < 90; tries++) {
      const a = rng() * 2 * Math.PI, d = Math.sqrt(rng()) * Math.max(R - r, 0);
      px = Math.cos(a) * d; py = Math.sin(a) * d;
      let ok = true;
      h.each(px, py, r + rmax + pad, (j) => {
        if (ok && Math.hypot(px - x[j], py - y[j]) < r + radii[j] + pad) ok = false;
      });
      if (ok) break;
    }
    x[i] = px; y[i] = py;
    h.add(px, py, i);
  }
  return { x, y };
}

function circleOverlap(r1, r2, d) {
  if (d >= r1 + r2) return 0;
  if (d <= Math.abs(r1 - r2)) return Math.PI * Math.min(r1, r2) ** 2;
  const a = r1 * r1 * Math.acos((d * d + r1 * r1 - r2 * r2) / (2 * d * r1));
  const b = r2 * r2 * Math.acos((d * d + r2 * r2 - r1 * r1) / (2 * d * r2));
  const c = 0.5 * Math.sqrt((-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (d + r1 + r2));
  return a + b - c;
}

/* ───────────────────────── network builder ───────────────────────── */
class Builder {
  constructor() {
    this.nx = []; this.ny = []; this.nz = [];
    this.a = []; this.b = []; this.g = []; this.v = []; this.e = []; this.k = []; this.r = []; this.w = [];
    this.inlet = []; this.outlet = [];
  }
  node(x, y, z) { this.nx.push(x); this.ny.push(y); this.nz.push(z); return this.nx.length - 1; }
  dist(i, j) { return Math.hypot(this.nx[i] - this.nx[j], this.ny[i] - this.ny[j], this.nz[i] - this.nz[j]); }
  bond(a, b, g, v, e, kind, r, w = r) {
    this.a.push(a); this.b.push(b); this.g.push(g); this.v.push(v); this.e.push(e);
    this.k.push(kind); this.r.push(r); this.w.push(w);
  }
  /** Poiseuille tube; `len` overrides the node-to-node distance. */
  tube(a, b, r, kind, len) {
    const l = Math.max(len ?? this.dist(a, b), 1e-9);
    this.bond(a, b, (Math.PI * r ** 4) / (8 * l), Math.PI * r * r * l, r, kind, r);
  }
  finish(meta) {
    const nN = this.nx.length, nB = this.a.length;
    const a = Int32Array.from(this.a), b = Int32Array.from(this.b);
    const deg = new Int32Array(nN + 1);
    for (let i = 0; i < nB; i++) { deg[a[i] + 1]++; deg[b[i] + 1]++; }
    for (let i = 0; i < nN; i++) deg[i + 1] += deg[i];
    const adjBond = new Int32Array(2 * nB), fill = deg.slice(0, nN);
    for (let i = 0; i < nB; i++) { adjBond[fill[a[i]]++] = i; adjBond[fill[b[i]]++] = i; }
    const v = Float64Array.from(this.v);
    let vol = 0;
    for (let i = 0; i < nB; i++) vol += v[i];
    return {
      nNodes: nN, nBonds: nB,
      x: Float64Array.from(this.nx), y: Float64Array.from(this.ny), z: Float64Array.from(this.nz),
      a, b, g: Float64Array.from(this.g), v, e: Float64Array.from(this.e),
      kind: Uint8Array.from(this.k), r: Float32Array.from(this.r), w: Float32Array.from(this.w),
      adjStart: deg, adjBond,
      inlet: Int32Array.from(this.inlet), outlet: Int32Array.from(this.outlet),
      meta: { ...meta, poreVolume: vol, phi: vol / meta.Vbulk },
    };
  }
}

/* ───────────────────────── generators ───────────────────────── */
function genBundle(p, geo) {
  const rng = mulberry32(p.seed * 7 + 1);
  const draw = makeRadius(p, rng);
  const { R, L } = geo;
  const radii = Array.from({ length: p.N }, draw);
  const { x, y } = placeDisks(radii, R, rng);
  const B = new Builder();
  for (let i = 0; i < radii.length; i++) {
    const top = B.node(x[i], y[i], L / 2), bot = B.node(x[i], y[i], -L / 2);
    B.tube(top, bot, radii[i], KIND.TUBE);
    B.inlet.push(top); B.outlet.push(bot);
  }
  const spacing = R * Math.sqrt(Math.PI / p.N);
  return B.finish({ model: "bundle", R, L, Lv: L, A: Math.PI * R * R, Rh: R, band: 1.5 * spacing, Vbulk: Math.PI * R * R * L });
}

function genLayered(p, geo, withGap) {
  // radii get their own stream so resizing the sample (which changes placement) keeps the same pores
  const rng = mulberry32(p.seed * 7 + (withGap ? 3 : 2));
  const draw = makeRadius(p, mulberry32(p.seed * 7 + 11));
  const { R, t } = geo;
  const b0 = withGap ? p.bRatio * p.d50 * voidScale(p) : 0;
  const L = p.layers * t + (p.layers - 1) * b0;
  const spacing = R * Math.sqrt(Math.PI / p.N);
  const B = new Builder();
  const layers = [];
  let zTop = L / 2;
  let rmaxSeen = 0;

  for (let k = 0; k < p.layers; k++) {
    const radii = Array.from({ length: p.N }, draw);
    const { x, y } = placeDisks(radii, R, rng);
    const top = new Int32Array(p.N), mid = new Int32Array(p.N), bot = new Int32Array(p.N);
    const hash = new Hash2(spacing);
    for (let i = 0; i < p.N; i++) {
      top[i] = B.node(x[i], y[i], zTop);
      mid[i] = B.node(x[i], y[i], zTop - t / 2);
      bot[i] = B.node(x[i], y[i], zTop - t);
      B.tube(top[i], mid[i], radii[i], KIND.TUBE);
      B.tube(mid[i], bot[i], radii[i], KIND.TUBE);
      hash.add(x[i], y[i], i);
      if (radii[i] > rmaxSeen) rmaxSeen = radii[i];
      if (k === 0) B.inlet.push(top[i]);
      if (k === p.layers - 1) B.outlet.push(bot[i]);
    }
    // idea 2a: horizontal holes drilled through the glass between neighbouring holes
    if (p.lateral > 0) {
      const linked = new Set();
      for (let i = 0; i < p.N; i++) {
        const n = poisson(rng, p.lateral);
        if (!n) continue;
        const near = [];
        hash.each(x[i], y[i], 2.5 * spacing, (j) => {
          if (j !== i && Math.hypot(x[i] - x[j], y[i] - y[j]) < 2.5 * spacing) near.push(j);
        });
        for (let m = 0; m < n && near.length; m++) {
          const j = near.splice(Math.floor(rng() * near.length), 1)[0];
          const key = i < j ? i * 100000 + j : j * 100000 + i;
          if (linked.has(key)) continue;
          linked.add(key);
          const rl = Math.min(draw(), 0.45 * t);
          const len = Math.max(Math.hypot(x[i] - x[j], y[i] - y[j]) - radii[i] - radii[j], rl);
          B.tube(mid[i], mid[j], rl, KIND.LATERAL, len);
        }
      }
    }
    layers.push({ x, y, radii, top, bot, hash, zBot: zTop - t });
    zTop -= t + b0;
  }

  // idea 1: holes in touching plates connect only where their circles overlap
  for (let k = 0; k + 1 < layers.length; k++) {
    const A = layers[k], Bl = layers[k + 1];
    for (let i = 0; i < p.N; i++) {
      Bl.hash.each(A.x[i], A.y[i], A.radii[i] + rmaxSeen, (j) => {
        const d = Math.hypot(A.x[i] - Bl.x[j], A.y[i] - Bl.y[j]);
        const ov = circleOverlap(A.radii[i], Bl.radii[j], d);
        if (ov <= 0) return;
        const ro = Math.sqrt(ov / Math.PI);
        // Sampson orifice: Q = ΔP r³ / 3μ
        B.bond(A.bot[i], Bl.top[j], (ro ** 3) / 3, 0, ro, KIND.ORIFICE, ro);
      });
    }
  }

  // Guido's idea: a thin, rough gap between plates that every hole drains into
  if (withGap) {
    const c = spacing;
    const n = Math.ceil(R / c);
    for (let k = 0; k + 1 < layers.length; k++) {
      const z = layers[k].zBot - b0 / 2;
      const id = new Map();
      const key = (i, j) => (i + 5000) * 10007 + (j + 5000);
      for (let i = -n; i <= n; i++)
        for (let j = -n; j <= n; j++)
          if (c * Math.hypot(i, j) <= R - 0.3 * c) id.set(key(i, j), B.node(i * c, j * c, z));
      for (let i = -n; i <= n; i++)
        for (let j = -n; j <= n; j++) {
          const here = id.get(key(i, j));
          if (here === undefined) continue;
          for (const [di, dj] of [[1, 0], [0, 1]]) {
            const there = id.get(key(i + di, j + dj));
            if (there === undefined) continue;
            const bb = b0 * Math.exp(p.bSigma * normal(rng) - 0.5 * p.bSigma * p.bSigma);
            // parallel-plate cell: q = b³/12μ per unit width; width = length = c
            B.bond(here, there, (bb ** 3) / 12, bb * c * c * 0.5, bb, KIND.GAP, bb, c * 0.62);
          }
        }
      const nearest = (x, y) => {
        let i = Math.round(x / c), j = Math.round(y / c);
        let best = id.get(key(i, j));
        if (best !== undefined) return best;
        let bd = Infinity;
        for (let di = -2; di <= 2; di++)
          for (let dj = -2; dj <= 2; dj++) {
            const q = id.get(key(i + di, j + dj));
            if (q === undefined) continue;
            const d = Math.hypot(B.nx[q] - x, B.ny[q] - y);
            if (d < bd) { bd = d; best = q; }
          }
        return best;
      };
      const connect = (L1, nodes) => {
        for (let i = 0; i < p.N; i++) {
          const q = nearest(L1.x[i], L1.y[i]);
          if (q === undefined) continue;
          const r = L1.radii[i];
          // radial spreading from the hole edge into the gap cell
          const g = (2 * Math.PI * b0 ** 3) / (12 * Math.log(Math.max(c / (2 * r), 1.25)));
          B.bond(nodes[i], q, g, 0, b0, KIND.JUNCTION, b0);
        }
      };
      connect(layers[k], layers[k].bot);
      connect(layers[k + 1], layers[k + 1].top);
    }
  }

  return B.finish({
    model: withGap ? "gap" : "plates", R, L, Lv: L, A: Math.PI * R * R, Rh: R,
    band: 1.5 * spacing, Vbulk: Math.PI * R * R * L, t, b0, spacing,
  });
}

function rotate(v, axis, ang) {
  // Rodrigues rotation of v about a unit axis
  const [x, y, z] = v, [u, w, s] = axis, c = Math.cos(ang), sn = Math.sin(ang);
  const dot = x * u + y * w + z * s;
  return [
    x * c + (w * z - s * y) * sn + u * dot * (1 - c),
    y * c + (s * x - u * z) * sn + w * dot * (1 - c),
    z * c + (u * y - w * x) * sn + s * dot * (1 - c),
  ];
}
function perpendicular(v, rng) {
  const q = [normal(rng), normal(rng), normal(rng)];
  const c = [v[1] * q[2] - v[2] * q[1], v[2] * q[0] - v[0] * q[2], v[0] * q[1] - v[1] * q[0]];
  const n = Math.hypot(...c) || 1;
  return c.map((t) => t / n);
}
function unit(v) { const n = Math.hypot(...v) || 1; return v.map((t) => t / n); }

function growTree(p, geo, nTrunks) {
  const rng = mulberry32(p.seed * 7 + 4);
  const draw = makeRadius(p, rng);
  const { R, L } = geo;
  const s = p.d50;
  const B = new Builder();
  const owner = [];
  const tips = [];
  const addNode = (x, y, z, id) => { owner.push(id); return B.node(x, y, z); };

  const radii = Array.from({ length: nTrunks }, () => Math.min(draw() * p.trunkFactor, 0.15 * R));
  const { x, y } = placeDisks(radii, R, rng, 0.5 * s);
  const nSeg = Math.max(2, Math.round(L / s));
  const ds = L / nSeg;
  let branchId = 0;
  const ang = (p.branchAngle * Math.PI) / 180;

  const grow = (from, dir, r, gen, id) => {
    let cur = from;
    for (let sg = 0; sg < p.segsPerGen; sg++) {
      dir = unit(rotate(dir, perpendicular(dir, rng), 0.2 * normal(rng)));
      const len = s * 0.85 ** gen;
      let nx = B.nx[cur] + dir[0] * len, ny = B.ny[cur] + dir[1] * len, nz = B.nz[cur] + dir[2] * len;
      if (nz < -L / 2) {
        const f = (B.nz[cur] + L / 2) / (B.nz[cur] - nz);
        if (f < 0.05) break;
        nx = B.nx[cur] + (nx - B.nx[cur]) * f; ny = B.ny[cur] + (ny - B.ny[cur]) * f; nz = -L / 2;
        const n = addNode(nx, ny, nz, id);
        B.tube(cur, n, r, KIND.BRANCH);
        B.outlet.push(n);
        return;
      }
      if (nz > L / 2 || Math.hypot(nx, ny) > R - r) { tips.push([cur, r, id]); return; }
      const n = addNode(nx, ny, nz, id);
      B.tube(cur, n, r, KIND.BRANCH);
      cur = n;
    }
    if (gen < p.generations) {
      const axis = perpendicular(dir, rng);
      for (const sgn of [1, -1]) {
        const rd = r * 2 ** (-1 / 3) * Math.exp(0.15 * normal(rng)); // Murray's law
        grow(cur, unit(rotate(dir, axis, sgn * ang)), rd, gen + 1, id);
      }
    } else tips.push([cur, r, id]);
  };

  for (let i = 0; i < nTrunks; i++) {
    const own = -1 - i;
    let prev = addNode(x[i], y[i], L / 2, own);
    B.inlet.push(prev);
    for (let k = 1; k <= nSeg; k++) {
      const n = addNode(x[i], y[i], L / 2 - k * ds, own);
      B.tube(prev, n, radii[i], KIND.TRUNK);
      prev = n;
      if (k === nSeg) { B.outlet.push(n); break; }
      if (rng() < p.pBranch) {
        const az = rng() * 2 * Math.PI, dip = ((-40 + 50 * rng()) * Math.PI) / 180;
        const dir = [Math.cos(dip) * Math.cos(az), Math.cos(dip) * Math.sin(az), Math.sin(dip)];
        grow(n, dir, radii[i] * p.sideRatio, 0, branchId++);
      }
    }
  }

  // anastomosis: a branch tip that meets another branch fuses with it, closing a loop
  if (p.anastomosis > 0) {
    const h = new Hash3(1.6 * s);
    for (let n = 0; n < B.nx.length; n++) h.add(B.nx[n], B.ny[n], B.nz[n], n);
    for (const [tip, r, id] of tips) {
      if (rng() >= p.anastomosis) continue;
      let best = -1, bd = 1.6 * s;
      h.each(B.nx[tip], B.ny[tip], B.nz[tip], 1.6 * s, (n) => {
        if (owner[n] === id || n === tip) return;
        const d = B.dist(tip, n);
        if (d < bd && d > 1e-9) { bd = d; best = n; }
      });
      if (best >= 0) B.tube(tip, best, r, KIND.BRANCH);
    }
  }

  return B.finish({
    model: "tree", R, L, Lv: L, A: Math.PI * R * R, Rh: R, band: 1.5 * s,
    Vbulk: Math.PI * R * R * L, trunks: nTrunks, tips: tips.length,
  });
}

function genTree(p, geo) {
  return growTree(p, geo, Math.max(4, Math.round(p.N / 5.5)));
}

function genLattice(p, geo) {
  const rng = mulberry32(p.seed * 7 + 5);
  const draw = makeRadius(p, rng);
  const keep = Math.min(1, p.coord / 6);
  const a = geo.scale * Math.sqrt((3 * keep * Math.PI * geo.Er2) / p.phi);
  const nz = Math.max(3, Math.round(geo.L / a) + 1);
  const L = (nz - 1) * a;
  const n = Math.ceil(geo.R / a);
  const B = new Builder();
  const cols = [];
  for (let i = -n; i <= n; i++)
    for (let j = -n; j <= n; j++) if (a * Math.hypot(i, j) <= geo.R - 0.3 * a) cols.push([i, j]);
  const id = new Map();
  const key = (i, j, k) => ((i + 3000) * 6007 + (j + 3000)) * 6007 + k;
  for (const [i, j] of cols)
    for (let k = 0; k < nz; k++) {
      const q = B.node(i * a, j * a, L / 2 - k * a);
      id.set(key(i, j, k), q);
      if (k === 0) B.inlet.push(q);
      if (k === nz - 1) B.outlet.push(q);
    }
  for (const [i, j] of cols)
    for (let k = 0; k < nz; k++) {
      const here = id.get(key(i, j, k));
      for (const [di, dj, dk] of [[1, 0, 0], [0, 1, 0], [0, 0, 1]]) {
        const there = id.get(key(i + di, j + dj, k + dk));
        if (there === undefined || rng() >= keep) continue;
        B.tube(here, there, Math.min(draw(), 0.45 * a), KIND.THROAT);
      }
    }
  const A = cols.length * a * a;
  return B.finish({
    model: "lattice", R: geo.R, L, Lv: L, A, Rh: Math.sqrt(A / Math.PI), band: 1.5 * a,
    Vbulk: A * L, spacing: a,
  });
}

function buildOnce(p, model, geo) {
  switch (model) {
    case "bundle": return genBundle(p, geo);
    case "plates": return genLayered(p, geo, false);
    case "gap": return genLayered(p, geo, true);
    case "tree": return genTree(p, geo);
    case "lattice": return genLattice(p, geo);
    default: throw new Error(`unknown model ${model}`);
  }
}

/**
 * φ is the total porosity: holes, gaps, side holes and branches all count.
 * Pore sizes are fixed by the soil, so the sample is resized until the total matches.
 */
export function buildNetwork(p, model = p.model) {
  let scale = 1, net = null;
  for (let pass = 0; pass < 6; pass++) {
    net = buildOnce(p, model, geometry(p, scale));
    const ratio = net.meta.phi / p.phi;
    if (Math.abs(ratio - 1) < 0.01) break;
    scale *= Math.min(Math.max(Math.sqrt(ratio), 0.6), 1.8);
  }
  return net;
}

/* ───────────────────────── linear solver ───────────────────────── */
/**
 * IC(0) of A = diag − offdiag, with the off-diagonal magnitudes in (col, val).
 * Returns z ← (L Lᵀ)⁻¹ r.
 */
function incompleteCholesky(n, rowStart, col, val, diag) {
  // lower-triangular rows, sorted by column
  const lStart = new Int32Array(n + 1);
  for (let i = 0; i < n; i++) {
    let c = 0;
    for (let k = rowStart[i]; k < rowStart[i + 1]; k++) if (col[k] < i) c++;
    lStart[i + 1] = lStart[i] + c;
  }
  const lCol = new Int32Array(lStart[n]), lVal = new Float64Array(lStart[n]);
  for (let i = 0; i < n; i++) {
    const row = [];
    for (let k = rowStart[i]; k < rowStart[i + 1]; k++) if (col[k] < i) row.push([col[k], -val[k]]);
    row.sort((p, q) => p[0] - q[0]);
    for (let m = 0; m < row.length; m++) { lCol[lStart[i] + m] = row[m][0]; lVal[lStart[i] + m] = row[m][1]; }
  }
  const Ld = new Float64Array(n);
  const pos = new Int32Array(n).fill(-1);
  for (let i = 0; i < n; i++) {
    for (let k = lStart[i]; k < lStart[i + 1]; k++) pos[lCol[k]] = k;
    let sq = 0;
    for (let k = lStart[i]; k < lStart[i + 1]; k++) {
      const j = lCol[k];
      let s = lVal[k];
      for (let m = lStart[j]; m < lStart[j + 1]; m++) {
        const q = pos[lCol[m]];
        if (q >= 0 && lCol[m] < j) s -= lVal[q] * lVal[m];
      }
      lVal[k] = s / Ld[j];
      sq += lVal[k] * lVal[k];
    }
    const dd = diag[i] - sq;
    Ld[i] = Math.sqrt(dd > 1e-12 * diag[i] ? dd : diag[i]);
    for (let k = lStart[i]; k < lStart[i + 1]; k++) pos[lCol[k]] = -1;
  }
  return (r, z) => {
    for (let i = 0; i < n; i++) {
      let s = r[i];
      for (let k = lStart[i]; k < lStart[i + 1]; k++) s -= lVal[k] * z[lCol[k]];
      z[i] = s / Ld[i];
    }
    for (let i = n - 1; i >= 0; i--) {
      z[i] /= Ld[i];
      const zi = z[i];
      for (let k = lStart[i]; k < lStart[i + 1]; k++) z[lCol[k]] -= lVal[k] * zi;
    }
  };
}

/**
 * Solve Kirchhoff's law on the wet bonds with fixed potentials `fixed` (NaN = free).
 * `P` is a warm start and is overwritten with the solution (NaN where disconnected).
 * Returns the total current leaving the nodes held at potential 1.
 */
export function solveFlow(net, wet, fixed, P, g = net.g) {
  const { nNodes: N, adjStart, adjBond, a, b } = net;
  const act = new Uint8Array(N), queue = new Int32Array(N);
  let qh = 0, qt = 0;
  for (let n = 0; n < N; n++) if (!Number.isNaN(fixed[n])) { act[n] = 1; queue[qt++] = n; }
  while (qh < qt) {
    const n = queue[qh++];
    for (let k = adjStart[n]; k < adjStart[n + 1]; k++) {
      const bd = adjBond[k];
      if (!wet[bd]) continue;
      const m = a[bd] === n ? b[bd] : a[bd];
      if (!act[m]) { act[m] = 1; queue[qt++] = m; }
    }
  }
  const idx = new Int32Array(N).fill(-1);
  const free = [];
  for (let n = 0; n < N; n++) if (act[n] && Number.isNaN(fixed[n])) { idx[n] = free.length; free.push(n); }
  const nf = free.length;

  const rowStart = new Int32Array(nf + 1);
  for (let i = 0; i < nf; i++) {
    const n = free[i];
    let c = 0;
    for (let k = adjStart[n]; k < adjStart[n + 1]; k++) {
      const bd = adjBond[k];
      if (!wet[bd]) continue;
      const m = a[bd] === n ? b[bd] : a[bd];
      if (idx[m] >= 0) c++;
    }
    rowStart[i + 1] = rowStart[i] + c;
  }
  const col = new Int32Array(rowStart[nf]), val = new Float64Array(rowStart[nf]);
  const diag = new Float64Array(nf), rhs = new Float64Array(nf);
  for (let i = 0; i < nf; i++) {
    const n = free[i];
    let c = rowStart[i];
    for (let k = adjStart[n]; k < adjStart[n + 1]; k++) {
      const bd = adjBond[k];
      if (!wet[bd]) continue;
      const m = a[bd] === n ? b[bd] : a[bd];
      diag[i] += g[bd];
      if (idx[m] >= 0) { col[c] = idx[m]; val[c++] = g[bd]; }
      else rhs[i] += g[bd] * fixed[m];
    }
  }

  const x = new Float64Array(nf);
  for (let i = 0; i < nf; i++) { const w = P[free[i]]; x[i] = Number.isFinite(w) ? w : 0.5; }
  const mul = (v, out) => {
    for (let i = 0; i < nf; i++) {
      let s = diag[i] * v[i];
      for (let c = rowStart[i]; c < rowStart[i + 1]; c++) s -= val[c] * v[col[c]];
      out[i] = s;
    }
  };
  if (nf) {
    // Conductances span ~10 decades (r⁴, b³, orifice slivers), so plain Jacobi CG crawls.
    // Incomplete Cholesky IC(0) is stable for this M-matrix and exact on tree-like parts.
    const pre = incompleteCholesky(nf, rowStart, col, val, diag);
    const r = new Float64Array(nf), z = new Float64Array(nf), d = new Float64Array(nf), Ad = new Float64Array(nf);
    mul(x, Ad);
    let bn = 0;
    for (let i = 0; i < nf; i++) { r[i] = rhs[i] - Ad[i]; bn += rhs[i] * rhs[i]; }
    bn = Math.sqrt(bn) || 1;
    pre(r, z);
    d.set(z);
    let rz = 0;
    for (let i = 0; i < nf; i++) rz += r[i] * z[i];
    const maxIt = Math.max(200, 2 * nf);
    for (let it = 0; it < maxIt; it++) {
      mul(d, Ad);
      let dAd = 0;
      for (let i = 0; i < nf; i++) dAd += d[i] * Ad[i];
      if (!(dAd > 0)) break;
      const al = rz / dAd;
      let rn = 0;
      for (let i = 0; i < nf; i++) { x[i] += al * d[i]; r[i] -= al * Ad[i]; rn += r[i] * r[i]; }
      if (Math.sqrt(rn) < 1e-11 * bn) break;
      pre(r, z);
      let rz2 = 0;
      for (let i = 0; i < nf; i++) rz2 += r[i] * z[i];
      const be = rz2 / rz;
      rz = rz2;
      for (let i = 0; i < nf; i++) d[i] = z[i] + be * d[i];
    }
  }

  for (let n = 0; n < N; n++) P[n] = act[n] ? (idx[n] >= 0 ? x[idx[n]] : fixed[n]) : NaN;
  let Q = 0;
  for (let n = 0; n < N; n++) {
    if (fixed[n] !== 1) continue;
    for (let k = adjStart[n]; k < adjStart[n + 1]; k++) {
      const bd = adjBond[k];
      if (!wet[bd]) continue;
      const m = a[bd] === n ? b[bd] : a[bd];
      Q += g[bd] * (1 - P[m]);
    }
  }
  return Q;
}

/** Conductance per unit thickness of a unit-conductivity disk between two rim electrodes. */
export function diskShapeFactor(band, theta0 = (40 * Math.PI) / 180, n = 121) {
  const h = 2 / (n - 1);
  const id = new Int32Array(n * n).fill(-1);
  const pos = [];
  const fixed = [];
  for (let i = 0; i < n; i++)
    for (let j = 0; j < n; j++) {
      const X = -1 + i * h, Y = -1 + j * h, rho = Math.hypot(X, Y);
      if (rho > 1) continue;
      id[i * n + j] = pos.length;
      pos.push([i, j]);
      const ang = Math.atan2(Y, X);
      if (rho >= 1 - band && Math.abs(ang) <= theta0) fixed.push(0);
      else if (rho >= 1 - band && Math.PI - Math.abs(ang) <= theta0) fixed.push(1);
      else fixed.push(NaN);
    }
  const B = new Builder();
  for (const [i, j] of pos) B.node(i, j, 0);
  for (const [i, j] of pos) {
    const here = id[i * n + j];
    if (i + 1 < n && id[(i + 1) * n + j] >= 0) B.bond(here, id[(i + 1) * n + j], 1, 0, 1, 0, 1);
    if (j + 1 < n && id[i * n + j + 1] >= 0) B.bond(here, id[i * n + j + 1], 1, 0, 1, 0, 1);
  }
  const net = B.finish({ Vbulk: 1 });
  const wet = new Uint8Array(net.nBonds).fill(1);
  return solveFlow(net, wet, Float64Array.from(fixed), new Float64Array(net.nNodes).fill(0.5));
}

/* ───────────────────────── non-capillary film water ───────────────────────── */
// Clay and silt surfaces hold water as adsorbed films that no set of capillary holes can
// represent. Retention follows the PDI model (Peters 2013; Iden & Durner 2014): full up to
// the capillary air entry h_a, then falling log-linearly to zero at oven dryness (pF 6.8).
// Flow along those films follows Tuller & Or (2001): van der Waals film thickness
// δ = (A / 6πρgh)^⅓ and per-unit-width flux δ³/3μ. With wetted area θ_nc/δ per unit
// volume and random orientation (factor 1/3), K_film = ρg/9μ · θ_nc(h) · δ². No free parameter
// beyond θ_nc, which comes from retention.
export const HAMAKER = 6e-20;       // J, magnitude of the solid–water–air Hamaker constant
export const H_OVEN = 10 ** 7.8;    // mm, pF 6.8
export function filmSaturation(h, ha, b = 0.3) {
  if (!(h > 0) || h <= 0) return 1;
  const x = Math.log10(h), xa = Math.log10(Math.max(ha, 1)), x0 = Math.log10(H_OVEN);
  const v = 1 + (x - xa + b * Math.log1p(Math.exp((xa - x) / b))) / (xa - x0);
  return Math.min(1, Math.max(0, v));
}
export function filmThickness(h) {
  const hm = Math.max(h, 1) * 1e-3; // m
  return Math.cbrt(HAMAKER / (6 * Math.PI * FLUID.rho * FLUID.g * hm));
}
export function filmK(h, thetaFilm) {
  const d = filmThickness(h);
  return ((FLUID.rho * FLUID.g) / (9 * FLUID.mu)) * thetaFilm * d * d;
}

/* ───────────────────────── structural macropores ───────────────────────── */
// Fine soils are not packed grains: aggregates, root channels, worm burrows and shrinkage
// cracks add a coarse second pore system (dual porosity, Durner 1994). It is modelled as a
// bundle of straight, randomly oriented tubes with a lognormal radius distribution by volume
// (median r_m, spread σ_m), in parallel with the grain network. Retention: the macropores with
// r < r* = HEAD_CONST/h are still full. Flow: a tube at angle β to the gradient carries
// ρg r²/8μ · cos²β per unit tube area, and ⟨cos²β⟩ = 1/3 for random orientation, so
// K_ma(h) = ρg/24μ · w · E_v[r²; r < r*]. For a lognormal by volume that expectation is
// r_m² e^{2σ²} Φ((ln(r*/r_m) − 2σ²)/σ). No free parameter beyond w, r_m and σ_m, all from retention.
export const MACRO_SIGMA = 0.5;
export function normCdf(z) {
  // Abramowitz & Stegun 7.1.26, |error| < 1.5e-7
  const t = 1 / (1 + 0.3275911 * Math.abs(z) / Math.SQRT2);
  const y = 1 - t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429)))) * Math.exp(-z * z / 2);
  return z >= 0 ? 0.5 * (1 + y) : 0.5 * (1 - y);
}
/** Fraction of the macropore volume still full at head h (mm); rm in mm. */
export function macroSaturation(h, rm, sig = MACRO_SIGMA) {
  if (!(h > 0)) return 1;
  return normCdf(Math.log(HEAD_CONST / h / rm) / sig);
}
/** Macropore conductivity (m/s) at head h (mm) for volume fraction w and median radius rm (mm). */
export function macroK(h, w, rm, sig = MACRO_SIGMA) {
  const r = rm * 1e-3, s2 = sig * sig;
  const f = h > 0 ? normCdf((Math.log(HEAD_CONST / h / rm) - 2 * s2) / sig) : 1;
  return ((FLUID.rho * FLUID.g) / (24 * FLUID.mu)) * w * r * r * Math.exp(2 * s2) * f;
}

/* ───────────────────────── drainage and conductivity ───────────────────────── */
function invade(net, rstar, access, out) {
  out.fill(0);
  const { e, nBonds } = net;
  if (!access) {
    for (let i = 0; i < nBonds; i++) if (e[i] > rstar) out[i] = 1;
    return;
  }
  const { nNodes, adjStart, adjBond, a, b, inlet } = net;
  const seen = new Uint8Array(nNodes), q = new Int32Array(nNodes);
  let qh = 0, qt = 0;
  for (const n of inlet) if (!seen[n]) { seen[n] = 1; q[qt++] = n; }
  while (qh < qt) {
    const n = q[qh++];
    for (let k = adjStart[n]; k < adjStart[n + 1]; k++) {
      const bd = adjBond[k];
      if (out[bd] || e[bd] <= rstar) continue;
      out[bd] = 1;
      const m = a[bd] === n ? b[bd] : a[bd];
      if (!seen[m]) { seen[m] = 1; q[qt++] = m; }
    }
  }
}

export function boundarySets(net) {
  const { nNodes, x, y, meta } = net;
  const fixedV = new Float64Array(nNodes).fill(NaN);
  for (const n of net.inlet) fixedV[n] = 1;
  for (const n of net.outlet) fixedV[n] = 0;
  const fixedH = new Float64Array(nNodes).fill(NaN);
  const th = (40 * Math.PI) / 180;
  const edge = meta.Rh - meta.band;
  for (let n = 0; n < nNodes; n++) {
    if (Math.hypot(x[n], y[n]) < edge) continue;
    const ang = Math.atan2(y[n], x[n]);
    if (Math.abs(ang) <= th) fixedH[n] = 0;
    else if (Math.PI - Math.abs(ang) <= th) fixedH[n] = 1;
  }
  return { fixedV, fixedH };
}

/** Mualem's prediction from the network's own retention curve, K_r = S^l · M(S)². */
export function mualemFromRetention(h, S, l = 0.5) {
  // h ascending (h[0] = 0 at saturation), S descending
  const n = h.length;
  const cum = new Float64Array(n); // ∫ from S_i up to 1 of dS / h
  for (let i = 1; i < n; i++) {
    const h0 = h[i - 1] > 0 ? h[i - 1] : h[i];
    cum[i] = cum[i - 1] + (S[i - 1] - S[i]) * 0.5 * (1 / h0 + 1 / h[i]);
  }
  const tail = S[n - 1] / (2 * h[n - 1]); // what is left drains beyond the last head
  const total = cum[n - 1] + tail;
  const M = new Float64Array(n), Kr = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    M[i] = (total - cum[i]) / total;
    Kr[i] = S[i] > 0 ? S[i] ** l * M[i] ** 2 : 0;
  }
  return { M, Kr };
}

export function computeCurves(net, p, { keepFields = true } = {}) {
  const B = net.nBonds;
  const { meta } = net;
  // drainage thresholds: largest entry radius down to the radius that empties 99 % of the volume
  const order = Array.from({ length: B }, (_, i) => i).filter((i) => net.v[i] > 0).sort((i, j) => net.e[j] - net.e[i]);
  let acc = 0, rLow = net.e[order[order.length - 1]];
  for (const i of order) { acc += net.v[i]; if (acc >= 0.995 * meta.poreVolume) { rLow = net.e[i]; break; } }
  const rTop = net.e[order[0]];
  const nS = p.steps;
  const rstar = [Infinity];
  for (let k = 0; k < nS; k++) rstar.push(rTop * 0.9999 * (rLow / rTop) ** (k / (nS - 1)));

  const { fixedV, fixedH } = boundarySets(net);
  const shape = diskShapeFactor(meta.band / meta.Rh);
  const PV = new Float64Array(net.nNodes).fill(0.5), PH = new Float64Array(net.nNodes).fill(0.5);
  const inv = new Uint8Array(B), wet = new Uint8Array(B);
  // Angular pores: a drained pore keeps water in its corners (Tuller & Or 2001). For a triangle
  // the corner area is CORNER_AREA·(r*/r)² of the pore and its conductance scales as (r*/r)⁴.
  const kappa = p.corner || 0, cExp = p.cornerExp ?? 4;
  const gEff = new Float64Array(B), all = new Uint8Array(B).fill(1);
  const steps = [];
  let Qv0 = 0, Qh0 = 0;
  const masks = keepFields ? new Uint8Array(rstar.length * B) : null;
  const fluxV = keepFields ? new Float32Array(rstar.length * B) : null;
  const fluxH = keepFields ? new Float32Array(rstar.length * B) : null;

  for (let s = 0; s < rstar.length; s++) {
    invade(net, rstar[s], p.access, inv);
    let water = 0;
    for (let i = 0; i < B; i++) {
      wet[i] = inv[i] ? 0 : 1;
      if (wet[i]) { water += net.v[i]; gEff[i] = net.g[i]; continue; }
      const q = kappa > 0 ? Math.min(1, rstar[s] / net.e[i]) : 0;
      gEff[i] = kappa * net.g[i] * q ** cExp;
      water += CORNER_AREA * q * q * net.v[i];
    }
    const mask = kappa > 0 ? all : wet;
    // a cluster touching only one boundary carries no current; round-off makes it ±1e-13
    let Qv = solveFlow(net, mask, fixedV, PV, gEff);
    let Qh = solveFlow(net, mask, fixedH, PH, gEff);
    if (s === 0) { Qv0 = Qv; Qh0 = Qh; }
    if (!(Qv > 1e-12 * Qv0)) Qv = 0;
    if (!(Qh > 1e-12 * Qh0)) Qh = 0;
    const theta = water / meta.Vbulk;
    steps.push({
      rstar: rstar[s],
      h: Number.isFinite(rstar[s]) ? HEAD_CONST / rstar[s] : 0,
      theta,
      S: water / meta.poreVolume,
      Kv: (K_UNIT * Qv * meta.Lv) / meta.A,
      Kh: (K_UNIT * Qh) / (shape * meta.L),
    });
    if (keepFields) {
      masks.set(wet, s * B);
      const o = s * B;
      for (let i = 0; i < B; i++) {
        if (!(gEff[i] > 0)) continue;
        const pa = PV[net.a[i]], pb = PV[net.b[i]];
        if (Number.isFinite(pa) && Number.isFinite(pb)) fluxV[o + i] = gEff[i] * (pa - pb);
        const ha = PH[net.a[i]], hb = PH[net.b[i]];
        if (Number.isFinite(ha) && Number.isFinite(hb)) fluxH[o + i] = gEff[i] * (ha - hb);
      }
    }
  }

  if (p.film > 0) {
    const ha = steps[1].h, tot = meta.poreVolume / meta.Vbulk + p.film;
    for (const d of steps) {
      const tf = p.film * filmSaturation(d.h, ha);
      d.theta += tf;
      d.S = d.theta / tot;
      const kf = filmK(d.h, tf);
      d.Kv += kf; d.Kh += kf;
    }
  }
  const h = steps.map((d) => d.h), S = steps.map((d) => d.S);
  const mu = mualemFromRetention(h, S);
  const Kv0 = steps[0].Kv, Kh0 = steps[0].Kh;
  // fit Mualem's exponent l to the network: ln K_r − 2 ln M = l ln S
  let sxy = 0, sxx = 0;
  for (let i = 1; i < steps.length; i++) {
    const kr = steps[i].Kv / Kv0;
    if (!(kr > 0) || S[i] < 0.03 || S[i] > 0.97 || mu.M[i] <= 0) continue;
    const X = Math.log(S[i]), Y = Math.log(kr) - 2 * Math.log(mu.M[i]);
    sxy += X * Y; sxx += X * X;
  }
  let Sc = 0;
  for (let i = 1; i < steps.length; i++) if (!(steps[i].Kv > 0)) { Sc = steps[i - 1].S; break; }
  // how far the network sits from Mualem ℓ = 0.5, in decades, where both are above 10⁻⁶
  let se = 0, ne = 0;
  for (let i = 1; i < steps.length; i++) {
    const kr = steps[i].Kv / Kv0;
    if (S[i] < 0.05 || mu.Kr[i] < 1e-6) continue;
    se += (Math.log10(Math.max(kr, 1e-6)) - Math.log10(mu.Kr[i])) ** 2; ne++;
  }
  const d = p.d50 * 1e-3;
  const KC = ((FLUID.rho * FLUID.g) / FLUID.mu) * (d * d / 180) * (meta.phi ** 3 / (1 - meta.phi) ** 2);
  return {
    steps, masks, fluxV, fluxH,
    mualem: Array.from(mu.Kr),
    summary: {
      model: meta.model, phi: meta.phi + (p.film || 0), Ksat: Kv0, KsatH: Kh0, aniso: Kv0 > 0 ? Kh0 / Kv0 : NaN,
      lFit: sxx > 0 ? sxy / sxx : NaN, Sc, nodes: net.nNodes, bonds: net.nBonds,
      hAE: steps[1].h, shape, KC, mualemErr: ne ? Math.sqrt(se / ne) : NaN,
    },
  };
}
