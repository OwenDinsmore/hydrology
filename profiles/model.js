/* Steady vertical flow above a water table, Brooks–Corey soil.
 *
 * z is height above the water table (cm, up is positive), s = −h is suction (cm),
 * q is the Darcy flux, positive upward (capillary rise), negative downward (percolation).
 *
 *   Θ_E(s) = 1                    s ≤ |h_A|
 *          = (s/|h_A|)^(−λ)        s > |h_A|
 *   K(s)   = K_s Θ_E^c,  c = 3 + 2/λ   →   K/K_s = (s/|h_A|)^(−λc)
 *
 *   q = −K (dh/dz + 1)   →   ds/dz = 1 + q/K(s)
 *
 * Every flux in the assignment is a multiple of K_s, so the profiles only need a = q/K_s
 * and K_s drops out of Θ_E(z).
 */

export const DEFAULTS = { hA: -30, lam: 1.2, Ks: 3.4e-4, L: 80 };

// Clapp and Hornberger (1978), Table 2: ψ_ae (cm), b = 1/λ, K_s (cm/s). Their K exponent is
// 2b + 3 = 3 + 2/λ, the same c as the assignment, so the classes drop straight into this model.
export const TEXTURES = [
  { key: "sand", name: "Sand", hA: -12.1, b: 4.05, Ks: 1.76e-2 },
  { key: "loamy-sand", name: "Loamy sand", hA: -9.0, b: 4.38, Ks: 1.56e-2 },
  { key: "sandy-loam", name: "Sandy loam", hA: -21.8, b: 4.90, Ks: 3.41e-3 },
  { key: "silt-loam", name: "Silt loam", hA: -78.6, b: 5.30, Ks: 7.2e-4 },
  { key: "loam", name: "Loam", hA: -47.8, b: 5.39, Ks: 6.95e-4 },
  { key: "sandy-clay-loam", name: "Sandy clay loam", hA: -29.9, b: 7.12, Ks: 6.30e-4 },
  { key: "silty-clay-loam", name: "Silty clay loam", hA: -35.6, b: 7.75, Ks: 1.70e-4 },
  { key: "clay-loam", name: "Clay loam", hA: -63.0, b: 8.52, Ks: 2.45e-4 },
  { key: "sandy-clay", name: "Sandy clay", hA: -15.3, b: 10.4, Ks: 2.17e-4 },
  { key: "silty-clay", name: "Silty clay", hA: -49.0, b: 10.4, Ks: 1.03e-4 },
  { key: "clay", name: "Clay", hA: -40.5, b: 11.4, Ks: 1.28e-4 },
].map((t) => ({ ...t, lam: 1 / t.b }));

/* The five fluxes from the assignment, as fractions of K_s. qRef is K at the surface under
 * hydrostatic conditions, K_s (L/|h_A|)^(−λc). */
export const CASES = [
  { key: "perc-half", label: "−½ K_s", kind: "percolation", frac: -0.5, ofRef: false },
  { key: "perc-quarter", label: "−¼ K_s", kind: "percolation", frac: -0.25, ofRef: false },
  { key: "hydro", label: "0", kind: "hydrostatic", frac: 0, ofRef: false },
  { key: "evap-half", label: "½ K(L)", kind: "capillary rise", frac: 0.5, ofRef: true },
  { key: "evap-97", label: ".97 K(L)", kind: "capillary rise", frac: 0.97, ofRef: true },
];

export function soil({ hA, lam, L }) {
  const sA = Math.abs(hA);
  const c = 3 + 2 / lam;
  return { sA, lam, c, eta: lam * c, L };
}

export const thetaE = (p, s) => (s <= p.sA ? 1 : Math.pow(s / p.sA, -p.lam));
export const kRel = (p, s) => (s <= p.sA ? 1 : Math.pow(s / p.sA, -p.eta));

/** K at the surface for a hydrostatic profile, over K_s. */
export const qRef = (p) => Math.pow(p.L / p.sA, -p.eta);

/** a = q/K_s for each assignment case. */
export function caseFlux(p, cs) {
  return cs.ofRef ? cs.frac * qRef(p) : cs.frac;
}

// Simpson on [0, 1] for ∫ dy / (1 + a y^η), which is smooth.
function lowerPart(a, eta, n = 400) {
  const f = (y) => 1 / (1 + a * Math.pow(y, eta));
  const h = 1 / n;
  let sum = f(0) + f(1);
  for (let i = 1; i < n; i++) sum += (i % 2 ? 4 : 2) * f(i * h);
  return (sum * h) / 3;
}

/** Height where suction goes to infinity for an upward flux a = q/K_s > 0. */
export function dryHeight(p, a) {
  // ∫_0^∞ dy / (1 + a y^η) = a^(−1/η) (π/η) / sin(π/η)
  const full = Math.pow(a, -1 / p.eta) * (Math.PI / p.eta) / Math.sin(Math.PI / p.eta);
  return p.sA / (1 + a) + p.sA * (full - lowerPart(a, p.eta));
}

/** Largest steady upward flux (over K_s) the column can lift to the surface. */
export function qMax(p) {
  let lo = -60, hi = 10; // ln a
  for (let i = 0; i < 200; i++) {
    const mid = (lo + hi) / 2;
    if (dryHeight(p, Math.exp(mid)) > p.L) lo = mid; else hi = mid;
  }
  return Math.exp((lo + hi) / 2);
}

/**
 * Θ_E(z) from the water table to the surface for a = q/K_s.
 * Returns points ordered by height plus summary numbers.
 */
export function profile(p, a) {
  const pts = [{ z: 0, s: 0, theta: 1 }];
  const out = { a, pts, zSat: null, sInf: null, dryAt: null, saturated: false };

  if (a <= -1) {
    // downward flux at or above K_s: friction carries all the weight, suction never builds
    out.saturated = true;
    out.zSat = p.L;
    pts.push({ z: p.L, s: 0, theta: 1 });
    return finish(p, out);
  }

  const zSat = p.sA / (1 + a);
  out.zSat = Math.min(zSat, p.L);
  if (zSat >= p.L) {
    out.saturated = true;
    pts.push({ z: p.L, s: p.sA * (p.L / zSat), theta: 1 });
    return finish(p, out);
  }
  pts.push({ z: zSat, s: p.sA, theta: 1 });

  if (a < 0) {
    // percolation: march up in z; s relaxes toward s*, where K(s*) = |q|
    out.sInf = p.sA * Math.pow(-1 / a, 1 / p.eta);
    const f = (s) => 1 + a * Math.pow(s / p.sA, p.eta);
    const dz = 0.05;
    let z = zSat, s = p.sA, k = 0;
    while (z < p.L - 1e-9) {
      const h = Math.min(dz, p.L - z);
      const k1 = f(s), k2 = f(s + (h / 2) * k1), k3 = f(s + (h / 2) * k2), k4 = f(s + h * k3);
      s += (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
      z += h;
      if (++k % 4 === 0 || z >= p.L - 1e-9) pts.push({ z, s, theta: thetaE(p, s) });
    }
    return finish(p, out);
  }

  // a ≥ 0: integrate z(s) = zSat + ∫ ds / (1 + a (s/s_A)^η) in u = ln(s/s_A), which handles s → ∞
  const g = (u) => (p.sA * Math.exp(u)) / (1 + a * Math.exp(p.eta * u));
  const du = 0.002;
  const uMax = Math.log(1e8) / p.lam; // Θ_E = 1e−8
  let z = zSat, u = 0, k = 0;
  while (u < uMax) {
    const dzStep = (du / 6) * (g(u) + 4 * g(u + du / 2) + g(u + du));
    if (z + dzStep >= p.L) {
      // land exactly on the surface
      const frac = (p.L - z) / dzStep;
      u += frac * du;
      const s = p.sA * Math.exp(u);
      pts.push({ z: p.L, s, theta: thetaE(p, s) });
      return finish(p, out);
    }
    z += dzStep;
    u += du;
    if (++k % 5 === 0) {
      const s = p.sA * Math.exp(u);
      pts.push({ z, s, theta: thetaE(p, s) });
    }
  }
  // the flux is more than the soil can lift: Θ_E → 0 below the surface, no steady state
  out.dryAt = a > 0 ? dryHeight(p, a) : z;
  pts.push({ z: out.dryAt, s: Infinity, theta: 0 });
  return finish(p, out);
}

function finish(p, out) {
  const top = out.pts.at(-1);
  out.sSurf = out.dryAt != null ? Infinity : top.s;
  out.thetaSurf = out.dryAt != null ? 0 : top.theta;
  return out;
}

/** Θ_E at a given height by linear interpolation along a profile. */
export function thetaAt(prof, z) {
  const pts = prof.pts;
  if (prof.dryAt != null && z >= prof.dryAt) return 0;
  for (let i = 1; i < pts.length; i++) {
    if (pts[i].z >= z) {
      const a = pts[i - 1], b = pts[i];
      const t = b.z === a.z ? 1 : (z - a.z) / (b.z - a.z);
      return a.theta + t * (b.theta - a.theta);
    }
  }
  return pts.at(-1).theta;
}

/** All five assignment profiles for one soil. */
export function solveAll(input) {
  const p = soil(input);
  const ref = qRef(p);
  const max = qMax(p);
  const profiles = CASES.map((cs) => ({ ...cs, ...profile(p, caseFlux(p, cs)) }));
  return { p, input, qRef: ref, qMax: max, profiles };
}
