import test from "node:test";
import assert from "node:assert/strict";
import * as M from "../profiles/model.js";

const p = M.soil(M.DEFAULTS);
const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg}: ${a} vs ${b}`);

test("c and λc from the notes", () => {
  close(p.c, 3 + 2 / 1.2, 1e-12, "c");
  close(p.eta, 5.6, 1e-12, "λc");
});

test("hydrostatic profile is s = z, so Θ_E = (z/|h_A|)^−λ above the fringe", () => {
  const pr = M.profile(p, 0);
  for (const z of [10, 30, 45, 60, 80]) {
    const exact = z <= 30 ? 1 : Math.pow(z / 30, -1.2);
    close(M.thetaAt(pr, z), exact, 2e-3, `z = ${z}`);
  }
  close(pr.thetaSurf, Math.pow(80 / 30, -1.2), 1e-5, "surface");
});

test("tension-saturated zone is |h_A| / (1 + q/K_s), capped at Θ_E = 1", () => {
  close(M.profile(p, -0.5).zSat, 60, 1e-9, "−½ K_s");
  close(M.profile(p, -0.25).zSat, 40, 1e-9, "−¼ K_s");
  for (const a of [-0.5, -0.25, 0, 1e-3]) {
    for (const pt of M.profile(p, a).pts) assert.ok(pt.theta <= 1 && pt.theta >= 0);
  }
});

test("percolation relaxes toward K(Θ*) = |q|, so Θ* = (|q|/K_s)^(1/c)", () => {
  for (const a of [-0.5, -0.25]) {
    const pr = M.profile(p, a);
    const thetaStar = Math.pow(-a, 1 / p.c);
    close(M.thetaE(p, pr.sInf), thetaStar, 1e-12, "asymptote");
    assert.ok(pr.thetaSurf > thetaStar && pr.thetaSurf - thetaStar < 0.01, `a = ${a} surface ${pr.thetaSurf}`);
  }
});

test("RK4 in z matches the quadrature z(s) = z_sat + ∫ ds / (1 + a (s/s_A)^η)", () => {
  const a = -0.25;
  const pr = M.profile(p, a);
  const pt = pr.pts.find((q) => q.z > 55);
  // Simpson on s from s_A to pt.s
  const n = 2000, h = (pt.s - p.sA) / n;
  const f = (s) => 1 / (1 + a * Math.pow(s / p.sA, p.eta));
  let sum = f(p.sA) + f(pt.s);
  for (let i = 1; i < n; i++) sum += (i % 2 ? 4 : 2) * f(p.sA + i * h);
  close(pr.zSat + (sum * h) / 3, pt.z, 1e-4, "height");
});

test("q_max: suction reaches ∞ exactly at the surface, about 1.33 K(L)", () => {
  const qm = M.qMax(p);
  close(M.dryHeight(p, qm), p.L, 1e-6, "dry height");
  close(qm / M.qRef(p), 1.3327, 1e-3, "ratio");
  assert.equal(M.profile(p, qm * 0.999).dryAt, null);
  assert.ok(M.profile(p, qm * 1.05).dryAt < p.L);
});

test("evaporation profiles from the notes stay wet at the surface", () => {
  const r = M.solveAll(M.DEFAULTS);
  const evap = r.profiles.find((x) => x.key === "evap-97");
  close(evap.thetaSurf, 0.223, 2e-3, "Θ_E at surface for .97 K(L)");
  // ordering: more upward flux → drier at every height
  const surf = r.profiles.map((x) => x.thetaSurf);
  for (let i = 1; i < surf.length; i++) assert.ok(surf[i] < surf[i - 1]);
});
