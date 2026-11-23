import { MODELS, buildNetwork, computeCurves } from "./model.js";

function pack(net) {
  return {
    nNodes: net.nNodes, nBonds: net.nBonds,
    x: Float32Array.from(net.x), y: Float32Array.from(net.y), z: Float32Array.from(net.z),
    a: net.a, b: net.b, kind: net.kind, r: net.r, w: net.w, meta: net.meta,
  };
}

function curvesOnly(out) {
  return { steps: out.steps, mualem: out.mualem, summary: out.summary };
}

self.onmessage = (ev) => {
  const { id, type, params } = ev.data;
  try {
    if (type === "run") {
      const t0 = performance.now();
      const net = buildNetwork(params);
      const out = computeCurves(net, params, { keepFields: true });
      const v = pack(net);
      self.postMessage(
        { id, type, net: v, ...curvesOnly(out), masks: out.masks, fluxV: out.fluxV, fluxH: out.fluxH, ms: performance.now() - t0 },
        [v.x.buffer, v.y.buffer, v.z.buffer, out.masks.buffer, out.fluxV.buffer, out.fluxH.buffer],
      );
    } else if (type === "compare") {
      const results = {};
      for (const m of MODELS) {
        self.postMessage({ id, type: "progress", model: m });
        results[m] = curvesOnly(computeCurves(buildNetwork(params, m), params, { keepFields: false }));
      }
      self.postMessage({ id, type, results });
    }
  } catch (err) {
    self.postMessage({ id, type: "error", message: String(err && err.stack || err) });
  }
};
