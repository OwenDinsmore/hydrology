import { DEFAULTS, solveAll } from "./model.js";

const COLORS = {
  "perc-half": "#1f5fa8",
  "perc-quarter": "#6b9fd8",
  hydro: "#5c5f66",
  "evap-half": "#e08a4c",
  "evap-97": "#b23a12",
};
const NAMES = {
  "perc-half": "−½ K<sub>s</sub>",
  "perc-quarter": "−¼ K<sub>s</sub>",
  hydro: "0 (hydrostatic)",
  "evap-half": "½ K(L)",
  "evap-97": "0.97 K(L)",
};

const res = solveAll(DEFAULTS);
const L = res.p.L;

// plot frame
const W = 640, H = 440, m = { t: 16, r: 16, b: 46, l: 56 };
const x = (v) => m.l + v * (W - m.l - m.r);
const y = (z) => H - m.b - (z / L) * (H - m.t - m.b);

const ns = "http://www.w3.org/2000/svg";
const svg = document.getElementById("plot");
function add(name, attrs, text) {
  const e = document.createElementNS(ns, name);
  for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v);
  if (text != null) e.textContent = text;
  svg.appendChild(e);
  return e;
}

for (let v = 0; v <= 1.0001; v += 0.2) {
  add("line", { x1: x(v), x2: x(v), y1: y(0), y2: y(L), class: "grid" });
  add("text", { x: x(v), y: y(0) + 18, "text-anchor": "middle", class: "tick" }, v.toFixed(1));
}
for (let z = 0; z <= L; z += 10) {
  add("line", { x1: x(0), x2: x(1), y1: y(z), y2: y(z), class: "grid" });
  add("text", { x: x(0) - 8, y: y(z) + 4, "text-anchor": "end", class: "tick" }, z);
}
add("path", { d: `M${x(0)} ${y(L)} V${y(0)} H${x(1)}`, class: "axis" });
add("text", { x: (x(0) + x(1)) / 2, y: H - 8, "text-anchor": "middle", class: "label" }, "effective saturation ΘE");
add("text", { transform: `translate(16 ${(y(0) + y(L)) / 2}) rotate(-90)`, "text-anchor": "middle", class: "label" }, "height above water table (cm)");
add("text", { x: x(0) + 6, y: y(L) + 14, class: "tick" }, "surface");

for (const pr of res.profiles) {
  const d = pr.pts.map((p) => `${x(p.theta).toFixed(1)} ${y(p.z).toFixed(1)}`).join(" L");
  add("path", {
    d: "M" + d,
    fill: "none",
    stroke: COLORS[pr.key],
    "stroke-width": 2.2,
    "stroke-dasharray": pr.key === "hydro" ? "6 4" : "none",
  });
}

document.getElementById("legend").innerHTML = res.profiles
  .map((pr) => `<span><i style="border-color:${COLORS[pr.key]};border-top-style:${pr.key === "hydro" ? "dashed" : "solid"}"></i>q = ${NAMES[pr.key]}</span>`)
  .join("");

document.getElementById("summary").innerHTML = res.profiles
  .map((pr) => `<tr><td>${NAMES[pr.key]}</td><td>${pr.zSat.toFixed(1)} cm</td><td>${pr.thetaSurf.toFixed(2)}</td></tr>`)
  .join("");
