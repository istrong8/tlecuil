const fs = require("fs");
const html = fs.readFileSync(process.argv[2], "utf8");
const src = html.split("// ENGINE-START")[1].split("// ENGINE-END")[0];
const m = { exports: {} }; new Function("module", src)(m); const FSP = m.exports;
const cases = JSON.parse(fs.readFileSync(process.argv[3], "utf8"));
let worst = 0;
for (const [name, c] of Object.entries(cases)) {
  const r = FSP.size(c.inputs);
  for (const [k, v] of Object.entries(c.results)) {
    const rel = Math.abs(r[k] - v) / Math.max(Math.abs(v), 1e-12);
    if (Math.abs(v) > 1e-9) worst = Math.max(worst, rel);
  }
  if (c.optimum_tc_k) { const [tc] = FSP.optimumTc(c.inputs); worst = Math.max(worst, Math.abs(tc - c.optimum_tc_k) / c.optimum_tc_k); }
  console.log(name.padEnd(28), "Pth", r.reactor_thermal_kwt.toFixed(3), "A", r.radiator_area_m2.toFixed(3));
}
console.log("worst relative difference vs Python:", worst.toExponential(2));
process.exit(worst < 1e-9 ? 0 : 1);
