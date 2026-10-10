#!/usr/bin/env node
const fs = require('fs');
const path = require('path');
const { createCanvas } = require('C:/Users/ASUS/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/@napi-rs/canvas');

const root = path.resolve(__dirname, '..');
const runDirs = [
  path.join(root, 'runs/sign-size-lambda-L4-5seed-20260918/results'),
];
const outDir = path.join(root, 'results/sign-size-lambda-L4-5seed-20260918/plots');
fs.mkdirSync(outDir, { recursive: true });

const grouped = new Map();
for (const run of runDirs) {
  for (const point of fs.readdirSync(run).filter(x => x.startsWith('point_'))) {
    const dir = path.join(run, point);
    for (const file of fs.readdirSync(dir).filter(x => x.endsWith('_results.json'))) {
      const d = JSON.parse(fs.readFileSync(path.join(dir, file), 'utf8'));
      const lam = Number(d.lambda.toFixed(10));
      if (!grouped.has(lam)) grouped.set(lam, []);
      grouped.get(lam).push(d);
    }
  }
}
const lambdas = [...grouped.keys()].sort((a, b) => a - b);
for (const x of lambdas) {
  if (grouped.get(x).length < 4) throw new Error(`lambda=${x} has ${grouped.get(x).length} chains`);
}

function combine(items, getter) {
  const a = items.map(getter), n = a.length;
  const mean = a.reduce((s, x) => s + x[0], 0) / n;
  const se = Math.sqrt(a.reduce((s, x) => s + x[1] * x[1], 0)) / n;
  return [mean, se];
}
function derived(d, key) {
  const x = d.vbs.blocking[0].derived[key];
  return [x.value[0], x.stderr_real];
}
function raw(d, key, component) {
  const x = d.vbs.blocking[0].raw[key];
  return [x.value[component], component ? x.stderr_imag : x.stderr_real];
}
function series(getter) {
  return lambdas.map(x => [x, ...combine(grouped.get(x), getter)]);
}

const W = 1800, H = 1280;
const palette = ['#2563eb', '#ea580c', '#16a34a', '#dc2626'];
function niceTicks(lo, hi, n = 5) {
  const span = Math.max(hi - lo, 1e-12);
  const raw = span / n;
  const mag = 10 ** Math.floor(Math.log10(raw));
  const r = raw / mag;
  const step = (r < 1.5 ? 1 : r < 3 ? 2 : r < 7 ? 5 : 10) * mag;
  const start = Math.floor(lo / step) * step, end = Math.ceil(hi / step) * step;
  const vals = [];
  for (let v = start; v <= end + step * 0.1; v += step) vals.push(v);
  return vals;
}
function fmt(v) {
  if (Math.abs(v) >= 100 || (Math.abs(v) > 0 && Math.abs(v) < 0.001)) return v.toExponential(1);
  return Number(v.toFixed(3)).toString();
}
function drawPanel(ctx, box, data, title, ylabel, color) {
  const m = { l: 118, r: 28, t: 70, b: 82 };
  const x0 = box.x + m.l, x1 = box.x + box.w - m.r;
  const y0 = box.y + m.t, y1 = box.y + box.h - m.b;
  const xmin = Math.min(...data.map(d => d[0])), xmax = Math.max(...data.map(d => d[0]));
  let ymin = Math.min(...data.map(d => d[1] - d[2])), ymax = Math.max(...data.map(d => d[1] + d[2]));
  const pad = Math.max((ymax - ymin) * 0.1, 1e-5); ymin -= pad; ymax += pad;
  const sx = x => x0 + (x - xmin) / (xmax - xmin) * (x1 - x0);
  const sy = y => y1 - (y - ymin) / (ymax - ymin) * (y1 - y0);

  ctx.strokeStyle = '#d1d5db'; ctx.lineWidth = 1;
  const yt = niceTicks(ymin, ymax, 5);
  ctx.font = '22px Arial'; ctx.fillStyle = '#1f2937';
  for (const t of yt) {
    const yy = sy(t); if (yy < y0 - 2 || yy > y1 + 2) continue;
    ctx.beginPath(); ctx.moveTo(x0, yy); ctx.lineTo(x1, yy); ctx.stroke();
    ctx.textAlign = 'right'; ctx.textBaseline = 'middle'; ctx.fillText(fmt(t), x0 - 12, yy);
  }
  const xt = [0, 0.5, 1, 1.5, 2, 2.5];
  for (const t of xt) {
    const xx = sx(t); ctx.beginPath(); ctx.moveTo(xx, y0); ctx.lineTo(xx, y1); ctx.stroke();
    ctx.textAlign = 'center'; ctx.textBaseline = 'top'; ctx.fillText(fmt(t), xx, y1 + 12);
  }
  ctx.strokeStyle = '#374151'; ctx.lineWidth = 2; ctx.strokeRect(x0, y0, x1 - x0, y1 - y0);
  ctx.strokeStyle = color; ctx.lineWidth = 3; ctx.beginPath();
  data.forEach((d, i) => { const xx = sx(d[0]), yy = sy(d[1]); i ? ctx.lineTo(xx, yy) : ctx.moveTo(xx, yy); }); ctx.stroke();
  ctx.lineWidth = 2;
  for (const d of data) {
    const xx = sx(d[0]), yy = sy(d[1]), ya = sy(d[1] - d[2]), yb = sy(d[1] + d[2]);
    ctx.beginPath(); ctx.moveTo(xx, ya); ctx.lineTo(xx, yb); ctx.moveTo(xx - 7, ya); ctx.lineTo(xx + 7, ya); ctx.moveTo(xx - 7, yb); ctx.lineTo(xx + 7, yb); ctx.stroke();
    ctx.fillStyle = color; ctx.beginPath(); ctx.arc(xx, yy, 6, 0, Math.PI * 2); ctx.fill();
  }
  ctx.fillStyle = '#111827'; ctx.font = 'bold 28px Arial'; ctx.textAlign = 'center'; ctx.textBaseline = 'top'; ctx.fillText(title, (x0 + x1) / 2, box.y + 20);
  ctx.font = '24px Arial'; ctx.fillText('λ', (x0 + x1) / 2, y1 + 48);
  ctx.save(); ctx.translate(box.x + 24, (y0 + y1) / 2); ctx.rotate(-Math.PI / 2); ctx.textBaseline = 'top'; ctx.fillText(ylabel, 0, 0); ctx.restore();
}
function makeFigure(filename, heading, panels, footnote) {
  const canvas = createCanvas(W, H), ctx = canvas.getContext('2d');
  ctx.fillStyle = '#ffffff'; ctx.fillRect(0, 0, W, H);
  ctx.fillStyle = '#111827'; ctx.font = 'bold 38px Arial'; ctx.textAlign = 'center'; ctx.textBaseline = 'top'; ctx.fillText(heading, W / 2, 24);
  const gap = 24, top = 82, bottom = footnote ? 55 : 20;
  const bw = (W - gap * 3) / 2, bh = (H - top - bottom - gap) / 2;
  panels.forEach((p, i) => drawPanel(ctx, { x: gap + (i % 2) * (bw + gap), y: top + Math.floor(i / 2) * (bh + gap), w: bw, h: bh }, p.data, p.title, p.ylabel, palette[i]));
  if (footnote) { ctx.fillStyle = '#4b5563'; ctx.font = '20px Arial'; ctx.textAlign = 'center'; ctx.fillText(footnote, W / 2, H - 35); }
  fs.writeFileSync(path.join(outDir, filename), canvas.toBuffer('image/png'));
}

makeFigure('order_parameters_vs_lambda.png', 'L=4 order parameters versus λ (5-seed scan)', [
  { title: 'VBS m_z²(q1)', ylabel: 'm_z²(q1)', data: series(d => derived(d, 'mz2_q1')) },
  { title: 'VBS m_z²(q3)', ylabel: 'm_z²(q3)', data: series(d => derived(d, 'mz2_q3')) },
  { title: 'Leading eigenvalue at q1', ylabel: 'leading eigenvalue', data: series(d => derived(d, 'leading_eigenvalue_q1')) },
  { title: 'Uniform displacement X_z', ylabel: 'X_z', data: series(d => [d.observables.X_z.value[0], d.observables.X_z.stderr_real]) },
], '1600 samples per point; λ=2.5 uses 4 valid chains. X_z is q=0, not finite-q VBS order.');

const f00 = series(d => raw(d, 'F_bare_dx0_dy0_AB', 1));
makeFigure('fermion_correlations_vs_lambda.png', 'L=4 Majorana correlations versus λ (5-seed scan)', [
  { title: 'Im F_AB(0,0)', ylabel: 'Im F_AB(0,0)', data: f00 },
  { title: 'Im F_BA(0,1)', ylabel: 'Im F_BA(0,1)', data: series(d => raw(d, 'F_bare_dx0_dy1_BA', 1)) },
  { title: 'Im F_AB(1,1)', ylabel: 'Im F_AB(1,1)', data: series(d => raw(d, 'F_bare_dx1_dy1_AB', 1)) },
  { title: '|F_AB(0,0)|', ylabel: '|F_AB(0,0)|', data: f00.map(d => [d[0], Math.abs(d[1]), d[2]]) },
], '1600 samples per point except λ=2.5 (1280); center equal-time matter-Majorana, not spin correlations.');

console.log(path.join(outDir, 'order_parameters_vs_lambda.png'));
console.log(path.join(outDir, 'fermion_correlations_vs_lambda.png'));
