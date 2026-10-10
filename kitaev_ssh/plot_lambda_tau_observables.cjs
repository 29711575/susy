#!/usr/bin/env node
const fs = require('fs');
const path = require('path');
const { createCanvas } = require('C:/Users/ASUS/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/@napi-rs/canvas');

const root = path.resolve(__dirname, '..');
const runDir = path.join(root, 'runs/sign-size-lambda-L4-tau-3seed-20260918/results');
const outDir = path.join(root, 'results/sign-size-lambda-L4-tau-3seed-20260918/plots');
fs.mkdirSync(outDir, { recursive: true });

// Read all three independent chains at each (lambda, tau) point.
const groups = new Map();
for (const point of fs.readdirSync(runDir).filter(x => /^point_\d+$/.test(x))) {
  const dir = path.join(runDir, point);
  for (const file of fs.readdirSync(dir).filter(x => x.endsWith('_results.json'))) {
    const d = JSON.parse(fs.readFileSync(path.join(dir, file), 'utf8'));
    const lambda = Number(d.lambda);
    const tau = Number(d.tau);
    const key = `${lambda.toPrecision(12)}|${tau.toPrecision(12)}`;
    if (!groups.has(key)) groups.set(key, { lambda, tau, items: [] });
    groups.get(key).items.push(d);
  }
}
for (const g of groups.values()) {
  if (g.items.length !== 3) throw new Error(`lambda=${g.lambda}, tau=${g.tau}: expected 3 chains, found ${g.items.length}`);
}
const lambdas = [...new Set([...groups.values()].map(g => g.lambda))].sort((a, b) => a - b);
const taus = [...new Set([...groups.values()].map(g => g.tau))].sort((a, b) => a - b);
if (lambdas.length !== 15 || taus.length !== 6) {
  throw new Error(`unexpected scan grid: ${lambdas.length} lambdas x ${taus.length} tau values`);
}

function valueAndError(d, kind, key) {
  if (kind === 'phase') {
    return [d.average_phase.value[0], d.average_phase.stderr_real || 0];
  }
  if (kind === 'derived') {
    const x = d.vbs.blocking[0].derived[key];
    return [x.value[0], x.stderr_real || 0];
  }
  if (kind === 'observable') {
    const x = d.observables[key];
    return [x.value[0], x.stderr_real || 0];
  }
  if (kind === 'fermion-imag') {
    const x = d.vbs.blocking[0].raw[key];
    return [x.value[1], x.stderr_imag || 0];
  }
  throw new Error(`unknown metric kind: ${kind}`);
}

function combine(items, kind, key) {
  const vals = items.map(d => valueAndError(d, kind, key));
  const n = vals.length;
  const mean = vals.reduce((s, x) => s + x[0], 0) / n;
  // Independent-chain standard errors combined for the mean.
  const se = Math.sqrt(vals.reduce((s, x) => s + x[1] * x[1], 0)) / n;
  return { y: mean, e: se };
}

function curves(kind, key) {
  return lambdas.map(lambda => ({
    lambda,
    points: taus.map(tau => {
      const g = groups.get(`${lambda.toPrecision(12)}|${tau.toPrecision(12)}`);
      const z = combine(g.items, kind, key);
      return { x: tau, y: z.y, e: z.e };
    }),
  }));
}

// A sequential multi-stop scale makes the lambda ordering immediately visible:
// small lambda is pale pink, then red, purple, and finally black.
const colorStops = [
  [0.00, '#f7c4d8'],
  [0.36, '#dc2626'],
  [0.68, '#7e22ce'],
  [1.00, '#111111'],
];
function hexColor(hex) {
  return [parseInt(hex.slice(1, 3), 16), parseInt(hex.slice(3, 5), 16), parseInt(hex.slice(5, 7), 16)];
}
function mixHex(a, b, t) {
  const ca = hexColor(a), cb = hexColor(b);
  const c = ca.map((x, i) => Math.round(x + (cb[i] - x) * t));
  return '#' + c.map(x => x.toString(16).padStart(2, '0')).join('');
}
function colorForLambda(i) {
  const t = lambdas.length === 1 ? 1 : i / (lambdas.length - 1);
  for (let k = 1; k < colorStops.length; k++) {
    if (t <= colorStops[k][0]) {
      const [t0, c0] = colorStops[k - 1], [t1, c1] = colorStops[k];
      return mixHex(c0, c1, (t - t0) / (t1 - t0));
    }
  }
  return colorStops[colorStops.length - 1][1];
}
const colors = lambdas.map((_, i) => colorForLambda(i));

function fmt(v) {
  if (Math.abs(v) >= 1000 || (Math.abs(v) > 0 && Math.abs(v) < 1e-4)) return v.toExponential(2);
  if (Math.abs(v) < 10) return Number(v.toFixed(4)).toString();
  return Number(v.toFixed(2)).toString();
}
function fmtLambda(v) {
  if (Math.abs(v - Math.round(v)) < 1e-10) return String(Math.round(v));
  return v.toFixed(2).replace(/0+$/, '').replace(/\.$/, '');
}
function niceTicks(lo, hi, n = 5) {
  const span = Math.max(hi - lo, 1e-12);
  const raw = span / n;
  const mag = 10 ** Math.floor(Math.log10(raw));
  const ratio = raw / mag;
  const step = (ratio < 1.5 ? 1 : ratio < 3 ? 2 : ratio < 7 ? 5 : 10) * mag;
  const start = Math.ceil(lo / step - 1e-9) * step;
  const end = Math.floor(hi / step + 1e-9) * step;
  const out = [];
  for (let x = start; x <= end + step * 0.1; x += step) out.push(Number(x.toPrecision(12)));
  return out.length ? out : [lo, hi];
}

function drawFigure(filename, title, ylabel, data, note) {
  const W = 1900, H = 1120;
  const canvas = createCanvas(W, H);
  const ctx = canvas.getContext('2d');
  ctx.fillStyle = '#fff'; ctx.fillRect(0, 0, W, H);
  const left = 122, top = 112, right = 330, bottom = 112;
  const x0 = left, x1 = W - right, y0 = top, y1 = H - bottom;
  const xmin = Math.min(...taus), xmax = Math.max(...taus);
  let ymin = Infinity, ymax = -Infinity;
  for (const curve of data) for (const p of curve.points) {
    ymin = Math.min(ymin, p.y - p.e); ymax = Math.max(ymax, p.y + p.e);
  }
  let span = ymax - ymin;
  if (!(span > 0)) span = Math.max(Math.abs(ymax), 1) * 1e-6;
  const pad = Math.max(span * 0.10, Math.abs(ymax) * 1e-6, 1e-7);
  ymin -= pad; ymax += pad;
  const sx = x => x0 + (x - xmin) / (xmax - xmin) * (x1 - x0);
  const sy = y => y1 - (y - ymin) / (ymax - ymin) * (y1 - y0);

  ctx.fillStyle = '#111827'; ctx.font = 'bold 38px Arial'; ctx.textAlign = 'center'; ctx.textBaseline = 'top';
  ctx.fillText(title, W / 2, 28);
  ctx.fillStyle = '#374151'; ctx.font = '22px Arial';
  const yt = niceTicks(ymin, ymax, 6);
  ctx.strokeStyle = '#dbe2ea'; ctx.lineWidth = 1;
  for (const t of yt) {
    const yy = sy(t); if (yy < y0 - 2 || yy > y1 + 2) continue;
    ctx.beginPath(); ctx.moveTo(x0, yy); ctx.lineTo(x1, yy); ctx.stroke();
    ctx.textAlign = 'right'; ctx.textBaseline = 'middle'; ctx.fillText(fmt(t), x0 - 14, yy);
  }
  for (const t of taus) {
    const xx = sx(t);
    ctx.beginPath(); ctx.moveTo(xx, y0); ctx.lineTo(xx, y1); ctx.stroke();
    ctx.textAlign = 'center'; ctx.textBaseline = 'top'; ctx.fillText(fmt(t), xx, y1 + 14);
  }
  ctx.strokeStyle = '#374151'; ctx.lineWidth = 2; ctx.strokeRect(x0, y0, x1 - x0, y1 - y0);
  ctx.save(); ctx.beginPath(); ctx.rect(x0, y0, x1 - x0, y1 - y0); ctx.clip();
  data.forEach((curve, i) => {
    const color = colors[i % colors.length];
    ctx.strokeStyle = color; ctx.lineWidth = 3; ctx.beginPath();
    curve.points.forEach((p, j) => { const xx = sx(p.x), yy = sy(p.y); j ? ctx.lineTo(xx, yy) : ctx.moveTo(xx, yy); });
    ctx.stroke();
    ctx.strokeStyle = color; ctx.lineWidth = 1.5;
    for (const p of curve.points) {
      const xx = sx(p.x), yy = sy(p.y), ya = sy(p.y - p.e), yb = sy(p.y + p.e);
      ctx.beginPath(); ctx.moveTo(xx, ya); ctx.lineTo(xx, yb); ctx.moveTo(xx - 5, ya); ctx.lineTo(xx + 5, ya); ctx.moveTo(xx - 5, yb); ctx.lineTo(xx + 5, yb); ctx.stroke();
      ctx.fillStyle = color; ctx.beginPath(); ctx.arc(xx, yy, 5, 0, Math.PI * 2); ctx.fill();
    }
  });
  ctx.restore();
  ctx.fillStyle = '#111827'; ctx.font = '26px Arial'; ctx.textAlign = 'center'; ctx.textBaseline = 'top'; ctx.fillText('imaginary time τ', (x0 + x1) / 2, y1 + 58);
  ctx.save(); ctx.translate(30, (y0 + y1) / 2); ctx.rotate(-Math.PI / 2); ctx.textAlign = 'center'; ctx.textBaseline = 'top'; ctx.fillText(ylabel, 0, 0); ctx.restore();
  // Compact legend for all lambda curves.
  const lx = x1 + 34, ly = y0 + 12;
  ctx.font = '21px Arial'; ctx.textAlign = 'left'; ctx.textBaseline = 'middle';
  ctx.fillStyle = '#111827'; ctx.font = 'bold 23px Arial'; ctx.fillText('λ', lx, ly);
  ctx.font = '21px Arial';
  data.forEach((curve, i) => {
    const yy = ly + 34 + i * 31;
    ctx.strokeStyle = colors[i % colors.length]; ctx.lineWidth = 4; ctx.beginPath(); ctx.moveTo(lx, yy); ctx.lineTo(lx + 25, yy); ctx.stroke();
    ctx.fillStyle = '#111827'; ctx.fillText(fmtLambda(curve.lambda), lx + 36, yy);
  });
  ctx.fillStyle = '#4b5563'; ctx.font = '19px Arial'; ctx.textAlign = 'center'; ctx.textBaseline = 'bottom'; ctx.fillText(note, W / 2, H - 22);
  fs.writeFileSync(path.join(outDir, filename), canvas.toBuffer('image/png'));
}

function drawTriptych(filename, title, panels, note) {
  const W = 2600, H = 1220;
  const canvas = createCanvas(W, H);
  const ctx = canvas.getContext('2d');
  ctx.fillStyle = '#fff'; ctx.fillRect(0, 0, W, H);
  const left = 104, top = 116, right = 320, bottom = 128, gap = 34;
  const panelW = (W - left - right - 2 * gap) / 3;
  const panelH = H - top - bottom;
  ctx.fillStyle = '#111827'; ctx.font = 'bold 40px Arial'; ctx.textAlign = 'center'; ctx.textBaseline = 'top';
  ctx.fillText(title, W / 2, 28);

  function drawSmallPanel(panel, index) {
    const box = { x: left + index * (panelW + gap), y: top, w: panelW, h: panelH };
    const m = { l: 92, r: 22, t: 70, b: 78 };
    const px0 = box.x + m.l, px1 = box.x + box.w - m.r;
    const py0 = box.y + m.t, py1 = box.y + box.h - m.b;
    let ymin = Infinity, ymax = -Infinity;
    for (const curve of panel.data) for (const p of curve.points) {
      ymin = Math.min(ymin, p.y - p.e); ymax = Math.max(ymax, p.y + p.e);
    }
    let span = ymax - ymin;
    if (!(span > 0)) span = Math.max(Math.abs(ymax), 1) * 1e-6;
    const pad = Math.max(span * 0.10, Math.abs(ymax) * 1e-6, 1e-7);
    ymin -= pad; ymax += pad;
    const sx = x => px0 + (x - taus[0]) / (taus[taus.length - 1] - taus[0]) * (px1 - px0);
    const sy = y => py1 - (y - ymin) / (ymax - ymin) * (py1 - py0);
    ctx.fillStyle = '#111827'; ctx.font = 'bold 27px Arial'; ctx.textAlign = 'center'; ctx.textBaseline = 'top';
    ctx.fillText(panel.title, (px0 + px1) / 2, box.y + 22);
    ctx.fillStyle = '#374151'; ctx.font = '20px Arial';
    ctx.strokeStyle = '#dbe2ea'; ctx.lineWidth = 1;
    for (const t of niceTicks(ymin, ymax, 5)) {
      const yy = sy(t); if (yy < py0 - 2 || yy > py1 + 2) continue;
      ctx.beginPath(); ctx.moveTo(px0, yy); ctx.lineTo(px1, yy); ctx.stroke();
      ctx.textAlign = 'right'; ctx.textBaseline = 'middle'; ctx.fillText(fmt(t), px0 - 12, yy);
    }
    for (const t of taus) {
      const xx = sx(t); ctx.beginPath(); ctx.moveTo(xx, py0); ctx.lineTo(xx, py1); ctx.stroke();
      ctx.textAlign = 'center'; ctx.textBaseline = 'top'; ctx.fillText(fmt(t), xx, py1 + 12);
    }
    ctx.strokeStyle = '#374151'; ctx.lineWidth = 2; ctx.strokeRect(px0, py0, px1 - px0, py1 - py0);
    ctx.save(); ctx.beginPath(); ctx.rect(px0, py0, px1 - px0, py1 - py0); ctx.clip();
    panel.data.forEach((curve, i) => {
      const color = colors[i];
      ctx.strokeStyle = color; ctx.lineWidth = 3; ctx.beginPath();
      curve.points.forEach((p, j) => { const xx = sx(p.x), yy = sy(p.y); j ? ctx.lineTo(xx, yy) : ctx.moveTo(xx, yy); });
      ctx.stroke();
      ctx.strokeStyle = color; ctx.lineWidth = 1.5;
      for (const p of curve.points) {
        const xx = sx(p.x), yy = sy(p.y), ya = sy(p.y - p.e), yb = sy(p.y + p.e);
        ctx.beginPath(); ctx.moveTo(xx, ya); ctx.lineTo(xx, yb); ctx.moveTo(xx - 5, ya); ctx.lineTo(xx + 5, ya); ctx.moveTo(xx - 5, yb); ctx.lineTo(xx + 5, yb); ctx.stroke();
        ctx.fillStyle = color; ctx.beginPath(); ctx.arc(xx, yy, 5, 0, Math.PI * 2); ctx.fill();
      }
    });
    ctx.restore();
    ctx.fillStyle = '#111827'; ctx.font = '23px Arial'; ctx.textAlign = 'center'; ctx.textBaseline = 'top';
    ctx.fillText('τ', (px0 + px1) / 2, py1 + 49);
    ctx.save(); ctx.translate(box.x + 23, (py0 + py1) / 2); ctx.rotate(-Math.PI / 2); ctx.textAlign = 'center'; ctx.textBaseline = 'top'; ctx.fillText(panel.ylabel, 0, 0); ctx.restore();
  }
  panels.forEach(drawSmallPanel);

  const lx = W - right + 26, ly = top + 72;
  ctx.fillStyle = '#111827'; ctx.font = 'bold 25px Arial'; ctx.textAlign = 'left'; ctx.textBaseline = 'middle'; ctx.fillText('λ', lx, ly - 34);
  ctx.font = '21px Arial';
  panels[0].data.forEach((curve, i) => {
    const yy = ly + i * 31;
    ctx.strokeStyle = colors[i]; ctx.lineWidth = 5; ctx.beginPath(); ctx.moveTo(lx, yy); ctx.lineTo(lx + 28, yy); ctx.stroke();
    ctx.fillStyle = '#111827'; ctx.fillText(fmtLambda(curve.lambda), lx + 39, yy);
  });
  ctx.fillStyle = '#4b5563'; ctx.font = '19px Arial'; ctx.textAlign = 'center'; ctx.textBaseline = 'bottom'; ctx.fillText(note, W / 2, H - 20);
  fs.writeFileSync(path.join(outDir, filename), canvas.toBuffer('image/png'));
}

const note = 'Each point combines 3 independent seeds; whiskers are combined 1σ chain errors. L=4, 150 warmup + 500 production sweeps.';
drawFigure('average_sign_vs_tau.png', 'L=4 average sign versus imaginary time', 'average sign  Re⟨phase⟩', curves('phase'), note);
drawFigure('mz2_q1_vs_tau.png', 'L=4 VBS order parameter m_z²(q1)', 'm_z²(q1)', curves('derived', 'mz2_q1'), note);
drawFigure('mz2_q3_vs_tau.png', 'L=4 VBS order parameter m_z²(q3)', 'm_z²(q3)', curves('derived', 'mz2_q3'), note);
drawFigure('uniform_Xz_vs_tau.png', 'L=4 uniform displacement X_z', 'X_z', curves('observable', 'X_z'), note);
drawFigure('fermion_F_AB00_vs_tau.png', 'L=4 fermion correlation Im F_AB(0,0)', 'Im F_AB(0,0)', curves('fermion-imag', 'F_bare_dx0_dy0_AB'), note);
drawFigure('fermion_F_BA01_vs_tau.png', 'L=4 fermion correlation Im F_BA(0,1)', 'Im F_BA(0,1)', curves('fermion-imag', 'F_bare_dx0_dy1_BA'), note);
drawFigure('fermion_F_AB11_vs_tau.png', 'L=4 fermion correlation Im F_AB(1,1)', 'Im F_AB(1,1)', curves('fermion-imag', 'F_bare_dx1_dy1_AB'), note);

drawTriptych('order_parameters_vs_tau_triptych.png', 'L=4 order parameters versus imaginary time', [
  { title: 'VBS m_z²(q1)', ylabel: 'm_z²(q1)', data: curves('derived', 'mz2_q1') },
  { title: 'VBS m_z²(q3)', ylabel: 'm_z²(q3)', data: curves('derived', 'mz2_q3') },
  { title: 'Uniform displacement', ylabel: 'X_z', data: curves('observable', 'X_z') },
], note);
drawTriptych('fermion_correlations_vs_tau_triptych.png', 'L=4 gauge-invariant fermion correlations versus imaginary time', [
  { title: 'Im F_string_AB(0,0)', ylabel: 'Im F_string_AB(0,0)', data: curves('fermion-imag', 'F_string_dx0_dy0_AB') },
  { title: 'Im F_string_BA(0,1)', ylabel: 'Im F_string_BA(0,1)', data: curves('fermion-imag', 'F_string_dx0_dy1_BA') },
  { title: 'Im F_string_AB(1,1)', ylabel: 'Im F_string_AB(1,1)', data: curves('fermion-imag', 'F_string_dx1_dy1_AB') },
], note + ' F_string includes the directed-u path factor and is locally gauge invariant.');

for (const f of ['average_sign_vs_tau.png', 'mz2_q1_vs_tau.png', 'mz2_q3_vs_tau.png', 'uniform_Xz_vs_tau.png', 'fermion_F_AB00_vs_tau.png', 'fermion_F_BA01_vs_tau.png', 'fermion_F_AB11_vs_tau.png', 'order_parameters_vs_tau_triptych.png', 'fermion_correlations_vs_tau_triptych.png']) {
  console.log(path.join(outDir, f));
}
