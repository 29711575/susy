#!/usr/bin/env node
const fs = require('fs');
const path = require('path');
const { createCanvas } = require('C:/Users/ASUS/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/@napi-rs/canvas');

const root = path.resolve(__dirname, '..');
const runDir = path.join(root, 'runs/sign-size-lambda-L4-tau-3seed-20260918/results');
const outDir = path.join(root, 'results/sign-size-lambda-L4-tau-3seed-20260918/plots');
fs.mkdirSync(outDir, { recursive: true });
const NC = 16;

function parseCsv(text) {
  const lines = text.trim().split(/\r?\n/);
  const headers = lines[0].split(',');
  return lines.slice(1).filter(Boolean).map(line => {
    const values = line.split(',');
    const row = {};
    headers.forEach((h, i) => { row[h] = values[i]; });
    return row;
  });
}

function addVec(a, b) {
  return a.map((z, i) => [z[0] + b[i][0], z[1] + b[i][1]]);
}
function scaleVec(a, s) {
  return a.map(z => [z[0] * s, z[1] * s]);
}
function subVec(a, b) {
  return a.map((z, i) => [z[0] - b[i][0], z[1] - b[i][1]]);
}
function zeroVec() {
  return Array.from({ length: NC + 1 }, () => [0, 0]);
}
function sumBlocks(blocks) {
  return blocks.reduce((sum, block) => addVec(sum, block), zeroVec());
}

function readChain(file) {
  const rows = parseCsv(fs.readFileSync(file, 'utf8'));
  const blocks = rows.map(row => {
    const vec = [];
    vec.push([Number(row.sum_phase_real), Number(row.sum_phase_imag)]);
    for (let q = 0; q < NC; q++) {
      vec.push([
        Number(row[`sum_phase_S_q${q}_zz_real`]),
        Number(row[`sum_phase_S_q${q}_zz_imag`]),
      ]);
    }
    return vec;
  });
  if (blocks.length < 2) throw new Error(`${file}: fewer than two blocks`);
  return { file, blocks, total: sumBlocks(blocks) };
}

function ratioReal(num, den) {
  const d = den[0] * den[0] + den[1] * den[1];
  if (!(d > 0)) throw new Error('unresolved phase denominator');
  return (num[0] * den[0] + num[1] * den[1]) / d;
}

function allQMax(sum) {
  const values = [];
  for (let q = 0; q < NC; q++) values.push(ratioReal(sum[q + 1], sum[0]) / NC);
  let q = 0;
  for (let i = 1; i < values.length; i++) if (values[i] > values[q]) q = i;
  return { value: values[q], q, values };
}

function estimate(chains) {
  const total = chains.reduce((sum, c) => addVec(sum, c.total), zeroVec());
  const point = allQMax(total);
  let variance = 0;
  for (const chain of chains) {
    const b = chain.blocks.length;
    const replicas = [];
    for (const block of chain.blocks) {
      const deleted = addVec(total, scaleVec(subVec(chain.total, scaleVec(block, b)), 1 / (b - 1)));
      replicas.push(allQMax(deleted).value);
    }
    const mean = replicas.reduce((s, x) => s + x, 0) / replicas.length;
    variance += (b - 1) / b * replicas.reduce((s, x) => s + (x - mean) ** 2, 0);
  }
  return { value: point.value, stderr: Math.sqrt(variance), q: point.q, q0: point.values[0], values: point.values };
}

const groups = new Map();
for (const point of fs.readdirSync(runDir).filter(x => /^point_\d+$/.test(x))) {
  const dir = path.join(runDir, point);
  const jsonFile = fs.readdirSync(dir).find(x => x.endsWith('_results.json'));
  const meta = JSON.parse(fs.readFileSync(path.join(dir, jsonFile), 'utf8'));
  const lambda = Number(meta.lambda), tau = Number(meta.tau);
  const key = `${lambda.toPrecision(12)}|${tau.toPrecision(12)}`;
  const csvFiles = fs.readdirSync(dir).filter(x => /_blocks\.csv$/.test(x)).sort();
  if (csvFiles.length !== 3) throw new Error(`${point}: expected 3 block files, found ${csvFiles.length}`);
  if (!groups.has(key)) groups.set(key, { lambda, tau, chains: [] });
  groups.get(key).chains.push(...csvFiles.map(f => readChain(path.join(dir, f))));
}
for (const group of groups.values()) if (group.chains.length !== 3) throw new Error(`lambda=${group.lambda}, tau=${group.tau}: expected 3 chains`);

const lambdas = [...new Set([...groups.values()].map(g => g.lambda))].sort((a, b) => a - b);
const taus = [...new Set([...groups.values()].map(g => g.tau))].sort((a, b) => a - b);
const results = new Map();
for (const group of groups.values()) results.set(`${group.lambda.toPrecision(12)}|${group.tau.toPrecision(12)}`, estimate(group.chains));

const colorStops = [[0.00, '#f7c4d8'], [0.36, '#dc2626'], [0.68, '#7e22ce'], [1.00, '#111111']];
function hexColor(hex) { return [parseInt(hex.slice(1, 3), 16), parseInt(hex.slice(3, 5), 16), parseInt(hex.slice(5, 7), 16)]; }
function mixHex(a, b, t) {
  const ca = hexColor(a), cb = hexColor(b);
  const c = ca.map((x, i) => Math.round(x + (cb[i] - x) * t));
  return '#' + c.map(x => x.toString(16).padStart(2, '0')).join('');
}
function colorFor(i) {
  const t = lambdas.length === 1 ? 1 : i / (lambdas.length - 1);
  for (let k = 1; k < colorStops.length; k++) if (t <= colorStops[k][0]) {
    const [t0, c0] = colorStops[k - 1], [t1, c1] = colorStops[k];
    return mixHex(c0, c1, (t - t0) / (t1 - t0));
  }
  return colorStops[colorStops.length - 1][1];
}
const colors = lambdas.map((_, i) => colorFor(i));

function fmt(v) {
  if (Math.abs(v) >= 1000 || (Math.abs(v) > 0 && Math.abs(v) < 1e-4)) return v.toExponential(2);
  if (Math.abs(v) < 10) return Number(v.toFixed(4)).toString();
  return Number(v.toFixed(2)).toString();
}
function fmtLambda(v) { return Math.abs(v - Math.round(v)) < 1e-10 ? String(Math.round(v)) : v.toFixed(2).replace(/0+$/, '').replace(/\.$/, ''); }
function niceTicks(lo, hi, n = 6) {
  const span = Math.max(hi - lo, 1e-12), raw = span / n;
  const mag = 10 ** Math.floor(Math.log10(raw)), ratio = raw / mag;
  const step = (ratio < 1.5 ? 1 : ratio < 3 ? 2 : ratio < 7 ? 5 : 10) * mag;
  const out = [];
  for (let x = Math.ceil(lo / step - 1e-9) * step; x <= Math.floor(hi / step + 1e-9) * step + step * .1; x += step) out.push(Number(x.toPrecision(12)));
  return out.length ? out : [lo, hi];
}

function curveData() {
  return lambdas.map(lambda => ({
    lambda,
    points: taus.map(tau => ({ tau, ...results.get(`${lambda.toPrecision(12)}|${tau.toPrecision(12)}`) })),
  }));
}

function drawLegend(ctx, x, y) {
  ctx.fillStyle = '#111827'; ctx.font = 'bold 25px Arial'; ctx.textAlign = 'left'; ctx.textBaseline = 'middle'; ctx.fillText('λ', x, y - 34);
  ctx.font = '21px Arial';
  lambdas.forEach((lambda, i) => {
    const yy = y + i * 31;
    ctx.strokeStyle = colors[i]; ctx.lineWidth = 5; ctx.beginPath(); ctx.moveTo(x, yy); ctx.lineTo(x + 28, yy); ctx.stroke();
    ctx.fillStyle = '#111827'; ctx.fillText(fmtLambda(lambda), x + 39, yy);
  });
}

function drawValueFigure() {
  const data = curveData(), W = 2100, H = 1120, left = 124, right = 310, top = 112, bottom = 116;
  const canvas = createCanvas(W, H), ctx = canvas.getContext('2d');
  ctx.fillStyle = '#fff'; ctx.fillRect(0, 0, W, H);
  const x0 = left, x1 = W - right, y0 = top, y1 = H - bottom;
  let ymin = Infinity, ymax = -Infinity;
  for (const c of data) for (const p of c.points) { ymin = Math.min(ymin, p.value - p.stderr); ymax = Math.max(ymax, p.value + p.stderr); }
  const span = Math.max(ymax - ymin, 1e-8), pad = Math.max(span * .1, 1e-7); ymin -= pad; ymax += pad;
  const sx = t => x0 + (t - taus[0]) / (taus[taus.length - 1] - taus[0]) * (x1 - x0);
  const sy = v => y1 - (v - ymin) / (ymax - ymin) * (y1 - y0);
  ctx.fillStyle = '#111827'; ctx.font = 'bold 38px Arial'; ctx.textAlign = 'center'; ctx.textBaseline = 'top'; ctx.fillText('L=4 all-q maximum of m_z²(q) versus imaginary time', W / 2, 26);
  ctx.fillStyle = '#374151'; ctx.font = '22px Arial'; ctx.strokeStyle = '#dbe2ea'; ctx.lineWidth = 1;
  for (const t of niceTicks(ymin, ymax)) { const yy = sy(t); if (yy < y0 - 2 || yy > y1 + 2) continue; ctx.beginPath(); ctx.moveTo(x0, yy); ctx.lineTo(x1, yy); ctx.stroke(); ctx.textAlign = 'right'; ctx.textBaseline = 'middle'; ctx.fillText(fmt(t), x0 - 14, yy); }
  for (const t of taus) { const xx = sx(t); ctx.beginPath(); ctx.moveTo(xx, y0); ctx.lineTo(xx, y1); ctx.stroke(); ctx.textAlign = 'center'; ctx.textBaseline = 'top'; ctx.fillText(fmt(t), xx, y1 + 14); }
  ctx.strokeStyle = '#374151'; ctx.lineWidth = 2; ctx.strokeRect(x0, y0, x1 - x0, y1 - y0);
  ctx.save(); ctx.beginPath(); ctx.rect(x0, y0, x1 - x0, y1 - y0); ctx.clip();
  data.forEach((curve, i) => {
    const color = colors[i]; ctx.strokeStyle = color; ctx.lineWidth = 3; ctx.beginPath();
    curve.points.forEach((p, j) => { const xx = sx(p.tau), yy = sy(p.value); j ? ctx.lineTo(xx, yy) : ctx.moveTo(xx, yy); }); ctx.stroke();
    ctx.strokeStyle = color; ctx.lineWidth = 1.5;
    for (const p of curve.points) { const xx = sx(p.tau), yy = sy(p.value), ya = sy(p.value - p.stderr), yb = sy(p.value + p.stderr); ctx.beginPath(); ctx.moveTo(xx, ya); ctx.lineTo(xx, yb); ctx.moveTo(xx - 5, ya); ctx.lineTo(xx + 5, ya); ctx.moveTo(xx - 5, yb); ctx.lineTo(xx + 5, yb); ctx.stroke(); ctx.fillStyle = color; ctx.beginPath(); ctx.arc(xx, yy, 5, 0, Math.PI * 2); ctx.fill(); }
  });
  ctx.restore();
  ctx.fillStyle = '#111827'; ctx.font = '26px Arial'; ctx.textAlign = 'center'; ctx.textBaseline = 'top'; ctx.fillText('imaginary time τ', (x0 + x1) / 2, y1 + 56);
  ctx.save(); ctx.translate(30, (y0 + y1) / 2); ctx.rotate(-Math.PI / 2); ctx.textAlign = 'center'; ctx.textBaseline = 'top'; ctx.fillText('all-q max  m_z²(q)=max_q[S_zz(q)/N_c]', 0, 0); ctx.restore();
  drawLegend(ctx, x1 + 35, y0 + 72);
  ctx.fillStyle = '#4b5563'; ctx.font = '19px Arial'; ctx.textAlign = 'center'; ctx.textBaseline = 'bottom'; ctx.fillText('3 seeds; error bars are chain-preserving block-jackknife errors; q=0 is included in the maximum.', W / 2, H - 20);
  fs.writeFileSync(path.join(outDir, 'mz2_all_q_max_vs_tau.png'), canvas.toBuffer('image/png'));
}

function drawPeakFigure() {
  const data = curveData(), W = 2100, H = 980, left = 124, right = 310, top = 112, bottom = 110;
  const canvas = createCanvas(W, H), ctx = canvas.getContext('2d'); ctx.fillStyle = '#fff'; ctx.fillRect(0, 0, W, H);
  const x0 = left, x1 = W - right, y0 = top, y1 = H - bottom;
  const sx = t => x0 + (t - taus[0]) / (taus[taus.length - 1] - taus[0]) * (x1 - x0);
  const sy = q => y1 - q / 15 * (y1 - y0);
  ctx.fillStyle = '#111827'; ctx.font = 'bold 38px Arial'; ctx.textAlign = 'center'; ctx.textBaseline = 'top'; ctx.fillText('Momentum selected by the all-q maximum', W / 2, 26);
  ctx.fillStyle = '#374151'; ctx.font = '22px Arial'; ctx.strokeStyle = '#dbe2ea'; ctx.lineWidth = 1;
  for (let q = 0; q <= 15; q += 3) { const yy = sy(q); ctx.beginPath(); ctx.moveTo(x0, yy); ctx.lineTo(x1, yy); ctx.stroke(); ctx.textAlign = 'right'; ctx.textBaseline = 'middle'; ctx.fillText(String(q), x0 - 14, yy); }
  for (const t of taus) { const xx = sx(t); ctx.beginPath(); ctx.moveTo(xx, y0); ctx.lineTo(xx, y1); ctx.stroke(); ctx.textAlign = 'center'; ctx.textBaseline = 'top'; ctx.fillText(fmt(t), xx, y1 + 14); }
  ctx.strokeStyle = '#374151'; ctx.lineWidth = 2; ctx.strokeRect(x0, y0, x1 - x0, y1 - y0);
  ctx.strokeStyle = '#9ca3af'; ctx.lineWidth = 2; ctx.setLineDash([8, 7]); ctx.beginPath(); ctx.moveTo(x0, sy(0)); ctx.lineTo(x1, sy(0)); ctx.stroke(); ctx.setLineDash([]);
  ctx.save(); ctx.beginPath(); ctx.rect(x0, y0, x1 - x0, y1 - y0); ctx.clip();
  data.forEach((curve, i) => { const color = colors[i]; ctx.strokeStyle = color; ctx.lineWidth = 3; ctx.beginPath(); curve.points.forEach((p, j) => { const xx = sx(p.tau), yy = sy(p.q); j ? ctx.lineTo(xx, yy) : ctx.moveTo(xx, yy); }); ctx.stroke(); for (const p of curve.points) { ctx.fillStyle = color; ctx.beginPath(); ctx.arc(sx(p.tau), sy(p.q), 5, 0, Math.PI * 2); ctx.fill(); } });
  ctx.restore();
  ctx.fillStyle = '#111827'; ctx.font = '26px Arial'; ctx.textAlign = 'center'; ctx.textBaseline = 'top'; ctx.fillText('imaginary time τ', (x0 + x1) / 2, y1 + 56);
  ctx.save(); ctx.translate(32, (y0 + y1) / 2); ctx.rotate(-Math.PI / 2); ctx.textAlign = 'center'; ctx.textBaseline = 'top'; ctx.fillText('selected flattened momentum index q*', 0, 0); ctx.restore();
  drawLegend(ctx, x1 + 35, y0 + 72);
  ctx.fillStyle = '#4b5563'; ctx.font = '19px Arial'; ctx.textAlign = 'center'; ctx.textBaseline = 'bottom'; ctx.fillText('q*=0 confirms a q=0 all-q maximum; nonzero q* points are shown explicitly.', W / 2, H - 18);
  fs.writeFileSync(path.join(outDir, 'mz2_all_q_peak_index_vs_tau.png'), canvas.toBuffer('image/png'));
}

const csv = ['lambda,tau,all_q_max_mz2,stderr,q_star,mz2_q0,q0_minus_max'];
for (const lambda of lambdas) for (const tau of taus) {
  const r = results.get(`${lambda.toPrecision(12)}|${tau.toPrecision(12)}`);
  csv.push([lambda, tau, r.value, r.stderr, r.q, r.q0, r.q0 - r.value].join(','));
}
fs.writeFileSync(path.join(outDir, 'mz2_all_q_max_by_lambda_tau.csv'), csv.join('\n') + '\n');
drawValueFigure();
drawPeakFigure();
console.log(path.join(outDir, 'mz2_all_q_max_vs_tau.png'));
console.log(path.join(outDir, 'mz2_all_q_peak_index_vs_tau.png'));
console.log(path.join(outDir, 'mz2_all_q_max_by_lambda_tau.csv'));
