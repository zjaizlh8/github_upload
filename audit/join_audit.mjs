// audit/join_audit.mjs —— 画像 join 轴的上游审计：对齐了但读数不同的行
// 问题：profile_v2.csv 的 LCS 只按 dev 相等配对；配对后未复核 tol；
//       且 none 的 4 行究竟是谁。本脚本按 (paper,ledger_line) 回读台账，逐行比较。
import { readFileSync } from 'node:fs';

function rowsCSV(text) {
  const lines = text.split(/\r?\n/).filter((s) => s.length > 0);
  const hdr = lines[0].split(',');
  return lines.slice(1).map((ln) => {
    const f = ln.split(',');
    while (f.length > hdr.length) {
      const tail = f.pop();
      f[f.length - 1] += ',' + tail;
    }
    const o = {};
    hdr.forEach((h, i) => (o[h.trim()] = (f[i] ?? '').trim()));
    return o;
  });
}
const num = (s) => (s === 'NA' || s === 'NA(implicit)' ? NaN : Number(s));
const eq = (a, b) => {
  const x = num(a), y = num(b);
  if (isNaN(x) && isNaN(y)) return true;
  if (isNaN(x) || isNaN(y)) return false;
  return x === y;
};

const prof = rowsCSV(readFileSync('D:/lw/github_upload/matrix/profile_v2.csv', 'utf8'));
const led = rowsCSV(readFileSync('D:/lw/github_upload/matrix/forcedness_ledger.csv', 'utf8'));
const ledBy = new Map();
for (const r of led) ledBy.set(`${r.paper}|${r.line}`, r);

let aligned = 0, devSame = 0, devDiff = [], tolDiff = [], none = [];
for (const r of prof) {
  if (r.join === 'none') { none.push(r); continue; }
  const l = ledBy.get(`${r.paper}|${r.ledger_line}`);
  if (!l) { none.push({ ...r, why: 'no ledger row' }); continue; }
  aligned++;
  if (eq(r.dev, l.dev)) devSame++;
  else devDiff.push({ paper: r.paper, seq: r.seq, cen: r.dev, led: l.dev, name: r.name.slice(0, 28), join: r.join });
  if (!eq(r.tol, l.tol)) tolDiff.push({ paper: r.paper, seq: r.seq, cen: r.tol, led: l.tol, name: r.name.slice(0, 28) });
}
console.log(`ROWS=${prof.length} ALIGNED=${aligned} NONE=${none.length}`);
console.log(`dev 相同=${devSame} · dev 不同=${devDiff.length} · tol 不同=${tolDiff.length}`);
console.log('\n== 对齐了但 dev 不同 ==');
for (const d of devDiff) console.log(`  ${d.paper}[${d.seq}] ${d.join} 输出=${d.cen} 台账=${d.led}  「${d.name}」`);
console.log('\n== 对齐了但 tol 不同（前 15）==');
for (const d of tolDiff.slice(0, 15)) console.log(`  ${d.paper}[${d.seq}] 输出=${d.cen} 台账=${d.led}  「${d.name}」`);
console.log(`  …总 ${tolDiff.length}`);
console.log('\n== 未对齐的 4 行 ==');
for (const n of none) console.log(`  ${n.paper}[${n.seq}] dev=${n.dev} tol=${n.tol} 「${(n.name || '').slice(0, 40)}」`);
