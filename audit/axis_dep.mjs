// audit/axis_dep.mjs —— 画像 v2 的「轴依赖图」：哪些列是别的列的函数
// 判据：B 是 A 的函数 ⟺ 每个 A 值只对应一个 B 值（|B values per A| == 1）
// 输出：方向图 + 最小生成轴集 + 每条独立轴的行内取值表（供 G54 三分）
import { readFileSync } from 'node:fs';

const P = 'D:/lw/github_upload/matrix/profile_v2.csv';
const L = 'D:/lw/github_upload/matrix/forcedness_ledger.csv';

function rowsCSV(text) {
  const lines = text.split(/\r?\n/).filter((s) => s.length > 0);
  const hdr = lines[0].split(',');
  return lines.slice(1).map((ln) => {
    const f = ln.split(',');
    // criterion 含逗号：多于表头长度时把尾部拼回最后一列
    while (f.length > hdr.length) {
      const tail = f.pop();
      f[f.length - 1] += ',' + tail;
    }
    const o = {};
    hdr.forEach((h, i) => (o[h.trim()] = (f[i] ?? '').trim()));
    return o;
  });
}

// ── dev 二型判定：R = 落在 k·2^-52 网格且 k=m/2^j（m≤64, j≤8）；A = 其他非零；ZERO / IMPLICIT
const ULP = Math.pow(2, -52);
function devType(s) {
  if (s === 'NA(implicit)') return 'IMPLICIT';
  const v = Number(s);
  if (!isFinite(v)) return 'BAD';
  if (v === 0) return 'ZERO';
  const k = v / ULP;
  for (let j = 0; j <= 8; j++) {
    const m = Math.abs(k) * Math.pow(2, j);
    const r = Math.abs(m - Math.round(m));
    if (Math.abs(m) <= 64 + 1e-9 && r <= 5e-6 * Math.max(1, Math.abs(m))) return 'R';
  }
  return 'A';
}
function dec(s) {
  if (!isFinite(Number(s))) return 'NA';
  const v = Math.abs(Number(s));
  if (v === 0) return 'zero';
  return String(Math.floor(Math.log10(v)));
}

const prof = rowsCSV(readFileSync(P, 'utf8'));
const led = rowsCSV(readFileSync(L, 'utf8'));
const ledBy = new Map();
for (const r of led) ledBy.set(`${r.paper}|${r.line}`, r);

// 把 profile 与台账合并成一张特征表
const F = prof.map((r) => {
  const l = ledBy.get(`${r.paper}|${r.ledger_line}`) || {};
  const nm = r.criterion_final || '';
  const nameState = nm.startsWith('NA(addr_dead)') ? 'no_object' : nm.startsWith('paper:') ? 'addr_dead' : 'named';
  return {
    paper: r.paper,
    seq: Number(r.seq),
    class: r.ledger_class,
    inv: r.ledger_inv,
    slice: (l.slice || '').replace(/;/g, '+'),
    devT: devType(r.dev),
    devDec: dec(r.dev),
    tolDec: dec(r.tol),
    nameSt: nameState,
    addr: r.addr,
    join: r.join,
  };
});

const AX = ['paper', 'class', 'inv', 'slice', 'devT', 'devDec', 'tolDec', 'nameSt', 'addr', 'join'];
const vals = (a) => new Set(F.map((r) => r[a]));
const nv = Object.fromEntries(AX.map((a) => [a, vals(a).size]));

// 函数依赖 A → B：按 A 分组，看 B 的取值个数是否恒为 1
function isFn(a, b) {
  const g = new Map();
  for (const r of F) {
    if (!g.has(r[a])) g.set(r[a], new Set());
    g.get(r[a]).add(r[b]);
  }
  for (const s of g.values()) if (s.size !== 1) return false;
  return true;
}

const deps = [];
for (const a of AX) for (const b of AX) if (a !== b && isFn(a, b)) deps.push([a, b]);

// 最小生成集：从全集中贪心去掉能被其他轴决定者
let keep = AX.filter((a) => !deps.some(([x, y]) => y === a && x !== 'seq'));
// 若 a→b 与 b→a 同时成立，两者等价，只留一个（保留基数大的）
const eq = deps.filter(([a, b]) => deps.some(([x, y]) => x === b && y === a));
for (const [a, b] of eq) {
  if (keep.includes(a) && keep.includes(b)) {
    const drop = nv[a] >= nv[b] ? b : a;
    keep = keep.filter((k) => k !== drop);
  }
}

console.log('== 轴的取值基数 ==');
for (const a of AX) console.log(`  ${a.padEnd(8)} ${nv[a]}`);
console.log('\n== 函数依赖（A → B，A 决定 B）==');
if (!deps.length) console.log('  （无）');
for (const [a, b] of deps) console.log(`  ${a} → ${b}`);
console.log('\n== 互相决定（等价轴对）==');
for (const [a, b] of eq) if (a < b) console.log(`  ${a} ⟺ ${b}`);
console.log('\n== 最小生成集 ==');
console.log('  ' + keep.join(' · ') + `   （${AX.length} → ${keep.length}）`);
console.log('\n== 独立轴的行内分布（供 G54 三分）==');
for (const a of keep) {
  const m = new Map();
  for (const r of F) m.set(r[a], (m.get(r[a]) || 0) + 1);
  const s = [...m.entries()].sort((x, y) => y[1] - x[1]).map(([k, v]) => `${k}:${v}`).join(' · ');
  console.log(`  [${a}] ${s}`);
}
console.log(`\nROWS=${F.length} AXES=${AX.length} DEPS=${deps.length} MINIMAL=${keep.length}`);
