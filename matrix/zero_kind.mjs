// matrix/zero_kind.mjs —— 把「dev=0」拆开：残差恰为零 vs 根本不是残差（违例计数/无残差可报）
// 判据：criterion 文本给出「不是残差」的自我声明 ⟹ 该行的 dev 不参与残差统计。
// 输出：三类计数 + 逐行清单（供 G54 三分与文档落盘）。
import { readFileSync } from 'node:fs';

function rowsCSV(text) {
  const lines = text.split(/\r?\n/).filter((s) => s.length);
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
const led = rowsCSV(readFileSync('D:/lw/github_upload/matrix/forcedness_ledger.csv', 'utf8'));
const isZero = (d) => d === '0.000000e+00' || d === '0.000e+00';
const isImpl = (d) => d === 'NA(implicit)';

// 「不是残差」的自我声明词汇（来自台账 criterion 文本，逐条可核）
const NOTRES = /违例|violation|计数|三类全对|无残差|not a residual|count/i;

const cat = { geo: [], count: [], impl: [], other: [] };
for (const r of led) {
  const d = r.dev, c = r.criterion || '';
  if (isImpl(d)) { cat.impl.push(r); continue; }
  if (!isZero(d)) { cat.other.push(r); continue; }
  if (NOTRES.test(c)) cat.count.push(r);
  else cat.geo.push(r);
}
const pr = (k, v) => console.log(`${k}\t${v.length}`);
pr('几何零（残差恰为 0）', cat.geo);
pr('计数零（dev 是违例数不是残差）', cat.count);
pr('隐含零（判决行无 dev token）', cat.impl);
pr('非零残差', cat.other);
console.log('合计', led.length);
console.log('\n== 计数零逐行 ==');
for (const r of cat.count) console.log(`  ${r.paper}:${r.line} ${r.class}/${r.invariance} 「${(r.criterion || '').slice(0, 70)}」`);
console.log('\n== 分篇：几何零 / 计数零 / 隐含零 ==');
const per = {};
for (const k of ['geo', 'count', 'impl']) for (const r of cat[k]) { per[r.paper] ??= { geo: 0, count: 0, impl: 0 }; per[r.paper][k]++; }
for (const p of ['bsd', 'hl', 'hodge', 'ns', 'pvsnp', 'rh', 'ym']) {
  const v = per[p] || { geo: 0, count: 0, impl: 0 };
  console.log(`  ${p.padEnd(6)} geo=${v.geo} count=${v.count} impl=${v.impl} sum=${v.geo + v.count + v.impl}`);
}
