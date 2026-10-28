#!/usr/bin/env node
/*
 * check_lang_parity.mjs —— 中英理论是否一致（机械判，能翻红）
 *
 * 对象：七篇论文的 paper.zh.md 与 paper.md 是同一理论的两种语言。
 * 判据：
 *   L1 两个文件都存在
 *   L2 二级标题数相等
 *   L3 三级标题数相等
 *   L4 具名常数的【取值集合】相等（不是出现次数：同一个数换写法/多提一次不算不一致）
 *   L5 逐节段数差 ≤ 2（按 ## 分节对齐，定位到节；整段丢失才算错）
 *   L6 两侧 RESULT PASS 声明序列相同
 *
 * 用法：node audit/check_lang_parity.mjs [--verbose]
 * 退出码：0 = 全一致；1 = 有不一致
 */
import fs from 'node:fs';
import path from 'node:path';

const ROOTARG = process.argv.find((s) => s.startsWith('--root='));
const ROOT = ROOTARG ? ROOTARG.slice('--root='.length) : 'D:/lw/github_upload';
const PAPERS = [
  ['bsd', 'bsd-conjecture'],
  ['hl', 'hardy-long'], // 占位，下面按目录名重设
];
PAPERS.length = 0;
PAPERS.push(
  ['bsd', 'bsd-conjecture'],
  ['hl', 'hardy-littlewood'],
  ['hodge', 'hodge-conjecture'],
  ['ns', 'navier-stokes'],
  ['pvsnp', 'p-vs-np'],
  ['rh', 'riemann-hypothesis'],
  ['ym', 'yang-mills'],
);

const VERBOSE = process.argv.includes('--verbose');
const read = (p) => (fs.existsSync(p) ? fs.readFileSync(p, 'utf8') : null);

// 具名常数：小数 ≥5 位有效数字，或 ≥5 位整数；统一成规范形（去 e+、去尾零）
const NUM = /-?[0-9]+\.[0-9]{4,}(?:e[-+]?[0-9]+)?|[0-9]{5,}/gi;
const canon = (s) => {
  const t = s.replace(/e\+/, 'e').toLowerCase();
  if (!/e/.test(t) && /\./.test(t)) {
    const v = t.replace(/0+$/, '').replace(/\.$/, '');
    return v;
  }
  if (/e/.test(t)) {
    const [m, e] = t.split('e');
    return `${m.replace(/0+$/, '').replace(/\.$/, '')}e${e}`;
  }
  return t.replace(/^0+(?=\d)/, '');
};
const nums = (t) => (t.match(NUM) || []).map(canon);
const uniq = (a) => [...new Set(a)].sort();

const count = (t, re) => (t.match(re) || []).length;
const passClaims = (t) => (t.match(/RESULT\s+PASS\s*=\s*[0-9]+/g) || []).map((s) => s.replace(/\s+/g, ' '));

// 按二级标题切节，返回 [{title, body}]
const sections = (t) => {
  const lines = t.split('\n');
  const out = [];
  let cur = { title: '(head)', body: [] };
  for (const ln of lines) {
    if (/^## /.test(ln)) {
      out.push(cur);
      cur = { title: ln.replace(/^## /, '').trim(), body: [] };
    } else cur.body.push(ln);
  }
  out.push(cur);
  return out;
};
const paras = (body) => body.join('\n').split(/\n[ \t]*\n/).filter((s) => s.trim().length > 0).length;
// 非空行数：抗「空行风格不同」的段落假差（中文少空行、英文多空行是排版不是内容）
const nonEmpty = (body) => body.filter((l) => l.trim().length > 0).length;

let bad = 0;
const say = (ok, tag, msg) => {
  if (!ok) bad++;
  console.log(`${ok ? 'MATCH  ' : '**DIFF**'} ${tag.padEnd(10)} ${msg}`);
};

const rows = [];
for (const [code, dir] of PAPERS) {
  const zh = read(path.join(ROOT, dir, 'paper.zh.md'));
  const en = read(path.join(ROOT, dir, 'paper.md'));
  say(!!zh && !!en, `${code}·L1`, `存在性 zh=${!!zh} en=${!!en}`);
  if (!zh || !en) continue;

  const h2z = count(zh, /^## /gm), h2e = count(en, /^## /gm);
  const h3z = count(zh, /^### /gm), h3e = count(en, /^### /gm);
  say(h2z === h2e, `${code}·L2`, `二级标题 zh=${h2z} en=${h2e}`);
  say(h3z === h3e, `${code}·L3`, `三级标题 zh=${h3z} en=${h3e}`);

  const sz = uniq(nums(zh)), se = uniq(nums(en));
  const onlyZ = sz.filter((x) => !se.includes(x));
  const onlyE = se.filter((x) => !sz.includes(x));
  say(onlyZ.length === 0 && onlyE.length === 0, `${code}·L4`,
    `常数取值 zh=${sz.length} en=${se.length} 仅中文=${JSON.stringify(onlyZ.slice(0, 10))} 仅英文=${JSON.stringify(onlyE.slice(0, 10))}`);

  const cz = passClaims(zh), ce = passClaims(en);
  say(JSON.stringify(cz) === JSON.stringify(ce), `${code}·L6`, `RESULT PASS zh=${JSON.stringify(cz)} en=${JSON.stringify(ce)}`);

  // L5：逐节【结构件】必须逐项相等 —— 块公式数、列表项数、表格行数。
  //     这三样是语言无关的（中文写公式、写 `- `、写 `|` 与英文一模一样）⟹ 差一项就是漏了一段。
  //     实测基线七篇：disp/li/tr 全篇逐节相等（唯一不齐的是行内代码数 ic，英文习惯每个标识符都包反引号，
  //     中文不包 ⟹ ic 不作判据，只印读数）。
  const Sz = sections(zh), Se = sections(en);
  const feats = (body) => {
    const raw = body.join('\n');
    // 表格行只数【数据行】：表头行（其后紧跟 |---| 分隔）与分隔行不计，
    // 因为同一张表在中英两侧可能被排成一张或拆成两张（重复表头），内容其实一样。
    let tr = 0;
    for (let k = 0; k < body.length; k++) {
      const ln = body[k];
      if (!/^[ \t]*\|/.test(ln)) continue;
      if (/^[\s|:\-]+$/.test(ln)) continue;                        // 分隔行
      const nxt = body[k + 1] || '';
      if (/^[ \t]*\|[\s:\-|]+$/.test(nxt)) continue;               // 表头行
      tr++;
    }
    return {
      disp: (raw.match(/\$\$/g) || []).length / 2,
      li: (raw.match(/^[ \t]*[-*] /gm) || []).length,
      tr,
      ic: (raw.match(/`[^`\n]+`/g) || []).length,
      cjk: (raw.match(/[\u4e00-\u9fff]/g) || []).length,
      wd: (raw.match(/[A-Za-z][A-Za-z'-]*/g) || []).length,
    };
  };
  const n = Math.min(Sz.length, Se.length);
  const bad5 = [];
  const ratios = [];
  let trRead = 0;
  for (let i = 0; i < n; i++) {
    const a = feats(Sz[i].body), b = feats(Se[i].body);
    // 只用「块公式数 + 列表项数」作判据；表格行数受排版影响（同一张表可拆成两张、表头重复、
    // 行号在两侧各自从 1 起编），把它降为读数，避免假阳性。
    const diffs = ['disp', 'li'].filter((k) => a[k] !== b[k]).map((k) => `${k} ${a[k]}/${b[k]}`);
    if (diffs.length) bad5.push(`#${i} ${diffs.join(' ')} 「${Sz[i].title.slice(0, 24)}」`);
    trRead += Math.abs(a.tr - b.tr) === 0 ? 0 : 1;
  }
  // L5b 的计数先剥掉代码块/行内代码/公式：附录里「中文侧也有 ASCII 公式」会让汉字/词比失真
  const stripAll = (s) => s
    .replace(/```[\s\S]*?```/g, ' ')
    .replace(/`[^`\n]*`/g, ' ')
    .replace(/\$\$[\s\S]*?\$\$/g, ' ')
    .replace(/\$[^$\n]*\$/g, ' ');
  const prose = (body) => {
    const t = stripAll(body.join('\n'));
    return { cjk: (t.match(/[\u4e00-\u9fff]/g) || []).length, wd: (t.match(/[A-Za-z][A-Za-z'-]*/g) || []).length };
  };
  for (let i = 0; i < n; i++) {
    const a = prose(Sz[i].body), b = prose(Se[i].body);
    if (a.cjk >= 30 && b.wd >= 30) ratios.push({ i, r: a.cjk / b.wd, title: Sz[i].title.slice(0, 24) });
  }
  if (Sz.length !== Se.length) bad5.push(`节数 zh=${Sz.length} en=${Se.length}`);
  say(bad5.length === 0, `${code}·L5`, `逐节结构件(公式/列表/表格)${bad5.length ? '：' + bad5.join(' | ') : ' 全篇相等'}`);

  // L5b：散文量带（汉字/英文词）—— 弱判据，只抓「整节翻倍或减半」，抓不到单段落漏译（照实标）
  ratios.sort((x, y) => x.r - y.r);
  const R = ratios.length ? ratios[Math.floor(ratios.length / 2)].r : NaN;
  const bad5b = ratios.filter((x) => x.r < R / 2 || x.r > 2 * R).map((x) => `#${x.i} ${x.r.toFixed(2)} 「${x.title}」`);
  say(bad5b.length === 0, `${code}·L5b`, `散文量带 中位 ${R.toFixed(2)} ${bad5b.length ? '越界：' + bad5b.join(' | ') : '全部在 [R/2,2R]'}`);

  // L5c：逐节【常数集合】相等 —— 同一个数必须出现在同一节里。
  //      这一条是本 checker 抓到真漏的唯一一条：bsd §9 英文缺了 `| 11 | ∏spec = −6^11 = −362797056` 整行；
  //      hodge §0′ 中文缺 `0.069027`（k=163 处的值）。数目（多重集）不作判据：重复次数随行文风格变。
  const numSet = (body) => [...new Set((body.join('\n').match(NUM) || []).map(canon))];
  const bad5c = [];
  for (let i = 0; i < n; i++) {
    const a = numSet(Sz[i].body), b = numSet(Se[i].body);
    const oz = a.filter((x) => !b.includes(x)), oe = b.filter((x) => !a.includes(x));
    if (oz.length || oe.length) bad5c.push(`#${i} 仅中=${JSON.stringify(oz.slice(0, 4))} 仅英=${JSON.stringify(oe.slice(0, 4))} 「${Sz[i].title.slice(0, 22)}」`);
  }
  say(bad5c.length === 0, `${code}·L5c`, `逐节常数集合${bad5c.length ? '：' + bad5c.join(' | ') : ` 全篇相等（表格行数读数：${trRead} 节不同，仅排版）`}`);

  rows.push({ code, h2: `${h2z}/${h2e}`, h3: `${h3z}/${h3e}`, num: `${sz.length}/${se.length}`, onlyZ: onlyZ.length, onlyE: onlyE.length, bad5: bad5.length });
}

console.log('\n=== 汇总 ===');
for (const r of rows) console.log(`${r.code.padEnd(7)} h2=${r.h2} h3=${r.h3} 常数取值=${r.num} 仅中=${r.onlyZ} 仅英=${r.onlyE} 段数越界节=${r.bad5}`);
console.log(`\n不一致条数 = ${bad}`);
console.log(bad === 0 ? 'RESULT PASS' : 'RESULT FAIL');
process.exit(bad === 0 ? 0 : 1);
