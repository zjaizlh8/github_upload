#!/usr/bin/env node
/*
 * tooth_lang_parity.mjs —— check_lang_parity.mjs 的判别力齿测（并发，各用独立副本）
 *
 * 每支齿对应一条判据，齿必须红，否则该判据是自比：
 *   T1 基线（不改）        ⟹ 期望 PASS
 *   T2 中文插入一条块公式  ⟹ 期望 L5 红（disp 不等）
 *   T3 英文删掉列表项      ⟹ 期望 L5 红（li 不等）
 *   T4 英文常数末位改动    ⟹ 期望 L4 红（常数取值集合不等）
 *   T5 英文三级标题降级    ⟹ 期望 L3 红（三级标题数不等）
 *   T6 英文删掉一个读数    ⟹ 期望 L5c 红（该节常数集合缺项）
 *   T7 中文散文扩到三倍    ⟹ 期望 L5b 红（散文量带越界）
 *
 * 用法：node audit/tooth_lang_parity.mjs
 */
import fs from 'node:fs';
import path from 'node:path';
import { execFile } from 'node:child_process';
import { promisify } from 'node:util';

const pexec = promisify(execFile);
const SRC = 'D:/lw/github_upload';
const BASE = 'D:/dpgo/cx/_tmp/lang_parity';
const DIRS = ['bsd-conjecture', 'hardy-littlewood', 'hodge-conjecture', 'navier-stokes', 'p-vs-np', 'riemann-hypothesis', 'yang-mills'];
const CHECK = path.join(SRC, 'audit', 'check_lang_parity.mjs');

const clone = (tag) => {
  const root = path.join(BASE, tag);
  for (const d of DIRS) {
    fs.mkdirSync(path.join(root, d), { recursive: true });
    for (const f of ['paper.md', 'paper.zh.md']) {
      const s = path.join(SRC, d, f);
      if (fs.existsSync(s)) fs.copyFileSync(s, path.join(root, d, f));
    }
  }
  return root;
};

const run = async (root) => {
  try {
    const { stdout } = await pexec(process.execPath, [CHECK, `--root=${root}`], { encoding: 'utf8', maxBuffer: 8 << 20 });
    return { code: 0, out: stdout };
  } catch (e) {
    return { code: e.code ?? 1, out: (e.stdout || '') + (e.stderr || '') };
  }
};

// 在指定文件的第 n 个二级节末尾插入/改写
const idxOfSection = (lines, nth) => {
  let seen = -1;
  for (let i = 0; i < lines.length; i++) {
    if (/^## /.test(lines[i])) {
      seen++;
      if (seen === nth) {
        let j = i + 1;
        while (j < lines.length && !/^## /.test(lines[j])) j++;
        return { start: i, end: j };
      }
    }
  }
  return null;
};

const mutate = {
  t1: () => {},
  t2: (root) => {                                            // zh 插入块公式
    const p = path.join(root, 'hodge-conjecture', 'paper.zh.md');
    const L = fs.readFileSync(p, 'utf8').split('\n');
    const s = idxOfSection(L, 3);
    L.splice(s.end, 0, '', '$$\\boxed{\\text{齿测插入的块公式}}$$', '');
    fs.writeFileSync(p, L.join('\n'));
  },
  t3: (root) => {                                            // en 删整节列表项（§4，该节有 4 条）
    const p = path.join(root, 'yang-mills', 'paper.md');
    const L = fs.readFileSync(p, 'utf8').split('\n');
    const s = idxOfSection(L, 7);
    const seg = L.slice(s.start, s.end);
    const kept = seg.filter((l) => !/^[ \t]*[-*] /.test(l));
    const removed = seg.length - kept.length;
    if (removed < 2) throw new Error(`T3 选错节：本节只有 ${removed} 条列表项，删不掉就测不出 L5 li`);
    fs.writeFileSync(p, [...L.slice(0, s.start), ...kept, ...L.slice(s.end)].join('\n'));
    return `删 ${removed} 条列表项`;
  },
  t4: (root) => {                                            // en 常数末位
    const p = path.join(root, 'yang-mills', 'paper.md');
    const t = fs.readFileSync(p, 'utf8');
    const m = t.match(/[0-9]+\.[0-9]{4,}/);
    const v = m[0];
    fs.writeFileSync(p, t.replace(v, v.slice(0, -1) + '9'));
    return `${v} → ${v.slice(0, -1)}9`;
  },
  t5: (root) => {                                            // en 三级标题降级
    const p = path.join(root, 'riemann-hypothesis', 'paper.md');
    const t = fs.readFileSync(p, 'utf8');
    const i = t.indexOf('\n### ');
    fs.writeFileSync(p, t.slice(0, i) + t.slice(i).replace('\n### ', '\n## '));
  },
  t6: (root) => {                                            // en 删掉某节的一个读数
    const p = path.join(root, 'hodge-conjecture', 'paper.md');
    const L = fs.readFileSync(p, 'utf8').split('\n');
    const k = L.findIndex((l) => l.includes('0.069027'));
    L[k] = L[k].replace('0.069027', '(removed for tooth test)');
    fs.writeFileSync(p, L.join('\n'));
  },
  t7: (root) => {                                            // zh 散文扩到三倍
    const p = path.join(root, 'navier-stokes', 'paper.zh.md');
    const L = fs.readFileSync(p, 'utf8').split('\n');
    const s = idxOfSection(L, 6);
    const filler = '齿测填充段：把本节的中文内容量推到中位比值的三倍以上，用来确认散文量带确实会越界；若它不红，该带就是自比。';
    L.splice(s.end, 0, '', filler + filler + filler, '');
    fs.writeFileSync(p, L.join('\n'));
  },
};

const TAGS = Object.keys(mutate);
const roots = {};
for (const tag of TAGS) roots[tag] = clone(tag);
const notes = {};
for (const tag of TAGS) notes[tag] = mutate[tag](roots[tag]) || '';
const results = await Promise.all(TAGS.map((tag) => run(roots[tag])));

const res = {};
TAGS.forEach((tag, i) => (res[tag] = results[i]));

let bad = 0;
const say = (ok, tag, msg) => {
  if (!ok) bad++;
  console.log(`${ok ? 'PASS' : '**FAIL**'} ${tag.padEnd(3)} ${msg}`);
};
const isRed = (tag, key) => res[tag].code === 1 && new RegExp(`${key}·L`).test(res[tag].out) && /\*\*DIFF\*\*/.test(res[tag].out);

say(res.t1.code === 0 && /RESULT PASS/.test(res.t1.out), 'T1', `基线 ⟹ PASS（exit=${res.t1.code}）`);
say(isRed('t2', 'hodge') && /hodge·L5 /.test(res.t2.out), 'T2', `zh 插入块公式 ⟹ 期望 L5 红（exit=${res.t2.code}）`);
say(isRed('t3', 'ym') && /ym·L5 /.test(res.t3.out), 'T3', `en 删列表项 ⟹ 期望 L5 红（exit=${res.t3.code}）`);
say(isRed('t4', 'ym') && /ym·L4/.test(res.t4.out), 'T4', `en 常数 ${notes.t4} ⟹ 期望 L4 红（exit=${res.t4.code}）`);
say(isRed('t5', 'rh') && /rh·L3/.test(res.t5.out), 'T5', `en 标题降级 ⟹ 期望 L3 红（exit=${res.t5.code}）`);
say(isRed('t6', 'hodge') && /hodge·L5c/.test(res.t6.out), 'T6', `en 删一个读数 ⟹ 期望 L5c 红（exit=${res.t6.code}）`);
say(isRed('t7', 'ns') && /ns·L5b/.test(res.t7.out), 'T7', `zh 散文三倍 ⟹ 期望 L5b 红（exit=${res.t7.code}）`);

for (const tag of TAGS) {
  const lines = res[tag].out.split('\n').filter((l) => l.includes('**DIFF**'));
  if (lines.length) console.log(`      ${tag}: ${lines.map((s) => s.replace(/\s+/g, ' ').trim()).join(' || ').slice(0, 260)}`);
}

fs.rmSync(BASE, { recursive: true, force: true });
console.log(`\n齿测 RESULT ${bad === 0 ? 'PASS' : 'FAIL'} FAIL=${bad}`);
process.exit(bad === 0 ? 0 : 1);
