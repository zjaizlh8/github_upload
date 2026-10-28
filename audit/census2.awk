#!/usr/bin/env awk -f
# census2.awk —— §177 判决行普查【带输出行号】，给台账造一根不会腐的锚
#
# 为什么需要它：台账第 2 列 `line` 指的是【论文 Markdown 的行号】，它随写作而腐
# （§172 实测 180 个不同地址里 68 个落空 = 37.8%）。而验证器输出的一行
# `PASS <名> dev=… tol=…` 是这台机器上【同一个二进制】的固定行位 ——
# 只要代码没改，它就不动；代码改了，行号变了也正是【该重新核账】的信号。
#
# 输出列：paper,seq,src,name,dev,tol
#   seq = 本篇第几条判决行（1 起）
#   src = 在 audit/outs/<paper>.txt 里的【物理行号】（1 起）
# 判据：抽出的行必须仍以判决词开头（PASS / [PASS] / **FAIL**），RESULT 汇总行不算。
BEGIN { FS = "," }
FNR == 1 {
  split(FILENAME, parts, "/"); p = parts[length(parts)]; sub(/\.txt$/, "", p); seq = 0
}
{
  line = $0
  if (line ~ /^[ \t]*RESULT/) next
  if (line !~ /^[ \t]*(\[PASS\]|PASS|\*\*FAIL\*\*)/) next
  seq++
  nm = line
  sub(/^[ \t]*(\[PASS\]|PASS|\*\*FAIL\*\*)[ \t]*/, "", nm)
  dev = "NA"; tol = "NA"
  if (match(nm, /dev=[-0-9.eE+]+[ \t]*\(?tol=[-0-9.eE+]+\)?/)) {
    pair = substr(nm, RSTART, RLENGTH)
    tmp = pair; sub(/^dev=/, "", tmp); sub(/[ \t].*/, "", tmp); dev = tmp
    tmp = pair; sub(/^.*tol=/, "", tmp); sub(/\).*$/, "", tmp); tol = tmp
    nm = substr(nm, 1, RSTART - 1)
  }
  gsub(/[ \t]+$/, "", nm); gsub(/^[ \t]+/, "", nm); gsub(/,/, ";", nm)
  printf "%s,%d,%d,%s,%s,%s\n", p, seq, FNR, nm, dev, tol
}
