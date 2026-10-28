#!/usr/bin/env awk -f
# build_anchor.awk —— §177 把台账的 198 行逐行接到【验证器判决行】上
#
# 输入（-v CLED=台账 -v CCEN=普查表）：
#   matrix/forcedness_ledger.csv   9 列：paper,line,class,dev,tol,ratio,slice,invariance,criterion
#   audit/census2.csv              6 列：paper,seq,src,name,dev,tol       （src = 输出物理行号）
#
# 输出：paper,ledger_line,ckind,n_cand,status,seq,src
#   ckind  = named（有成形的判据名）| addr_fail（criterion 只写「行号指向空行」）| no_object（NA(zero-dev)）
#   chance = 同篇内与本行 (dev,tol) 相同的判决行条数（决定能不能被迫对上）
#   status = forced（恰 1 条 ⟹ 被迫）| ambiguous（≥2 条 ⟹ 不许挑）| no_key（0 条 ⟹ 表有问题）
#
# 为什么这样分：§173.9 已算过「按位置补名字」不成立（台账顺序 ≠ 实跑顺序：bsd 37/60 ·
# hl 23/35 · pvsnp 9/27 · rh 2/22），§177 又算过「按 P 记号匹配」只有 5/77 唯一
# （同一 Pn 在 part1/part2 各印一遍）。⟹ **能被迫的只有 (dev,tol) 恰唯一的那些行**，
# 其余照实写 ambiguous，不挑、不猜。
BEGIN { FS = "," ; OFS = "," }
function norm(s) { if (s ~ /^NA/) return "NA"; return sprintf("%.6e", s + 0) }
FNR == 1 { stage++ }
stage == 1 {                                  # ---- 普查表：建 (paper|dev|tol) 计数 ----
  if (FNR == 1) next
  k = $1 "|" norm($5) "|" norm($6); cnt[k]++; cand[k] = cand[k] " " $2 ":" $3
  next
}
stage == 2 {                                  # ---- 台账：逐行定位 ----
  if (FNR == 1) { print "paper,ledger_line,ckind,n_cand,status,seq,src"; next }
  if ($9 ~ /^NA\(zero-dev\)/)      ck = "no_object"
  else if ($9 ~ /points to/)       ck = "addr_fail"
  else                             ck = "named"
  k = $1 "|" norm($4) "|" norm($5); c = cnt[k] + 0
  if (c == 1) { st = "forced"; split(substr(cand[k], 2), pp, ":"); seq = pp[1]; src = pp[2] }
  else if (c >= 2) { st = "ambiguous"; seq = "NA"; src = "NA" }
  else { st = "no_key"; seq = "NA"; src = "NA" }
  print $1, $2, ck, c, st, seq, src
}
