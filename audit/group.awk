#!/usr/bin/env awk -f
# group.awk —— §177 把 198 条判决行按【块 + 标签 + 名称骨架】归成判据组
#
# 为什么必须做这一步：
#   1. `Pn` 在各篇里【不是标识符】—— 同一篇内 part1 与 part2 各自从 P1 编号，
#      实测 46 组「同号不同物」（§177.3）。所以 `P5` 单独一个记号指不到对象。
#   2. `dev` 也不是标识符 —— 零残差在所有恒等判据上相同（§158.3 已证）。
#   3. 于是能当键的只有【块 + 标签 + 骨架】：块 = 输出里的 `[part k/m]` 段，
#      骨架 = 名称去掉数字与运算符后的文字。同块同骨架 ⟹ 同一个判据在多个参数点上跑。
#
# 输出（stdout，TSV）：paper,part,label,skeleton,rows,dev_set
#   rows    = 该组用了几行（= 参数点数）
#   dev_set = 该组出现过的不同 dev 值的个数（1 ⟹ 各参数点读数相同 ⟹ 需查是否恒真）
#
# 末尾向 stderr 打一行总账，供 verify 脚本取用。
BEGIN { FS = ","; hdr = 1 }
function skel(s,   o) {
  o = s; gsub(/[0-9]/, "", o); gsub(/[=+*\/()（）·、,;:.\-—–_{}^$|]/, " ", o)
  gsub(/ +/, " ", o); sub(/^ /, "", o); sub(/ $/, "", o)
  if (o == "") o = "(no-name)"
  return o
}
FNR == 1 { split(FILENAME, q, "/"); p = q[length(q)]; sub(/\.txt$/, "", p); part = 0; next }
{
  if ($0 ~ /part [0-9]+ *\/ *[0-9]+/) { part++; next }            # 块头 ⟹ 新块
  if ($0 ~ /^[ \t]*RESULT/) next
  if ($0 !~ /^[ \t]*(\[PASS\]|PASS|\*\*FAIL\*\*)/) next
  nm = $0
  sub(/^[ \t]*(\[PASS\]|PASS|\*\*FAIL\*\*)[ \t]*/, "", nm)
  dev = "NA"
  if (match(nm, /dev=[-0-9.eE+]+/)) { dev = substr(nm, RSTART + 4, RLENGTH - 4); nm = substr(nm, 1, RSTART - 1) }
  if (match(nm, /^P[0-9]+[a-z]?/) || match(nm, /^S[0-9]+/)) lab = substr(nm, RSTART, RLENGTH)
  else lab = "-"
  k = p "|" part "|" lab "|" skel(nm)
  if (!(k in rows)) { key[++ng] = k }
  rows[k]++
  if (!((k, dev) in seen)) { seen[k, dev] = 1; nd[k]++ }
}
END {
  if (hdr) { printf "paper\tpart\tlabel\tskeleton\trows\tdev_set\n"; hdr = 0 }
  for (i = 1; i <= ng; i++) {
    k = key[i]; split(k, a, "|")
    printf "%s\t%s\t%s\t%s\t%d\t%d\n", a[1], a[2], a[3], a[4], rows[k], nd[k]
    tot += rows[k]
  }
  printf "GROUPS %d ROWS %d\n", ng, tot > "/dev/stderr"
}
