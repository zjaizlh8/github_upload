#!/usr/bin/env bash
# verify_anchor.sh —— §177 判决行锚表的自证与三条新轴（A1–A6）
#
# 背景：台账第 2 列 `line` 指的是【论文 Markdown 的行号】，随写作而腐（§172：180 个不同地址
# 里 68 个落空）。本支建的是另一根锚：**验证器输出的物理行号**（audit/census2.csv 的 `src`），
# 只要二进制不变它就不动；代码一改行号变，正好是该重核账的信号。
#
# A1 census2 行数 198 且逐篇 = RESULT PASS
# A2 **自证**：census2 每行的 src 行，重新从 outs 读出来，必须仍以判决词开头、
#    且含本行的 dev 字面（或 dev=NA）、且名字前 6 个字符一致 ⟹ 表能证明自己（本条是最硬的）
# A3 `paper#seq` 唯一标识 198 个不重复
# A4 判据组数：198 行按 (篇,块,标签,骨架) 归组 = 183 ⟹ 参数点重复 15 行（能翻红）
# A5 同号不同物：同篇同一 Pn 落在 ≥2 条判决行且骨架不同 —— **rh 必须 0**（§177.4 已改 S 前缀），
#    其余六篇照实报读数（不假装已修）
# A6 锚表：198 行 = forced 61 + ambiguous 137，且每个 forced 行的 dev 与 census 该行逐字相同
#
# 用法：bash audit/verify_anchor.sh      退出码 0 = 全 PASS
set -u
cd "$(dirname "$0")/.." || exit 2

fails=0
pass() { printf '  [PASS] %s  dev=0.000000e+00\n' "$1"; }
fail() { printf '  [FAIL] %s  dev=%d\n' "$1" "$2"; fails=$((fails + 1)); }

# ---------- A1 / A2 / A3 / A5：单支 awk 过 census2 与七份 outs ----------
ANCHOR_AWK='
function skel(s,   o) {
  o = s; gsub(/[0-9]/, "", o); gsub(/[=+*\/()（）·、,;:.\-—–_{}^$|]/, " ", o)
  gsub(/ +/, " ", o); sub(/^ /, "", o); sub(/ $/, "", o)
  if (o == "") o = "(no-name)"
  return o
}
BEGIN { FS = "," }
NR == 1 { next }
{
  paper = $1; seq = $2; src = $3; nm = $4; dev = $5; tol = $6
  n++
  per[paper]++
  id = paper "#" seq
  if (id in seenId) dupId++; else seenId[id] = 1
  # ---- A2 自证：重读 outs 的 src 行 ----
  f = "audit/outs/" paper ".txt"
  line = ""
  if ((f, src) in LINECACHE) line = LINECACHE[f, src]
  else {
    c = 0
    while ((getline l < f) > 0) { c++; if (c == src) { line = l; break } }
    close(f)
    LINECACHE[f, src] = line
  }
  if (line == "") { badsrc++; if (badsrc <= 3) printf "    [A2] %s src=%s 取不到行\n", paper, src }
  else {
    if (line !~ /^[ \t]*(\[PASS\]|PASS|\*\*FAIL\*\*)/) { badsrc++; if (badsrc <= 3) printf "    [A2] %s src=%s 不是判决行\n", paper, src }
    else if (dev != "NA" && index(line, "dev=" dev) == 0) { badev++; if (badev <= 3) printf "    [A2] %s src=%s 不含 dev=%s\n", paper, src, dev }
    else if (substr(nm, 1, 6) != "" && index(line, substr(nm, 1, 6)) == 0) { badnm++; if (badnm <= 3) printf "    [A2] %s src=%s 名字前6字不符\n", paper, src }
  }
  # ---- A5 撞号：同篇同标签、骨架不同 ----
  if (match(nm, /^P[0-9]+[a-z]?/)) {
    lab = substr(nm, RSTART, RLENGTH); k = paper "|" lab
    sk = skel(nm)
    if (!((k, sk) in seenSk)) { seenSk[k, sk] = 1; nsk[k]++ ; names[k] = names[k] " | " substr(nm, 1, 26) }
  }
}
END {
  printf "ROWS %d DUPID %d BADSRC %d BADEV %d BADNM %d\n", n, dupId + 0, badsrc + 0, badev + 0, badnm + 0
  for (k in nsk) { split(k, a, "|"); if (nsk[k] > 1) { coll++; COLL[a[1]]++; if (nsk[k] > 1 && COLL[a[1]] <= 2) printf "    [A5] %s %s 骨架 %d 种%s\n", a[1], a[2], nsk[k], names[k] } }
  printf "COLL %d RHCOLL %d\n", coll + 0, COLL["rh"] + 0
  printf "PER"; for (p in per) printf " %s=%d", p, per[p]; print ""
}'

out=$(awk "$ANCHOR_AWK" audit/census2.csv)
rows=$(echo "$out" | awk '/^ROWS/{print $2}')
dupid=$(echo "$out" | awk '/^ROWS/{print $4}')
badsrc=$(echo "$out" | awk '/^ROWS/{print $6}')
badev=$(echo "$out" | awk '/^ROWS/{print $8}')
badnm=$(echo "$out" | awk '/^ROWS/{print $10}')
coll=$(echo "$out" | awk '/^COLL/{print $2}')
rhcoll=$(echo "$out" | awk '/^COLL/{print $4}')

echo "  ── 锚表自证 ──"
[ "$rows" = "198" ] && pass "A1 census2 行数 = 198" || fail "A1 census2 行数 = 198" "$((rows - 198))"
[ "$dupid" = "0" ] && pass "A3 paper#seq 唯一（198 个标识）" || fail "A3 paper#seq 有重复" "$dupid"
s=$((badsrc + badev + badnm))
[ "$s" = "0" ] && pass "A2 自证：198 行 src 重读后仍是判决行且 dev/名字一致" || fail "A2 自证失败行数" "$s"

echo "  ── 两条新轴（读数 + 一条能翻红）──"
grp=$(for p in bsd hl hodge ns pvsnp rh ym; do awk -f audit/group.awk "audit/outs/$p.txt" 2>/dev/null; done | awk -F'\t' 'NR>1 && $1!="paper"{g++; s+=$5} END{printf "%d %d", g, s}')
gn=$(echo "$grp" | awk '{print $1}'); gs=$(echo "$grp" | awk '{print $2}')
[ "$gs" = "198" ] && pass "A4 组表覆盖 198 行（组数 $gn）" || fail "A4 组表行数与 198 不符" "$((gs - 198))"
[ "$((198 - gn))" = "15" ] && pass "A4 参数点重复 = 198 − $gn = 15 行" || fail "A4 参数点重复应为 15" "$((198 - gn - 15))"
[ "$rhcoll" = "0" ] && pass "A5 rh 同号不同物 = 0（已改 S 前缀）" || fail "A5 rh 仍有同号不同物" "$rhcoll"
printf '  [READ] A5 其余六篇同号不同物 = %d 组（照实报读数，未修）\n' "$coll"

echo "  ── 锚表强制匹配 ──"
an=matrix/ledger_anchor.csv
[ -f "$an" ] || { fail "A6 缺 matrix/ledger_anchor.csv" 1; an=/dev/null; }
st=$(awk -F, 'NR>1{if($5=="forced")f++; else if($5=="ambiguous")a++; else n++} END{printf "%d %d %d", f+0, a+0, n+0}' "$an" 2>/dev/null || echo "0 0 0")
f=$(echo "$st" | awk '{print $1}'); a=$(echo "$st" | awk '{print $2}'); nn=$(echo "$st" | awk '{print $3}')
[ "$((f + a + nn))" = "198" ] && pass "A6 锚表 198 行 = 强制 $f + 多义 $a" || fail "A6 锚表行数" "$((f + a + nn - 198))"
# 强制行必须真的读数相同：anchor 的第 2 列是台账行号、第 6 列是 census 的 seq
dv=$(awk -F, '
  FNR == 1 { stage++; next }
  stage == 1 { CN[$1 "|" $2] = $5; next }                 # census2: paper|seq -> dev
  stage == 2 { LG[$1 "|" $2] = $4; next }                 # ledger:  paper|line -> dev
  stage == 3 {
    if ($5 != "forced") next
    n++
    a = LG[$1 "|" $2]; b = CN[$1 "|" $6]
    if (a == b) ok++; else { bad++; if (bad <= 3) printf "  [A6] %s:%s 台账 dev=%s ≠ census dev=%s\n", $1, $2, a, b }
  }
  END { printf "A6SUM %d %d\n", ok + 0, bad + 0 }
' audit/census2.csv matrix/forcedness_ledger.csv "$an" 2>/dev/null)
dv=$(echo "$dv" | awk '/A6SUM/{print $2, $3}')
dok=$(echo "$dv" | awk '{print $1}'); dbad=$(echo "$dv" | awk '{print $2}')
[ "$dbad" = "0" ] && pass "A6 强制 $dok 行的 dev 与 census 逐字相同" || fail "A6 强制行 dev 不一致" "$dbad"

echo
echo "RESULT FAIL=$fails"
[ "$fails" -eq 0 ] && { echo "EXIT=0"; exit 0; } || { echo "EXIT=1"; exit 1; }
