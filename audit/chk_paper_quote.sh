#!/usr/bin/env bash
# chk_paper_quote.sh —— 论文里【转录的实跑输出块】与当前实跑的分档对账（§176）
#
# 为什么要它：§159 抓到「三个 194 是同一个常量打三遍」；§176 又抓到
# `bsd` 的复现块写着 `RESULT PASS=12`（实跑 15）、`pvsnp` 写着 `PASS=26`（实跑 27）。
# 这类错【不会自己报】：论文照常渲染、代码照常跑，只有人眼对才发现。
#
# 三类被检查的行（都在 ``` 围栏内）：
#   A 转录行：以 [PASS] / [sanity] / RESULT / EXIT / k= / res= 开头 —— 必须能在
#     audit/outs/<p>.txt 的【同标签行】里找到（原文，或数字序列完全相同）（译文/改词）。
#   B 拼接行：`RESULT … EXIT=0` 写在一行 —— 先按标签切段，逐段按 A 判。
#   C 期望注记：以 # 开头且提到 RESULT PASS=n / N · M 小计 —— n 与两个小计必须是本篇
#     RESULT 行里真出现过的值。**这条是 §176 抓到 pvsnp PASS=26 的那一处。**
#
# 分档：exact（原文命中）/ num（同标签下数字序列相同 ⟹ 译文）/ drift（找不到 ⟹ 错）。
# 判定压在 drift：译文不影响复现，数字对不上就影响。
#
# 用法：bash audit/chk_paper_quote.sh      退出码 0 = drift 0 段
set -u
cd "$(dirname "$0")/.." || exit 2

OUTS=audit/outs
PROG='
function label(s) {
  if (match(s, /^(\[PASS\]|\[sanity\]|RESULT|EXIT|k=|res=)/)) return substr(s, RSTART, RLENGTH)
  return ""
}
function numseq(s,   out, k) {
  out = ""
  while (match(s, /[-+]?[0-9]+(\.[0-9]+)?([eE][-+]?[0-9]+)?/)) {
    out = out (out == "" ? "" : ",") substr(s, RSTART, RLENGTH)
    s = substr(s, RSTART + RLENGTH)
  }
  return out
}
BEGIN { ns = 0 ; nd = 0 ; nn = 0 }
FNR == 1 { stage++ }
stage == 1 {                      # ---- 读实跑输出 ----
  line = $0
  sub(/^[ \t]+/, "", line)
  L = label(line)
  if (L != "") { EX[L "\t" line] = 1 ; SQ[L "\t" numseq(line)] = 1 }
  if (L == "RESULT") {
    t = line
    while (match(t, /PASS=[0-9]+/)) { PV[substr(t, RSTART+5, RLENGTH-5)] = 1 ; t = substr(t, RSTART+RLENGTH) }
  }
  next
}
stage == 2 {                      # ---- 读论文 ----
  if ($0 ~ /^```/) { inq = !inq ; next }
  if (!inq) next
  line = $0
  if (line ~ /^#/) {              # C 期望注记
    t = line ; bad = 0
    while (match(t, /RESULT PASS=[0-9]+/)) {
      v = substr(t, RSTART + 12, RLENGTH - 12)
      if (!(v in PV)) { bad = 1 ; printf "    [DRIFT-C] %s: PASS=%s 不在本篇 RESULT 里\n", FILE, v }
      t = substr(t, RSTART + RLENGTH)
    }
    if (line ~ /(subtotals?|小计)/) {
      t = line
      while (match(t, /[0-9]+ *· *[0-9]+/)) {
        seg = substr(t, RSTART, RLENGTH) ; gsub(/·/, " ", seg) ; split(seg, pp, " ")
        for (i in pp) if (pp[i] != "" && !(pp[i] in PV)) { bad = 1 ; printf "    [DRIFT-C] %s: 小计 %s 不在本篇 RESULT 里\n", FILE, pp[i] }
        t = substr(t, RSTART + RLENGTH)
      }
    }
    if (bad) nd++
    next
  }
  sub(/^[ \t]+/, "", line)
  while (line != "") {            # A/B 转录行（B 会切成多段）
    L = label(line)
    if (L == "") break
    rest = substr(line, RLENGTH + 1)
    if (match(rest, /[ \t](RESULT|EXIT|k=|res=|\[PASS\]|\[sanity\])/)) {
      seg = L substr(rest, 1, RSTART - 1) ; line = substr(rest, RSTART + 1)
    } else { seg = L rest ; line = "" }
    gsub(/[ \t]+$/, "", seg)
    ns++
    if ((seg) in EX || ((L "\t" seg) in EX)) { ne++ ; continue }
    if ((L "\t" numseq(seg)) in SQ) { nn++ ; printf "    [num ] %s: %s\n", FILE, substr(seg, 1, 88) ; continue }
    nd++ ; printf "    [DRIFT] %s: %s\n", FILE, substr(seg, 1, 88)
  }
}
END { printf "%-22s %6d %6d %5d %5d\n", FILE, ns, ne, nn, nd }
'
printf '%-22s %6s %6s %5s %5s\n' "file" "segment" "exact" "num" "drift"
T=0; E=0; N=0; D=0
for pair in "bsd-conjecture:bsd" "hardy-littlewood:hl" "hodge-conjecture:hodge" \
            "navier-stokes:ns" "p-vs-np:pvsnp" "riemann-hypothesis:rh" "yang-mills:ym"; do
  dir=${pair%%:*}; key=${pair##*:}
  for lang in paper.zh.md paper.md; do
    f="$dir/$lang"; [ -f "$f" ] || continue
    out=$(awk -v FILE="$dir/$lang" "$PROG" "$OUTS/$key.txt" "$f")
    echo "$out" | tail -1
    echo "$out" | sed '$d'
    v=$(echo "$out" | tail -1)
    T=$((T + $(echo "$v" | awk '{print $2}')))
    E=$((E + $(echo "$v" | awk '{print $3}')))
    N=$((N + $(echo "$v" | awk '{print $4}')))
    D=$((D + $(echo "$v" | awk '{print $5}')))
  done
done
sed_clean() { :; }
echo
echo "合计：段 $T · 逐字 $E · 数字序列同 $N · 漂移 $D"
[ "$D" -eq 0 ] && { echo "EXIT=0"; exit 0; } || { echo "EXIT=1"; exit 1; }
