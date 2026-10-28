#!/usr/bin/env bash
# verify_ledger.sh —— 台账全轴单命令复算器（§175）
#
# 为什么要它：dev / tol / ratio / class / name / line 这几条轴，此前是【六支一次性脚本】
# （census.awk · reconcile.sh · fix_dev.awk · fix_tol.awk · fix_ratio.awk · chk_addr.sh）
# 分别查过的。一次性脚本的问题：**它查完就没了，下一刀改了台账没人再查。**
# 本支把六条轴合成一条命令，每轴给 PASS / FAIL。
#
# 两栏制（本支的骨架）：
#   【机器支撑】= 能从台账 + outs 现算出来，算错就翻红
#       X1 行数 · X2 (dev,tol) 对能在本篇 outs 同一行找到 · X3 ratio=dev/tol
#       X4 class 跟随 dev · X5 四条守恒式（class/invariance/dev/name）
#   【仅声明】= 推不出来，只有一个数在 README 里写着；本支只做一件事：
#       **把 README 的声明数字读出来与台账现算比对**（X6，读 README 里的 LEDGER-AXES 块）
#       ⟹ 声明可以假，但**声明与台账不许对不上**，这才是能翻红的部分。
#
# 用法：bash audit/verify_ledger.sh          退出码 0 = 全 PASS
set -u
cd "$(dirname "$0")/.." || exit 2

LED=${1:-matrix/forcedness_ledger.csv}
RD=${2:-README.md}
OUTS=audit/outs
fails=0

say_fail() { printf '  [FAIL] %-36s dev=%s\n' "$1" "$2"; fails=$((fails+1)); }
say_pass() { printf '  [PASS] %-36s dev=0.000000e+00\n' "$1"; }

# ---------- 一趟 awk：X1-X5 全部现算 ----------
readarray -t R < <(awk -F, -v OUTS="$OUTS" '
BEGIN{
  split("bsd hl hodge ns pvsnp rh ym", P, " ")
  for (pi in P) {
    p = P[pi]; f = OUTS "/" p ".txt"
    while ((getline ln < f) > 0) {
      d = ""; t = ""
      if (match(ln, /dev=[-0-9.eE+]+/)) { d = substr(ln, RSTART+4, RLENGTH-4) }
      if (match(ln, /tol=[-0-9.eE+]+/)) { t = substr(ln, RSTART+4, RLENGTH-4) }
      if (d != "" && t != "") SEEN[p "|" (d+0) "|" (t+0)] = 1
      if (d != "") SEEN_D[p "|" (d+0)] = 1
    }
    close(f)
  }
}
FNR == 1 { next }
{
  n++
  dev = $4; tol = $5; cls = $3
  zero = (dev == "NA(implicit)" || (dev+0) == 0) ? 1 : 0

  if ((cls == "FORCED") != (zero == 1)) badcls++

  if (dev == "NA(implicit)") { implied++ }
  else {
    tt = tol + 0; rr = $6 + 0
    if (tol == "NA(implicit)") badratio++
    else if (tt == 0) { if (rr != 0) badratio++; }
    else { c = (dev+0)/tt; d = (c>0)?(rr-c)/c:0; if (d<0) d=-d
           if (d >= 1e-6) badratio++ }
    # tol=NA ⟹ 过程面/读数行（TREND·FLOW）：不要求 tol 对得上，但仍要求 dev 在那篇输出里逐字出现
    if (tol == "NA") { if ((($1 "|" (dev+0)) in SEEN_D)) pair_ok++; else { pair_bad++; if (pair_show < 5) { printf "    pair-miss %s %s dev=%s tol=NA(仅校 dev)\n", $1, $2, dev > "/dev/stderr"; pair_show++ } } }
    else if ((($1 "|" (dev+0) "|" (tt+0)) in SEEN)) pair_ok++
    else { pair_bad++; if (pair_show < 5) { printf "    pair-miss %s %s dev=%s tol=%s\n", $1, $2, dev, tol > "/dev/stderr"; pair_show++ } }
  }

  c_cls[$3]++; c_inv[$8]++
  if (dev == "NA(implicit)") {}
  else if ((dev+0) == 0) { if ((tol+0) == 0) devz_ex++; else devz_sl++ }
  else devz_re++
  if ($9 ~ /points to/) c_addr++
  else if ($9 == "" || $9 ~ /^NA\(/) c_noobj++
  else c_named++
  if ($9 == "" || $9 ~ /^NA\(/) unnamed++
}
END{
  printf "SUMMARY %d %d %d %d %d %d\n", n, pair_ok, pair_bad, badratio, badcls, unnamed
  for (k in c_cls) printf "CLS %s %d\n", k, c_cls[k]
  for (k in c_inv) printf "INV %s %d\n", k, c_inv[k]
  printf "DEVZ %d %d %d %d\n", devz_ex, devz_sl, devz_re, implied
  # 规范序：与 README 的 LEDGER-AXES 块必须逐位同序，否则比的是排序不是对象
  printf "CLSCSV %d,%d,%d,%d,%d,%d\n", c_cls["FORCED"], c_cls["FORCED_eps"], c_cls["FORCED_trunc"], c_cls["FLOW"], c_cls["DUP"], c_cls["TREND"]
  printf "INVCSV %d,%d,%d,%d,%d,%d\n", c_inv["sample"], c_inv["theorem"], c_inv["cited"], c_inv["dup"], c_inv["na"], c_inv["WALL"]
  printf "DEVCSV %d,%d,%d,%d\n", devz_ex, devz_sl, devz_re, implied
  printf "NAMECSV %d,%d,%d\n", c_named, c_addr, c_noobj
}' "$LED" 2>/tmp/vl_pair.txt)

summary=$(printf '%s\n' "${R[@]}" | awk '$1=="SUMMARY"{print $2,$3,$4,$5,$6,$7}')
n=$(echo "$summary"      | cut -d' ' -f1)
pair_ok=$(echo "$summary"| cut -d' ' -f2)
pair_bad=$(echo "$summary"|cut -d' ' -f3)
badratio=$(echo "$summary"|cut -d' ' -f4)
badcls=$(echo "$summary" | cut -d' ' -f5)
unnamed=$(echo "$summary"| cut -d' ' -f6)
implied=$(( n - pair_ok - pair_bad ))
devz=$(printf '%s\n' "${R[@]}" | awk '$1=="DEVZ"{printf "%d,%d,%d,%d", $2,$3,$4,$5}')
cls_csv=$(printf '%s\n' "${R[@]}" | awk '$1=="CLSCSV"{print $2}')
inv_csv=$(printf '%s\n' "${R[@]}" | awk '$1=="INVCSV"{print $2}')
devz_csv=$(printf '%s\n' "${R[@]}" | awk '$1=="DEVCSV"{print $2}')
name_csv=$(printf '%s\n' "${R[@]}" | awk '$1=="NAMECSV"{print $2}')
c_named=$(echo "$name_csv" | cut -d, -f1); c_addr=$(echo "$name_csv" | cut -d, -f2); c_noobj=$(echo "$name_csv" | cut -d, -f3)
# 两条新轴（§177）：组数（按 篇/块/标签/骨架 归组的判据条数）与锚表强制/多义
grp_n=$(for p in bsd hl hodge ns pvsnp rh ym; do awk -f "$(dirname "$0")/group.awk" "$OUTS/$p.txt" 2>/dev/null; done | awk -F'\t' 'NR>1 && $1!="paper"{g++} END{print g+0}')
anc_csv=$(awk -F, 'NR>1{if($5=="forced")f++; else if($5=="ambiguous")a++} END{printf "%d,%d", f+0, a+0}' matrix/ledger_anchor.csv 2>/dev/null)
cls_sum=$(printf '%s\n' "${R[@]}" | awk '$1=="CLS"{s+=$3}END{print s+0}')
inv_sum=$(printf '%s\n' "${R[@]}" | awk '$1=="INV"{s+=$3}END{print s+0}')
devz_sum=$(echo "$devz" | awk -F, '{print $1+$2+$3+$4}')

echo "=========================================================="
echo " 台账全轴复算（容器 = $LED · 空间 = 6 条机器轴 + 1 条声明轴）"
echo "=========================================================="
echo
echo "[机器栏] 现算"
printf '  行数            %s\n' "$n"
printf '  (dev,tol) 对    命中 %s · 未命中 %s · 无 dev=（implicit）%s\n' "$pair_ok" "$pair_bad" "$implied"
printf '  ratio=dev/tol  不一致 %s\n' "$badratio"
printf '  class 跟随 dev  违反 %s\n' "$badcls"
printf '  class 五档      %s ⟹ %s\n' "$cls_csv" "$cls_sum"
printf '  invariance      %s ⟹ %s\n' "$inv_csv" "$inv_sum"
printf '  dev 四类        %s ⟹ %s\n' "$devz" "$devz_sum"
printf '  name 三分       有名 %s · 地址落空 %s · 无对象 %s\n' "$c_named" "$c_addr" "$c_noobj"
printf '  §177 判据组       %s 组 / %s 行 ⟹ 参数点重复 %s\n' "$grp_n" "$n" "$((n-grp_n))"
printf '  §177 锚表        %s（forced,ambiguous）\n' "$anc_csv"
sed -n '1,5p' /tmp/vl_pair.txt 2>/dev/null

echo
echo "[声明栏] 读 README 的 LEDGER-AXES 块，与上面现算逐项比"
# 取块内某键的值：只认行首的「键=」，块以 <!-- LEDGER-AXES 起、--> 止
decl_get() {
  awk -v k="$1" '/<!-- LEDGER-AXES/{f=1;next} /-->/{if(f)exit} f { if (index($0, k "=") == 1) { v=substr($0, length(k)+2); gsub(/\r/,"",v); print v; exit } }' "$RD"
}
cmp_decl() {  # $1 键 $2 现算
  local k="$1" got="$2" want
  want=$(decl_get "$k")
  if [ -z "$want" ]; then say_fail "X6 README 未声明 $k" 1; return; fi
  if [ "$want" = "$got" ]; then say_pass "X6 $k = $got（README 与台账一致）"
  else say_fail "X6 $k：README=$want 台账=$got" 1; fi
}
cmp_decl rows       "$n"
cmp_decl class      "$cls_csv"
cmp_decl invariance "$inv_csv"
cmp_decl dev        "$devz_csv"
cmp_decl name       "$name_csv"
cmp_decl groups     "$grp_n"
cmp_decl anchor     "$anc_csv"
# line 轴推不出来（§172：180 个不同地址里 68 个落空）⟹ 只查 README 里的声明还在，不查值
line_decl=$(decl_get line)
if [ -n "$line_decl" ]; then say_pass "X6 line 已声明 = $line_decl（推不出来，见 §172/§175）"
else say_fail "X6 README 未声明 line" 1; fi

echo
echo "[判定] 机器栏"
[ "$n" = "198" ] && say_pass "X1 行数 = 198" || say_fail "X1 行数 = 198" "$((n-198))"
[ "$pair_bad" = "0" ] && say_pass "X2 (dev,tol) 对全部可在 outs 找到" || say_fail "X2 (dev,tol) 对全部可在 outs 找到" "$pair_bad"
[ "$badratio" = "0" ] && say_pass "X3 ratio 全部 = dev/tol" || say_fail "X3 ratio 全部 = dev/tol" "$badratio"
[ "$badcls" = "0" ] && say_pass "X4 class 全部跟随 dev" || say_fail "X4 class 全部跟随 dev" "$badcls"
[ "$cls_sum" = "$n" ] && say_pass "X5a class 五档合计 = 行数" || say_fail "X5a class 五档合计 = 行数" "$((cls_sum-n))"
[ "$inv_sum" = "$n" ] && say_pass "X5b invariance 合计 = 行数" || say_fail "X5b invariance 合计 = 行数" "$((inv_sum-n))"
[ "$devz_sum" = "$n" ] && say_pass "X5c dev 四类合计 = 行数" || say_fail "X5c dev 四类合计 = 行数" "$((devz_sum-n))"
[ "$unnamed" = "106" ] && say_pass "X5d 无名判据 = 106（三分：有名 $c_named / 地址落空 $c_addr / 无对象 $c_noobj）" || say_fail "X5d 无名判据 = 106" "$((unnamed-106))"
[ "$((c_named + c_addr + c_noobj))" = "$n" ] && say_pass "X5e name 三分合计 = 行数" || say_fail "X5e name 三分合计 = 行数" "$((c_named + c_addr + c_noobj - n))"
[ "$((n - grp_n))" = "15" ] && say_pass "X5f §177 参数点重复 = 198 − $grp_n = 15 行" || say_fail "X5f §177 参数点重复应为 15" "$((n - grp_n - 15))"
[ -n "$anc_csv" ] && say_pass "X5g §177 锚表可读 = $anc_csv" || say_fail "X5g §177 锚表不可读" 1

echo
echo "[判定] 锚表自证（§177，转 audit/verify_anchor.sh）"
if [ -f "$(dirname "$0")/verify_anchor.sh" ]; then
  an_out=$(bash "$(dirname "$0")/verify_anchor.sh" 2>&1); an_rc=$?
  printf '%s\n' "$an_out" | grep -E '^\s+\[(PASS|FAIL|READ)\]' | sed 's/^/ /'
  [ "$an_rc" -eq 0 ] && say_pass "X6 §177 锚表 A1–A6 全绿" || say_fail "X6 §177 锚表有 FAIL" "$an_rc"
else say_fail "X6 缺 audit/verify_anchor.sh" 1; fi

echo
echo "[判定] 画像 v2（§179：保序最大匹配 + 判据文本三分）"
# 作者确认: 2026-10-29 用户令「完成所有解决不了的遗留问题」（AUDIT §179）——
# 本段只【新增】X7a–X7d 四条判项（画像行数 / 两个三分守恒 / 可复现），
# 既有 X1–X6 的判定对象、阈值与口径一字未动。
PROF=matrix/profile_v2.csv
prof_n=$(awk 'NR>1' "$PROF" 2>/dev/null | wc -l | tr -d ' ')
prof_join=$(awk -F, 'NR>1{c[$7]++} END{printf "%d,%d,%d", c["pair"]+0, c["lcs"]+0, c["none"]+0}' "$PROF" 2>/dev/null)
prof_crit=$(awk -F, 'NR>1{
  s = $13; for (i = 14; i <= NF; i++) s = s "," $i;   # criterion 文本里可能含逗号：从 13 列起拼回整串
  if (s ~ /^paper:/) t = "p";
  else if (s ~ /^NA\(/) t = "d";
  else t = "n";
  c[t]++ } END{printf "%d,%d,%d", c["n"]+0, c["p"]+0, c["d"]+0}' "$PROF" 2>/dev/null)
prof_pap=$(awk -F, 'NR>1{p[$1]=1} END{n=0; for(k in p) n++; printf "%d", n}' "$PROF" 2>/dev/null)
cmp_decl profile   "$prof_join"
cmp_decl criterion "$prof_crit"
prof_sum=$(echo "$prof_join" | awk -F, '{print $1+$2+$3}')
if [ "$prof_n" = "198" ] && [ "$prof_sum" = "198" ]
then say_pass "X7a 画像 198 行 · join 三分合计 = 198（pair,lcs,none = $prof_join）"
else say_fail "X7a 画像行数或 join 三分不守恒" "$((198-prof_n))"; fi
prof_csum=$(echo "$prof_crit" | awk -F, '{print $1+$2+$3}')
if [ "$prof_csum" = "198" ]
then say_pass "X7b 判据文本三分合计 = 198（sigma/paper/NA_dead = $prof_crit）"
else say_fail "X7b 判据文本三分不守恒" "$((prof_csum-198))"; fi
if [ "$prof_pap" = "7" ]; then say_pass "X7c 画像覆盖 7 篇"
else say_fail "X7c 画像覆盖篇数 = $prof_pap" "$((7-prof_pap))"; fi
md5_now=$(md5sum "$PROF" 2>/dev/null | cut -d' ' -f1 | cut -c1-12)
md5_decl=$(decl_get profile_md5)
if [ -n "$md5_now" ] && [ "$md5_now" = "$md5_decl" ]
then say_pass "X7d 画像指纹 = README 声明（$md5_now）"
else say_fail "X7d 画像指纹与 README 不符（now=${md5_now:-空} decl=${md5_decl:-未声明}）" 1; fi
bash "$(dirname "$0")/build_profile.sh" >/dev/null 2>&1
md5_after=$(md5sum "$PROF" 2>/dev/null | cut -d' ' -f1 | cut -c1-12)
if [ "$md5_now" = "$md5_after" ]
then say_pass "X7e 画像可复现（重建 md5 不变）"
else say_fail "X7e 画像不可复现（$md5_now → $md5_after）" 1; fi
# X7f：相位偏移（自由参数）的声明与现算逐位比 —— 声明格式 cnt,max（非零相位的篇数, 最大相位）
PHF=matrix/profile_phase.csv
ph_decl=$(decl_get phase)
ph_now=$(awk -F, 'NR>1{ if($2+0!=0){ n++; if($2+0>m) m=$2+0 } } END{ printf "%d,%d", n, m }' "$PHF" 2>/dev/null)
if [ -n "$ph_decl" ] && [ "$ph_decl" = "$ph_now" ]
then say_pass "X7f 相位声明 = 现算（${ph_now}；产物 $PHF）"
else say_fail "X7f 相位声明与现算不符（now=${ph_now:-空} decl=${ph_decl:-未声明}）" 1; fi

# X9：台账 tol 与输出 tol 在【两侧都是数值】的行上必须相等（对象面不许背离）
# 作者确认(2026-10-30): 只【新增】X9 一条判项；X1–X8 的判定对象、阈值与口径一字未动。
# 依据（AUDIT §184）：census2 陈旧曾让 bsd 三条「代码按 1e-10 判、台账记 1e-08」藏了一整轮，
# 而 X2 只自证 dev、不证 tol ⟹ 必须有一条判据专门盯两侧 tol 是否是同一个界。
x9_out=$(awk -F, -f audit/x9_tol.awk audit/census2.csv "$LED" "$PROF" 2>/dev/null)
x9_bad=$(printf '%s\n' "$x9_out" | awk '/^BAD/{print $2}')
x9_tot=$(printf '%s\n' "$x9_out" | awk '/^BAD/{print $4}')
if [ "${x9_bad:-1}" = "0" ]
then say_pass "X9 台账 tol = 输出 tol（两侧都是数值的 $x9_tot 行，对象面 0 处背离）"
else say_fail "X9 台账 tol ≠ 输出 tol（$x9_bad 处）" "$x9_bad"; fi

# X10：画像 v3 的三条划分（面 / 零的种类 / 残差型）都必须把同一 198 行分完，
#      且「残差型 R 的条数」必须等于 C 探针 resid_type 自己印出来的条数（双路交叉）。
# 作者确认(2026-10-30): 只【新增】X10；X1–X9 的判定对象、阈值与口径一字未动。
P3=matrix/profile_v3.csv
p3_n=$(awk 'END{print NR}' "$P3" 2>/dev/null)
p3_sums=$(awk -F, '{f[$8]++; z[$9]++; r[$10]++}
  END{ for(k in f) s1+=f[k]; for(k in z) s2+=z[k]; for(k in r) s3+=r[k];
       printf "%d,%d,%d", s1+0, s2+0, s3+0 }' "$P3" 2>/dev/null)
p3_face=$(awk -F, '{c[$8]++} END{printf "%d,%d,%d", c["object"]+0, c["process"]+0, c["dup"]+0}' "$P3" 2>/dev/null)
p3_ra=$(awk -F, '{c[$10]++} END{printf "%d,%d,%d", c["R"]+0, c["A"]+0, c["X"]+0}' "$P3" 2>/dev/null)
if [ "$p3_n" = "198" ] && [ "$p3_sums" = "198,198,198" ]
then say_pass "X10a 画像 v3 三条划分各自加回 198（face=$p3_face · ra=$p3_ra）"
else say_fail "X10a 画像 v3 行数或三条划分不守恒（rows=${p3_n:-空} sums=${p3_sums:-空}）" 1; fi
rt_out=${RESID_OUT:-../dpgo/cx/resid_type_out.txt}
[ -f "$rt_out" ] || rt_out=/d/dpgo/cx/resid_type_out.txt
rt_r=$(grep -o 'R 型（.*）= [0-9]* / 60' "$rt_out" 2>/dev/null | sed 's/.*= \([0-9]*\) .*/\1/' | head -1)
p3_r=$(echo "$p3_ra" | cut -d, -f1)
if [ -n "$rt_r" ] && [ "$rt_r" = "$p3_r" ]
then say_pass "X10b 画像 v3 的 R 条数 = C 探针 resid_type 的 R 条数（$p3_r；双路交叉）"
else say_fail "X10b R 条数两路不等（画像=${p3_r:-空} 探针=${rt_r:-空}）" 1; fi
# X10c：README 里声明的 v3 三条划分与两侧 tol 一致行数，必须与现算逐位相同
cmp_decl v3_face "$p3_face"
p3_zk=$(awk -F, '{c[$9]++} END{printf "%d,%d,%d,%d", c["geo"]+0, c["count"]+0, c["implicit"]+0, c["nonzero"]+0}' "$P3" 2>/dev/null)
cmp_decl v3_zk "$p3_zk"
cmp_decl v3_ra "$p3_ra"
cmp_decl tol_agree "${x9_tot:-0}"

# X8：`tol` 有数值 ⟺ 该行在【对象面】（2026-10-29 · §179.7 用户定案 · 判据在 audit/x8_face.awk）
#   对象面 = FORCED / FORCED_eps / FORCED_trunc  ⟹ tol 必须数值
#   过程面 = FLOW / TREND                        ⟹ tol 必须 NA
#   DUP                                          ⟹ 重复标记，不属两面，只计数
# 作者确认: 用户 2026-10-29 定案「tol 有值 ⟺ 对象面」，本次为【新增】X8 三子项，不动任何旧判据（G22）
X8_LINE=$(awk -f "$(dirname "$0")/x8_face.awk" "$LED" 2>/tmp/vl_x8.err)
read -r _ x8_obj x8_objn x8_prc x8_prcn x8_dup x8_unk x8_bad <<< "$X8_LINE"
x8_rows=$(( $(wc -l < "$LED") - 1 ))
x8_cov=$(( x8_obj + x8_prc + x8_dup ))
if [ "$x8_objn" = "$x8_obj" ]
then say_pass "X8a 对象面 tol 齐全（$x8_objn/$x8_obj 行有数值）"
else say_fail "X8a 对象面缺 tol（$x8_objn/$x8_obj）" "$((x8_obj-x8_objn))"; fi
if [ "$x8_prcn" = "$x8_prc" ]
then say_pass "X8b 过程面 tol 全 NA（$x8_prcn/$x8_prc 行）"
else say_fail "X8b 过程面带 tol（$x8_prcn/$x8_prc 行是 NA）" "$((x8_prc-x8_prcn))"; fi
if [ "$x8_unk" = "0" ] && [ "$x8_bad" = "0" ]
then say_pass "X8c 无未归面桶 / 无违规行（bad=$x8_bad · unk=$x8_unk）"
else say_fail "X8c 未归面桶或违规行（bad=$x8_bad unk=$x8_unk）" "$((x8_bad+x8_unk))"; fi
if [ "$x8_cov" = "$x8_rows" ]
then say_pass "X8d 三分覆盖 $x8_obj+$x8_prc+$x8_dup = $x8_rows 行（缺口 $((x8_rows-x8_cov))）"
else say_fail "X8d 三分不覆盖（$x8_cov ≠ $x8_rows）" "$((x8_rows-x8_cov))"; fi

echo
echo "RESULT FAIL=$fails"
if [ "$fails" -eq 0 ]; then echo "EXIT=0"; exit 0; else echo "EXIT=1"; exit 1; fi
