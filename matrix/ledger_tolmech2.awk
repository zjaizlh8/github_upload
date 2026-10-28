# ledger_tolmech2.awk — 台账第二轮：机制界回写 + FLOW 清 tol + rank 趋零降 TREND（2026-10-29）
# 用法：awk -F',' -f ledger_tolmech2.awk forcedness_ledger.csv > out
# 只按 $1:$2 精确命中，命中不在表里就报 MISS（防「改了 0 行还说改完」）。
BEGIN {
  OFS=","; nchg=0; nmiss=0
  # 键 = 行号标识；值 = "新class|新tol|新criterion后缀"
  T["bsd:42"]="|1e-10|〔tol 由机制导出：Ω 积分累加界 ε·NPTS·max|f|·4h/3 ≈ 1e-10（弃圆整 1e-8）〕"
  T["bsd:72"]="|1e-10|〔同 bsd:42〕"
  T["bsd:44"]="|1.53e-10|〔tol=(16·½·10⁻¹⁰+2·1e-10)/Ω=1.53e-10：锚 10 位截断 + Ω 累加界，两源独立相加〕"
  T["bsd:73"]="|1.53e-10|〔同 bsd:44〕"
  T["bsd:46"]="|8.0e-10|〔残差预测 = 16·(L真−L锚) ⟹ 上界 16·½·10⁻¹⁰=8.0e-10；实测 4.316e-10 与预测 4.315e-10 四位吻合 ⟹ 残差全来自锚截断〕"
  T["bsd:74"]="|8.0e-10|〔同 bsd:46〕"
  T["bsd:50#dup(sep)"]="TREND|NA|〔rank≥1 的「Abel 和趋零」是过程面（走到哪儿）⟹ 降 TREND；结构面由 part3 整数域判据承担（S1）〕"
  T["bsd:51#dup(sep)"]="TREND|NA|〔同 bsd:50〕"
  seen["bsd:42"]=0; seen["bsd:72"]=0; seen["bsd:44"]=0; seen["bsd:73"]=0
  seen["bsd:46"]=0; seen["bsd:74"]=0; seen["bsd:50#dup(sep)"]=0; seen["bsd:51#dup(sep)"]=0
}
NR==1 { print; next }
{
  key=$1":"$2
  fl=$3
  # FLOW = 有限窗流体量（过程面）⟹ tol 无结构含义，一律 NA
  if ($3=="FLOW") {
    $5="NA"; $6="NA"; $9=$9 " 〔FLOW 是过程面读数（S1），tol 无结构含义 ⟹ 清 NA〕"; fl="FLOW"
    print; nchg++; next
  }
  if (key in T) {
    seen[key]=1
    split(T[key], v, "|")
    if (v[1]!="") $3=v[1]
    $5=v[2]; $6="NA"
    if ($4+0 > 0 && v[2] != "NA") $6=sprintf("%.6g", ($4+0)/(v[2]+0))
    $9=$9 " " v[3]
    print; nchg++; next
  }
  print
}
END {
  for (k in seen) if (seen[k]==0) { printf "MISS %s（键不在表里）\n", k > "/dev/stderr"; nmiss++ }
  printf "CHANGED=%d MISS=%d\n", nchg, nmiss > "/dev/stderr"
}
