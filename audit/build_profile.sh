#!/usr/bin/env bash
# build_profile.sh —— 结构数据画像 v2 的唯一驱动
# 1) 逐篇把台账 line 列解析成地址状态（valid/blank/oob）
# 2) 用保序最大匹配把「输出侧」与「声明侧」对齐，输出 matrix/profile_v2.csv
# 全部为只读计算，不改台账、不改论文。
set -u
cd "$(dirname "$0")/.." || exit 1

# 0) 【必须先做】从当前实跑的 outs 重新生成 census2 —— 它是画像的「输出侧」。
#    旧版直接吃仓库里那份 census2.csv，它由一轮手工命令生成后就没人再生成 ⟹
#    代码按机制收紧 tol（§181 把 bsd 三条改成 1e-10/1.53e-10/8.0e-10）之后
#    census2 里仍是旧的 1e-08/1e-09 ⟹ 画像与台账对不上而无人报（X2 只自证 dev，
#    不证 tol）。这一步让「陈旧」在结构上不可能发生。
#    注意：census2.awk 只吐数据行、不吐表头（A1 用 NR==1 跳表头）⟹ 表头在这里补。
{ echo "paper,seq,src,name,dev,tol"; awk -f audit/census2.awk audit/outs/bsd.txt audit/outs/hl.txt \
    audit/outs/hodge.txt audit/outs/ns.txt audit/outs/pvsnp.txt audit/outs/rh.txt \
    audit/outs/ym.txt; } > audit/census2.csv

declare -A DIR=( [bsd]="bsd-conjecture" [hl]="hardy-littlewood" [hodge]="hodge-conjecture" \
                 [ns]="navier-stokes" [pvsnp]="p-vs-np" [rh]="riemann-hypothesis" [ym]="yang-mills" )
OUT=matrix/addr_state.csv
: > "$OUT"
for p in bsd hl hodge ns pvsnp rh ym; do
    f="${DIR[$p]}/paper.zh.md"
    [ -f "$f" ] || { echo "  缺 $f" >&2; exit 1; }
    awk -F, -f audit/addr_state.awk "$f" matrix/forcedness_ledger.csv | awk -F, -v p="$p" '$1==p'
done >> "$OUT"

awk -F, -v PHASEFILE=matrix/profile_phase.csv -f audit/profile_v2.awk audit/census2.csv "$OUT" matrix/forcedness_ledger.csv > matrix/profile_v2.csv

printf "[addr] 台账地址状态："
awk -F, '{c[$3]++} END{for(k in c) printf "%s=%d ", k, c[k]; print ""}' "$OUT"
printf "[profile] 行数（含表头）= %s\n" "$(wc -l < matrix/profile_v2.csv)"

# 画像 v3：把【面 / 零的种类 / 残差型】三条划分并进一张表（§184）。
# 与 v2 一样是【派生物】——每次都由 v2 + 台账现算，绝不手改。
awk -F, -f audit/profile_v3.awk matrix/profile_v2.csv matrix/forcedness_ledger.csv > matrix/profile_v3.csv
printf "[profile-v3] "
awk -F, '{f[$8]++; z[$9]++; r[$10]++}
  END{ printf "face["; for(k in f) printf "%s=%d ", k, f[k];
       printf "] zk[";  for(k in z) printf "%s=%d ", k, z[k];
       printf "] ra[";  for(k in r) printf "%s=%d ", k, r[k]; printf "]\n" }' matrix/profile_v3.csv
