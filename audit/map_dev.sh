#!/usr/bin/env bash
# §173c 台账 dev → 实跑原值 的映射（只认【唯一】匹配：同一篇内，按 4 位有效数字对齐）
set -u
cd "$(dirname "$0")/.." || exit 1
for p in rh ym hl pvsnp; do
  echo "=== $p"
  awk -F, -v P="$p" 'NR>1&&$1==P{  # 台账
      v=$4+0;
      if (v==0) next
      key=sprintf("%.3e",v); led[key]=led[key] " " $2
    }
    END{ for(k in led) print "L\t" k "\t" substr(led[k],2) }' matrix/forcedness_ledger.csv | sort > /tmp/lm_$p.txt
  awk -F, -v P="$p" '$1==P && $3!="NA"{v=$3+0; if(v==0) next; printf "C\t%.3e\t%.9e\t%s\n", v, v, $2}' audit/census.csv | sort -u > /tmp/cm_$p.txt
  # 每个台账 key 在普查里出现的行数
  awk -F'\t' 'NR==FNR{if($1=="C"){n[$2]++; val[$2]=val[$2]" "$3}} {next}' /tmp/cm_$p.txt /dev/null 2>/dev/null
  cut -f2 /tmp/cm_$p.txt | sort | uniq -c | awk '{print $2"\t"$1"\t"$3}' > /tmp/cn_$p.txt
  join -t$'\t' -1 2 -2 1 /tmp/lm_$p.txt /tmp/cn_$p.txt 2>/dev/null | awk -F'\t' '{printf "  %s -> %s  (普查同键 %s 条)  台账行: %s\n", $1, $5, $3, $2}'
done
