#!/usr/bin/env awk -f
# §173 从【实跑输出】重建权威普查表：paper,name,dev,tol
# 三种已知格式都要吃得下：
#   A `PASS <name>  dev=X tol=Y`            （rh / ym / hodge 段二）
#   B `[PASS] <name>  dev=X tol=Y`          （pvsnp / ns 段二）
#   C `PASS <name>  dev=X (tol=Y)`          （bsd / hl）
#   D `PASS <name>`（无 dev 行，如 pvsnp 的 P6/P7/P8a-c）
BEGIN{ }
FNR==1{ split(FILENAME, parts, "/"); p=parts[length(parts)]; sub(/\.txt$/,"",p) }
{
  line=$0
  if (line ~ /^[ \t]*RESULT/) next
  if (line !~ /^[ \t]*(\[PASS\]|PASS|\*\*FAIL\*\*)/) next
  nm=line
  sub(/^[ \t]*(\[PASS\]|PASS|\*\*FAIL\*\*)[ \t]*/,"",nm)
  dev="NA"; tol="NA"
  if (match(nm, /dev=[-0-9.eE+]+[ \t]*\(?tol=[-0-9.eE+]+\)?/)) {
    pair=substr(nm, RSTART, RLENGTH)
    tmp=pair; sub(/^dev=/,"",tmp); sub(/[ \t].*/,"",tmp); dev=tmp
    tmp=pair; sub(/^.*tol=/,"",tmp); sub(/\).*$/,"",tmp); tol=tmp
    nm=substr(nm, 1, RSTART-1)
  }
  gsub(/[ \t]+$/,"",nm); gsub(/^[ \t]+/,"",nm); gsub(/,/,";",nm)
  printf "%s,%s,%s,%s\n", p, nm, dev, tol
}
