#!/usr/bin/env bash
# run7.sh —— Makefile `run` 的等价复现（本机无 make）
# 编译旗标与 Makefile 逐字一致：-Wall -Wextra -O2 -D_GNU_SOURCE + -lm
set -u
cd "$(dirname "$0")" || exit 2
CFLAGS="-Wall -Wextra -O2 -D_GNU_SOURCE"
mkdir -p bin audit/outs
fail=0
for src in bsd-conjecture/code/verify_bsd.c \
           hardy-littlewood/code/verify_hl.c \
           hodge-conjecture/code/verify_hodge.c \
           navier-stokes/code/verify_ns.c \
           p-vs-np/code/verify_pvsnp.c \
           riemann-hypothesis/code/verify_rh.c \
           yang-mills/code/verify_ym.c ; do
  name=$(basename "$src" .c); p=${name#verify_}
  out_bin="bin/$name"
  if ! gcc $CFLAGS -o "$out_bin" "$src" -lm 2>"/tmp/build_$p.log"; then
    echo "$src  BUILD_FAIL"; fail=1; continue
  fi
  warn=$(wc -l < "/tmp/build_$p.log")
  if out=$("./$out_bin" 2>&1); then
    printf '%s' "$out" > "audit/outs/$p.txt"
    res=$(printf '%s\n' "$out" | grep -oE 'RESULT PASS=[0-9]+ FAIL=[0-9]+' | tail -1)
    printf '%-34s warn=%s %s OK\n' "$src" "$warn" "$res"
  else
    printf '%-34s warn=%s EXIT!=0\n' "$src" "$warn"; fail=1
  fi
done
echo "RUN_FAIL=$fail"
exit $fail
