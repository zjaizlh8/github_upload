# 一键复算 / One-command recompute
#
#   make                 编译 7 个验证器到 bin/         (7 verifiers)
#   make run             编译 + 逐支实跑 7 个验证器      (build + run the 7)
#   make clean           清理 bin/
#
# 编译旗标固定 -Wall -Wextra -O2 -D_GNU_SOURCE，力求零告警。
# 七篇各一支验证器：一支 .c / 一支 .exe，一一对应。

ifeq ($(origin CC),default)
CC      := gcc
endif
CFLAGS  ?= -Wall -Wextra -O2 -D_GNU_SOURCE
LDLIBS  ?= -lm
BIN     := bin

# —— 七个猜想，七支验证器（顺序 = README 表序）——
VERIFS  := \
	bsd-conjecture/code/verify_bsd.c \
	hardy-littlewood/code/verify_hl.c \
	hodge-conjecture/code/verify_hodge.c \
	navier-stokes/code/verify_ns.c \
	p-vs-np/code/verify_pvsnp.c \
	riemann-hypothesis/code/verify_rh.c \
	yang-mills/code/verify_ym.c

BINS    := $(patsubst %.c,$(BIN)/%,$(VERIFS))


.PHONY: all run clean
all: $(BINS)

$(BIN)/%: %.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

run: all
	@mkdir -p audit/outs; fail=0; for b in $(BINS); do \
	  p=$$(basename $$b | sed 's/^verify_//'); \
	  printf '%-40s ' "$$b"; \
	  if out=$$(./$$b 2>&1); then \
	    echo "$$out" > audit/outs/$$p.txt; \
	    echo "$$(printf '%s\n' "$$out" | grep -oE 'RESULT PASS=[0-9]+ FAIL=[0-9]+' | tail -1) OK"; \
	  else echo "EXIT!=0"; fail=1; fi; \
	done; exit $$fail
clean:
	rm -rf $(BIN)

# ---- 审计器（不在七支之列；EXIT=0 表示被审对象无缺陷且标层完备）----
# ★ 必须先 `make run`（刷新 audit/outs/），审计器会活读那份实跑输出与台账/README 对账。
.PHONY: audit
audit:
	@mkdir -p $(BIN)
	$(CC) $(CFLAGS) -o $(BIN)/verify_audit audit/verify_audit.c $(LDLIBS)
	@$(BIN)/verify_audit matrix/audit_layer_gap.csv
