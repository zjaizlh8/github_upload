# 一键复算 / One-command recompute
#
#   make                 编译 7 个验证器到 bin/         (7 verifiers)
#   make run             编译 + 逐支实跑 7 个验证器      (build + run the 7)
#   make tools           （可选）编 5 支归档检查工具      (optional: archive tools)
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

# —— 可选归档工具（不属于七篇，不参与 make run）——
TOOLS_SRC := $(wildcard shared/tools/*.c)
TOOLS_BIN := $(patsubst shared/tools/%.c,$(BIN)/tools/%,$(TOOLS_SRC))

.PHONY: all run tools clean
all: $(BINS)

$(BIN)/%: %.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

run: all
	@fail=0; for b in $(BINS); do \
	  printf '%-40s ' "$$b"; \
	  if out=$$(./$$b 2>&1); then \
	    echo "$$(printf '%s\n' "$$out" | grep -oE 'RESULT PASS=[0-9]+ FAIL=[0-9]+' | tail -1) OK"; \
	  else echo "EXIT!=0"; fail=1; fi; \
	done; exit $$fail

tools: $(TOOLS_BIN)

$(BIN)/tools/%: shared/tools/%.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

clean:
	rm -rf $(BIN)
