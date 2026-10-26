# 一键复算：make          → 编译全部 verify_*.c / chk_*.c 到 bin/
#           make run      → 跑全部验证器
#           make gates    → 用 4 道落档门禁扫全部 .c
#           make all-verify → 编译 + 跑 + 门禁
# 编译旗标固定 -Wall -Wextra -O2（与仓库口径一致），力求零告警。

# make 内置默认 CC=cc 会让 ?= 失效（Windows 上无 cc），故只在"默认来源"时改成 gcc
ifeq ($(origin CC),default)
CC      := gcc
endif
CFLAGS  ?= -Wall -Wextra -O2 -D_GNU_SOURCE
LDLIBS  ?= -lm
BIN     := bin
SRC     := $(shell find . -name '*.c' | sort)
BINS    := $(patsubst ./%,$(BIN)/%,$(basename $(SRC)))
GATES   := check_onestep chk_loops chk_res chk_scan
# 验证器 = 全部 .c 去掉 4 道门禁源码 = 20 支（含 shared/code/verify_simlayer.c 共用模拟器）
VERIFS  := $(filter-out $(addprefix $(BIN)/shared/code/,$(GATES)),$(BINS))

.PHONY: all all-verify run gates clean
all: $(BINS)

$(BIN)/%: ./%.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

run: all
	@fail=0; for b in $(VERIFS); do \
	  printf '%-34s ' "$$b"; \
	  if out=$$(timeout 300 ./$$b 2>&1); then \
	    echo "$$(printf '%s\n' "$$out" | grep -oE 'PASS=[0-9]+ FAIL=[0-9]+' | head -1) OK"; \
	  else echo "EXIT!=0"; fail=1; fi; \
	done; exit $$fail

gates: all
	@fail=0; for g in $(GATES); do \
	  printf '%-16s ' "$$g"; \
	  if ./$(BIN)/shared/code/$$g $(SRC) >/dev/null 2>&1; then echo "PASS (exit 0)"; else echo "FAIL"; fail=1; fi; \
	done; exit $$fail

all-verify: run gates

clean:
	rm -rf $(BIN)
