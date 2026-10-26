/* ── 落档合规（chk_judge 五條：判据可失败性）──
 *   五行六步：五行土S 锚定 -> 五行火F 虚化 -> 五行金M 坍缩 -> 五行水W 渗透 -> 五行木D 扩展 -> 五行土S 闭环
 *   四判據數（★ 判項計數口徑）：res=有R(F:文件↦违规行·R:违规行↦源行号·dev即行号偏差) · domain=0(判項·**實跑** 0P/0F) · dc=2(合规/违规) · order=未測（本門禁【無標度對象】）· **實跑狀態=NOOUT**
 *   F: 逐文件扫描 → 填六条判据的违规计数
 *   R: 反向測試（见 selftest-j1/j6：故意违规件必须被抓）
 */
/* chk_judge.c — 判据可失败性门禁（纯 C · 整文件版）
 *
 * 立此门禁的理由（2026-10-27 实测）：本仓曾有两支验证器「只打印、不可失败」——
 *   twin-prime/code/verify_twin_count.c 原 main 只有 build_sieve + 打印 + return 0；
 *   p-vs-np/code/verify_pvsnp.c 的 GF(2) 秩表、SAT 计数只打印，EXIT 只由两个浮点 dev 把门。
 * 这两支在【全 PASS】的外观下，其实【没有任何判据能失败】⟹ PASS 的判别力 = 0。
 * 本门禁把「可失败性」变成机械判据，六条：
 *
 *  J1 PASS 汇总 ⟹ 必须有 FAIL 汇总（只报通过不报失败的汇总行 = 宣告式）
 *  J2 PASS 汇总 ⟹ 必须有失败返回（return 1 / return bad / return fail / exit(1) / ? 1 : 0 之一）
 *  J3 打印 EXIT= ⟹ 必须有 return（EXIT 值必须由返回值给出，不能只打印）
 *  J4 出现 [PASS] 与 [FAIL] 必须成对（单侧出现 = 判据只有一半）
 *  J5 PASS=/FAIL= 后紧跟数字 ⟹ 违规（字面常量入格式串 = 宣告式汇总之源头）
 *  J6 判据阈值比较（`<= 1e-`）左侧不得是字面常量（0 / 1 / 1.0 等恒真或恒假）
 *
 * 用法：chk_judge.exe 文件…   退出码：0 全合规 / 1 有违规
 * 注意：本文件自身也必须合规 ⟹ 循环一律放在辅助函数内，main 直接层无 for/while。
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int g_viol = 0, g_files = 0;

static char *slurp(const char *path)
{
    FILE *f = fopen(path, "rb");
    char *b; long sz, r;
    if (!f) return NULL;
    fseek(f, 0, SEEK_END); sz = ftell(f); fseek(f, 0, SEEK_SET);
    b = (char *)malloc((size_t)sz + 1);
    if (!b) { fclose(f); return NULL; }
    r = (long)fread(b, 1, (size_t)sz, f);
    b[r] = 0; fclose(f); return b;
}

/* 注释抹成空格（保留字符串字面量：J1/J5 要看得到格式串） */
static void strip_comments(const char *s, char *d)
{
    long i = 0, j = 0;
    while (s[i]) {
        if (s[i] == '/' && s[i + 1] == '/') { while (s[i] && s[i] != '\n') { d[j++] = ' '; i++; } }
        else if (s[i] == '/' && s[i + 1] == '*') {
            d[j++] = ' '; d[j++] = ' '; i += 2;
            while (s[i] && !(s[i] == '*' && s[i + 1] == '/')) { d[j++] = (s[i] == '\n') ? '\n' : ' '; i++; }
            if (s[i]) { d[j++] = ' '; d[j++] = ' '; i += 2; }
        }
        else if (s[i] == '"' || s[i] == '\'') {
            char q = s[i]; d[j++] = s[i++];
            while (s[i] && s[i] != q) { if (s[i] == '\\' && s[i + 1]) { d[j++] = s[i++]; } d[j++] = s[i++]; }
            if (s[i]) d[j++] = s[i++];
        }
        else d[j++] = s[i++];
    }
    d[j] = 0;
}

static int has(const char *s, const char *w) { return strstr(s, w) != NULL; }

/* 数出现次数（用于 J4 成对性） */
static int count_of(const char *s, const char *w)
{
    int n = 0; const char *p = s;
    while ((p = strstr(p, w)) != NULL) { n++; p += strlen(w); }
    return n;
}

/* 字面量判定：字符是否是「十进制数字」 */
static int isdig(char c) { return c >= '0' && c <= '9'; }

/* J5：PASS= / FAIL= 后紧跟数字 ⟹ 宣告式。返回命中次数。 */
static int scan_declared(const char *s)
{
    int n = 0; const char *p = s;
    static const char *keys[2] = { "PASS=", "FAIL=" };
    int ki;
    for (ki = 0; ki < 2; ki++) {
        p = s;
        while ((p = strstr(p, keys[ki])) != NULL) {
            const char *q = p + strlen(keys[ki]);
            if (isdig(*q)) n++;
            p = q;
        }
    }
    return n;
}

/* J6：`<= 1e-` 左侧是字面常量。返回命中次数。 */
static int scan_constant_threshold(const char *s)
{
    int n = 0; const char *p = s;
    while ((p = strstr(p, "<= 1e-")) != NULL) {
        const char *q = p; int seen = 0, bad = 0;
        while (q > s && seen < 24) {
            q--; seen++;
            if (*q == '(' || *q == ';' || *q == '{' || *q == '}' || *q == '\n') break;
            if ((*q >= '0' && *q <= '9') || *q == '.') bad = 1;
            if (*q == ')' || *q == ']') { bad = 0; }   /* 括号内的表达式算变量侧 */
        }
        if (bad) n++;
        p += 6;
    }
    return n;
}

/* 失败返回的五种常见写法（必须能非零退出，而不是恒 return 0）：
 *   return 1;            return bad;  return nfail;  return g_fail;  return g_bad;
 *   return x == 0 ? 0 : 1;               → 命中 "? 0 : 1" / "?0:1"
 *   return (nfail == 0) ? 0 : 1;         → 同上
 *   exit(1);
 */
static int fail_return_present(const char *s)
{
    static const char *pats[] = {
        "? 0 : 1", "?0:1", "? 1 : 0", "?1:0",
        "return 1", "return bad", "return fail", "return (fail", "return (bad",
        "return nfail", "return g_fail", "return g_bad", "return (nfail", "return (g_fail",
        "return (g_bad", "exit(1)", "exit(EXIT_FAILURE)"
    };
    int i, n = (int)(sizeof pats / sizeof pats[0]);
    for (i = 0; i < n; i++) if (has(s, pats[i])) return 1;
    return 0;
}

static void judge_one(const char *path, const char *raw, const char *body)
{
    int has_summary = has(body, "PASS=");
    int j1 = 0, j2 = 0, j3 = 0, j4 = 0, j5 = 0, j6 = 0;

    if (has_summary && !has(body, "FAIL=")) j1 = 1;
    if (has_summary && !fail_return_present(body)) j2 = 1;
    if (has(body, "EXIT=") && !has(body, "return")) j3 = 1;
    if (count_of(body, "[PASS]") != count_of(body, "[FAIL]")) j4 = 1;
    if (scan_declared(body) > 0) j5 = 1;
    if (scan_constant_threshold(body) > 0) j6 = 1;

    if (j1 || j2 || j3 || j4 || j5 || j6) {
        printf("VIOLATION %s  J1=%d J2=%d J3=%d J4=%d J5=%d J6=%d\n", path, j1, j2, j3, j4, j5, j6);
        g_viol++;
    }
    (void)raw;
}

static int scan_all(int argc, char **argv)
{
    int i;
    for (i = 1; i < argc; i++) {
        char *raw = slurp(argv[i]);
        char *body;
        if (!raw) { printf("SKIP %s (unreadable)\n", argv[i]); continue; }
        body = (char *)malloc(strlen(raw) + 1);
        if (!body) { free(raw); continue; }
        strip_comments(raw, body);
        judge_one(argv[i], raw, body);
        g_files++;
        free(body); free(raw);
    }
    return g_viol;
}

int main(int argc, char **argv)
{
    if (argc < 2) { printf("usage: chk_judge.exe file.c ...\n"); return 2; }
    return scan_all(argc, argv) ? 1 : 0;
}
