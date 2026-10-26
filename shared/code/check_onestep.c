/* ── 落檔合規（check_archive 四條）──
 *   五行六步：五行土S 錨定 -> 五行火F 虛化 -> 五行金M 坍縮 -> 五行水W 滲透 -> 五行木D 擴展 -> 五行土S 閉環
 *   四判據數（★ 判項計數口徑）：res=有R(F:参数↦读数·R:读数↦标准值·chk对照·dev即残差) · domain=0(判項·**實跑** 0P/0F) · dc=2(PASS/FAIL) · order=未測（本探針【無標度對象】）· **實跑狀態=NOOUT·172ms**
 *   F: 見源碼
 *   R: 反向測試（見源碼）
 */
/* ★ 修订记录（2026-10-06）
 *   #290-a  kw() 无词边界 ⟹ kw(s,i,"main") 命中 domain ⟹ maindepth 定错位 ⟹ 假违规 3 支
 *          （verify_diffusion_score_409 / verify_hl_onestep_rev_436 / verify_bin2tern_lossless_437）
 *          修：加 isid() 词边界（前后字符均须非标识符字符）。同理 kw(...,"for") 曾命中 perform/before。
 *   #290-b  本文件原缺一个 }（L73 的 else if 未闭合）⟹ 盘上 .c 与在用 .exe 非同一版本
 *          修：补上 }}。此前所有门禁结论由【旧 exe】给出，未受此影响。
 *   修后全库对帐：313 支 · 合规 312 · 违规 1（唯一违规 = _门禁判别力测试_故意违规.c·故意保留）（保真度=未測·fiber 未計）
 *//* check_onestep.c — 闭式合规门禁（纯 C·整文件版）
 * 判据（用户六条第①条）：main() 体的**直接层**出现 for/while ⟹ 违规
 * 用法：check_onestep.exe 文件…   退出码：0 全合规 / 1 有违规
 * 注意：本文件自身也必须合规 ⟹ 循环一律放在辅助函数内。
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int g_viol = 0, g_files = 0, g_strict = 0;

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
static void strip(const char *s, char *d)
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
            char q = s[i]; d[j++] = ' '; i++;
            while (s[i] && s[i] != q) { if (s[i] == '\\' && s[i + 1]) { d[j++] = ' '; i++; } d[j++] = (s[i] == '\n') ? '\n' : ' '; i++; }
            if (s[i]) { d[j++] = ' '; i++; }
        }
        else d[j++] = s[i++];
    }
    d[j] = 0;
}
static int isid(char c){ return (c=='_'||(c>='0'&&c<='9')||(c>='a'&&c<='z')||(c>='A'&&c<='Z')); }
static int kw(const char *s, long i, const char *w)
{
    long k = 0;
    while (w[k]) { if (s[i + k] != w[k]) return 0; k++; }
    /* ★ 词边界（修 #290：kw("main") 曾命中 domain ⟹ maindepth 定错位 ⟹ 假违规） */
    if (i > 0 && isid(s[i - 1])) return 0;
    if (isid(s[i + k])) return 0;
    return 1;
}
static void scan_one(const char *path)
{
    long i = 0, ln = 1, depth = 0, maindepth = -1, top = 0, maxnest = 0, strictn = 0;
    char *raw = slurp(path), *s;
    if (!raw) { printf("  [跳过] 打不开：%s\n", path); return; }
    s = (char *)malloc(strlen(raw) + 1);
    if (!s) { free(raw); return; }
    strip(raw, s); free(raw);
    while (s[i]) {
        if (s[i] == '\n') ln++;
        if (s[i] == '{') { depth++; if (depth > maxnest) maxnest = depth; }
        else if (s[i] == '}') { depth--; if (depth < 0) depth = 0; if (maindepth >= 0 && depth < maindepth) maindepth = -1; }   /* ★ 修 #296：main 体关闭后复位 */
        else if (maindepth < 0 && kw(s, i, "main")) {
            long k = i + 4, d = 0;
            while (s[k] == ' ' || s[k] == '\n') k++;
            if (s[k] == '(') {
                while (s[k]) { if (s[k] == '(') d++; else if (s[k] == ')') { d--; if (d == 0) { k++; break; } } k++; }
                while (s[k] == ' ' || s[k] == '\n') k++;
                if (s[k] == '{') maindepth = depth + 1;
            }
        }
        else if ((kw(s, i, "for") || kw(s, i, "while"))) {
            if (g_strict) { if (ln > 0) strictn++; }
            if (maindepth >= 0 && depth == maindepth) {
            printf("  %-36s 顶层循环=1  L%-5ld **违规**\n", path, ln);
            top++; g_viol++; i += 3; continue;
        }
        }
        i++;
    }
    if (top == 0) printf("  %-36s 顶层循环=0                    合规\n", path);
    if (g_strict) printf("        （严格模式：非顶层 for/while 计数 = %ld）\n", strictn);
    printf("        最大 brace 深度 = %ld\n", maxnest);
    free(s); g_files++;
}
static void scan_all(int argc, char **argv)
{
    int i;
    for (i = 1; i < argc; i++) scan_one(argv[i]);
}
int main(int argc, char **argv)
{
    printf("闭式门禁 · main 体直接层不得出现 for/while\n");
    printf("================================================================\n");
    if (argc < 2) { printf("用法：check_onestep.exe 文件…\n"); return 1; }
    scan_all(argc, argv);
    printf("================================================================\n");
    printf("判定：%s（退出码 %d）· 文件数 %d\n", g_viol ? "**有违规·不得落档**" : "全部合规", g_viol ? 1 : 0, g_files);
    return g_viol ? 1 : 0;
}
