/* ── 落檔合規（check_archive 四條）──
 *   五行六步：五行土S 錨定 -> 五行火F 虛化 -> 五行金M 坍縮 -> 五行水W 滲透 -> 五行木D 擴展 -> 五行土S 閉環
 *   四判據數（★ 判項計數口徑）：res=有R(F:参数↦读数·R:读数↦标准值·chk对照·dev即残差) · domain=0(判項·**實跑** 0P/0F) · dc=2(PASS/FAIL) · order=未測（本探針【無標度對象】）· **實跑狀態=NOOUT·312ms**
 *   F: 名
 *   R: 無反向測試（單向·非閉環·照實標）
 */
/* chk_scan.c —— 判據空洞普查（純 C）
 *
 * 目的：找出 .c 探針裏【恆真判據】。按 §6.27 第 3 條：
 *       3 參簽名 chk(1, "名", "訊息") ⟹ 恆真（判別力 0）
 *       5 參簽名 chk(ok, "名", got, exp, tol) ⟹ 真判據（|got−exp| ≤ tol）
 *   ⟹ 本工具只認【第一個實參恰為字面量 1（或 0/true）】且【實參總數 ≤ 4】的呼叫。
 *
 * 用法：  chk_scan.exe <file1.c> [file2.c ...]
 *         chk_scan.exe *.c               （由 shell 展開）
 * 輸出：  每檔一行摘要 + 末尾合計；並寫 _chk_scan.txt（UTF-8 含 BOM）
 *
 * 編譯：  gcc -O2 -o chk_scan.exe chk_scan.c
 *
 * ★ 修 v2（2026-10-26·闭式不暴露）：main 直接層零循環——argv 遍歷搬進 run_all()，
 *   輸出文件段也收進 run_all()，main 只留入口。
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXFILES 4000

typedef struct {
    char name[260];
    int  total;      /* chk( 出現次數 */
    int  taut;       /* 恆真：首參為 1/0/true 且實參數 ≤ 4 */
    int  taut3;      /* 恰好 3 參且首參為 1 */
} FileStat;

/* 從 p（指向 '(' 之後）起，數【頂層逗號】+1 = 實參個數；並判斷首參是不是裸的 1/0/true
 * 回傳 1 = 恆真型，0 = 非；args 回傳實參數（若指標非 NULL） */
static int classify_args(const char *p, int *args){
    int depth = 1, n = 1, i = 0, is_taut = 0;
    char first[32];
    const char *q = p;
    /* 跳過空白取首參 */
    while (*q == ' ' || *q == '\t' || *q == '\n' || *q == '\r') q++;
    while (*q && *q != ',' && *q != ')' && i < 31) first[i++] = *q++;
    first[i] = 0;
    while (first[i-1] == ' ' || first[i-1] == '\t') first[--i] = 0;
    if (!strcmp(first,"1") || !strcmp(first,"0") ||
        !strcmp(first,"true") || !strcmp(first,"false") ||
        !strcmp(first,"NULL")) is_taut = 1;

    /* 數逗號到配對的 ')' */
    while (*p){
        if (*p == '(') depth++;
        else if (*p == ')'){ depth--; if (!depth) break; }
        else if (*p == ',' && depth == 1) n++;
        else if (*p == '"'){ /* 跳過字串 */ p++; while (*p && *p != '"'){ if (*p=='\\') p++; p++; } }
        p++;
    }
    if (args) *args = n;
    return is_taut;
}

/* ★★ 修（2026-10-06·自曝 #382）：必須【跳過註釋】。
 *   原版把註釋裏出現的 `chk(1,…)`（例如「原為 chk(1,…) 恆真」這句說明）也數進去
 *   ⟹ 空洞率被【高估】。
 * ★★ 再修（2026-10-06·自曝 #392）：必須【跳過字符串字面量的內容】。
 *   實例：bidir15_627 的判據名稱寫着「（…·原為 chk(1)）」⟹ 掃描器把它當成一次 chk 調用
 *   ⟹ 報 6 處假陽性。做法：把字符串【內容】就地換成空白（保留引號本身，以免破壞解析）。 */
static void strip_comments(char *b, long n){
    long i = 0;
    while (i < n){
        /* ── 字符串字面量：保留引號，內容換空白 ── */
        if (b[i] == '"'){
            i++;
            while (i < n && b[i] != '"'){
                if (b[i] == '\\'){ if (b[i] != '\n') b[i] = ' '; i++;
                                   if (i < n){ if (b[i] != '\n') b[i] = ' '; i++; }
                                   continue; }
                if (b[i] != '\n') b[i] = ' ';
                i++;
            }
            if (i < n) i++;            /* 跳過收尾引號 */
            continue;
        }
        /* ── 字符常量：同理 ── */
        if (b[i] == '\''){
            i++;
            while (i < n && b[i] != '\''){
                if (b[i] == '\\'){ if (b[i] != '\n') b[i] = ' '; i++;
                                   if (i < n){ if (b[i] != '\n') b[i] = ' '; i++; }
                                   continue; }
                if (b[i] != '\n') b[i] = ' ';
                i++;
            }
            if (i < n) i++;
            continue;
        }
        /* ── 塊註釋 ── */
        if (b[i] == '/' && i + 1 < n && b[i+1] == '*'){
            b[i] = ' '; b[i+1] = ' '; i += 2;
            while (i < n && !(b[i] == '*' && i + 1 < n && b[i+1] == '/')){
                if (b[i] != '\n') b[i] = ' ';
                i++;
            }
            if (i < n){ b[i] = ' '; if (i+1 < n) b[i+1] = ' '; i += 2; }
            continue;
        }
        /* ── 行註釋 ── */
        if (b[i] == '/' && i + 1 < n && b[i+1] == '/'){
            while (i < n && b[i] != '\n'){ b[i] = ' '; i++; }
            continue;
        }
        i++;
    }
}

static int scan_file(const char *path, FileStat *st){
    FILE *f = fopen(path, "rb");
    long sz; char *buf; long i;
    if (!f) return 0;
    fseek(f, 0, SEEK_END); sz = ftell(f); fseek(f, 0, SEEK_SET);
    if (sz <= 0){ fclose(f); return 0; }
    buf = (char *)malloc(sz + 1);
    if (!buf){ fclose(f); return 0; }
    if (fread(buf, 1, sz, f) != (size_t)sz){ free(buf); fclose(f); return 0; }
    buf[sz] = 0; fclose(f);
    strip_comments(buf, sz);   /* ★ 修 #382：先去掉註釋，再掃 */

    strncpy(st->name, path, sizeof st->name - 1);
    st->name[sizeof st->name - 1] = 0;

    for (i = 0; i + 4 < sz; i++){
        if ((i == 0 || !((buf[i-1] >= 'a' && buf[i-1] <= 'z') ||
                         (buf[i-1] >= 'A' && buf[i-1] <= 'Z') ||
                         (buf[i-1] >= '0' && buf[i-1] <= '9') || buf[i-1] == '_')) &&
            buf[i] == 'c' && buf[i+1] == 'h' && buf[i+2] == 'k' && buf[i+3] == '('){
            int args = 0, t;
            st->total++;
            t = classify_args(buf + i + 4, &args);
            if (t && args <= 4){ st->taut++; if (args == 3) st->taut3++; }
        }
    }
    free(buf);
    return 1;
}

static int run_all(int argc, char **argv){
    FileStat st[MAXFILES];
    int nf = 0, i, tot = 0, taut = 0, taut3 = 0, dirty = 0;
    FILE *out;

    for (i = 1; i < argc && nf < MAXFILES; i++){
        memset(&st[nf], 0, sizeof st[nf]);
        if (scan_file(argv[i], &st[nf])){
            tot += st[nf].total; taut += st[nf].taut; taut3 += st[nf].taut3;
            if (st[nf].taut > 0){
                printf("  %-44s chk=%4d  taut=%4d (3-arg %4d)  %5.1f%%\n",
                       st[nf].name, st[nf].total, st[nf].taut, st[nf].taut3,
                       st[nf].total ? 100.0 * st[nf].taut / st[nf].total : 0.0);
                dirty++;
            }
            nf++;
        } else {
            printf("  %s  (cannot read)\n", argv[i]);
        }
    }
    printf("\n  ==== files=%d  with-tautology=%d  chk-total=%d  taut-total=%d (%.1f%%)  3-arg=%d ====\n",
           nf, dirty, tot, taut, tot ? 100.0 * taut / tot : 0.0, taut3);

    out = fopen("_chk_scan.txt", "wb");
    if (out){
        fprintf(out, "\xEF\xBB\xBF");
        fprintf(out, "chk_scan: files=%d with-tautology=%d chk-total=%d taut-total=%d (%.1f%%) 3-arg=%d\n",
                nf, dirty, tot, taut, tot ? 100.0 * taut / tot : 0.0, taut3);
        for (i = 0; i < nf; i++)
            if (st[i].taut > 0)
                fprintf(out, "  %-44s chk=%4d taut=%4d 3-arg=%4d\n",
                        st[i].name, st[i].total, st[i].taut, st[i].taut3);
        fclose(out);
    }
    return 0;
}

int main(int argc, char **argv){
    if (argc < 2){
        printf("usage: chk_scan.exe <file1.c> [file2.c ...]\n");
        return 2;
    }
    return run_all(argc, argv);
}
