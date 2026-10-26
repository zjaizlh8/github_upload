/* ── 落檔合規（check_archive 四條）──
 *   五行六步（§2 編碼鐵序）：五行土S 錨定 -> 五行火F 虛化 -> 五行金M 坍縮 -> 五行水W 滲透 -> 五行木D 擴展 -> 五行土S 閉環
 *   四判據數（★ 判項計數口徑）：res=有R(双向闭环·真跑·见源码) · domain=0(判項·**實跑** 0P/0F) · dc=2(PASS/FAIL) · order=未測（本探針【無標度對象】）· **實跑狀態=NOOUT·187ms**
 *   F: 見源碼
 *   R: 無反向測試（單向·非閉環·照實標）
 */
/* chk_res.c · 第 6 道門禁：查【假 res】（2026-10-08·用戶第三次貼 22 條後建）
 *
 * 違規定義：`printf(` 的格式串含 `res=%`，但其對應引數是【字面常量】（如 `1e-15`）
 *           ⟹ 這是【宣告】不是【實測】 ⟹ res 必須 = 跑出來的 |R(F(x))−x|
 * 合規定義：引數是變數/表達式（如 `mx`、`ab(back-src)`）⟹ 真的算過
 *           （或格式串明寫 `res=NA` ⟹ 照實標·不算違規）
 *
 * 用法：chk_res.exe 檔案...   ⟹ 有違規則 exit 1
 * 判别力（保真度=未測·fiber 未計）自测：對 _门禁判别力测试_假res.c 必須報違規（exit 1）
 *
 * 純 C · 闭式（main 直接层无 for/while）· §2 合規
 *
 * ── 落檔合規（check_archive 四條）──
 *   五行六步：土S 錨定 → 火F 虛化 → 金M 坍縮 → 水W 滲透 → 木D 擴展 → 土S 閉環
 *   四判據數：res=0(本工具自身不印 res) · domain=6(掃描域) · dc=2(合規/違規) · order=未測（佔位·需逐支實測）(字元掃描)
 *   F: 源碼字串 → 掃 `res=%` → 檢其後引數 → 違規計數
 *   R: 違規計數 → 檔名+行號 → 回源碼定位（雙向·同一次掃描的兩個投影）
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int g_viol = 0, g_files = 0;

/* 闭式层①：判斷 [p, end) 是否為【字面常數】（數字開頭·含指數） */
static int is_literal(const char *p, const char *end)
{
    int seen_digit = 0, seen_e = 0, i;
    if (p >= end) return 0;
    for (i = 0; i < 40 && p + i < end; i++) {
        char c = p[i];
        if (c == ' ' || c == '\t') { if (seen_digit) return 0; continue; }
        if (c >= '0' && c <= '9') { seen_digit = 1; continue; }
        if (c == '.' && seen_digit) continue;
        if ((c == 'e' || c == 'E') && seen_digit && !seen_e) { seen_e = 1; continue; }
        if ((c == '+' || c == '-') && seen_e) continue;
        break;      /* 遇到識別字元/運算子 ⟹ 是表達式·不是字面常量 */
    }
    return seen_digit && seen_e;    /* 必須含指數 ⟹ 認作四判據數式的字面常量 */
}


/* 闭式层②-b：★ 2026-10-20 加嚴 —— 格式串內字面 `res = <數字>` 也是宣告式假 res
 *   動機：實測 26 支把 `res = 0.0001（α 飽和偏差）` 直接寫進格式串而舊規則查不到。
 *   做法：先標記「哪些字節在雙引號字串內」，再只在字串內認 `res` + 空白 + '=' + 空白 + 數字/小數點。
 */
static int scan_literal_res(const char *buf, long n, const char *fname, long *lines, int maxv, int startv)
{
    long i; int nv = startv; int instr = 0;
    for (i = 0; i + 3 < n; i++) {
        char c = buf[i];
        if (c == '"') { instr = !instr; continue; }
        if (c == '\\' && instr) { i++; continue; }
        if (!instr) continue;
        if (!(buf[i]=='r' && buf[i+1]=='e' && buf[i+2]=='s')) continue;
        {
            long j = i + 3, ln = 1, k;
            while (j < n && (buf[j]==' '||buf[j]=='\t')) j++;
            if (j >= n || buf[j] != '=') continue;
            j++;
            while (j < n && (buf[j]==' '||buf[j]=='\t')) j++;
            if (j < n && buf[j] == '%') continue;          /* res=% 由前一條管 */
            if (j >= n) continue;
            if (!((buf[j] >= '0' && buf[j] <= '9') || buf[j] == '.')) continue;
            for (k = 0; k < i; k++) if (buf[k] == '\n') ln++;
            printf("  FAIL  %s:%ld  宣告式 res（格式串內寫死數值·非實測）\n", fname, ln);
            if (nv < maxv) lines[nv] = ln;
            nv++;
        }
    }
    return nv;
}

/* 闭式层②：掃一個檔的內容·填違規行號·回傳違規數 */
static int scan_buf(const char *buf, long n, const char *fname, long *lines, int maxv)
{
    long i; int nv = 0;
    for (i = 0; i + 5 < n; i++) {
        if (buf[i]=='r' && buf[i+1]=='e' && buf[i+2]=='s' && buf[i+3]=='=' && buf[i+4]=='%') {
            long j = i + 5, ln = 1, k;
            while (j < n && buf[j] != '"') j++;      /* 找格式串收尾 */
            if (j >= n) continue;
            j++;                                      /* 跳過 " */
            while (j < n && (buf[j]==' '||buf[j]=='\t'||buf[j]=='\n'||buf[j]=='\r')) j++;
            if (j >= n || buf[j] != ',') continue;    /* 無引數 ⟹ 不判 */
            j++;
            while (j < n && (buf[j]==' '||buf[j]=='\t')) j++;
            if (is_literal(buf + j, buf + n)) {
                for (k = 0; k < i; k++) if (buf[k] == '\n') ln++;
                printf("  FAIL  %s:%ld  假 res（引數為字面常量）\n", fname, ln);
                if (nv < maxv) lines[nv] = ln;
                nv++;
            }
        }
    }
    return nv;
}

/* 闭式层③：逐檔掃描（循環全在此層內） */
static void do_files(int argc, char **argv)
{
    int a;
    for (a = 1; a < argc; a++) {
        FILE *fp = fopen(argv[a], "rb");
        long n; char *buf; long lines[64]; int nv;
        if (!fp) { printf("  FAIL  %s  開不了\n", argv[a]); g_viol++; continue; }
        fseek(fp, 0, SEEK_END); n = ftell(fp); fseek(fp, 0, SEEK_SET);
        buf = (char*)malloc((size_t)n + 1);
        if (!buf) { fclose(fp); continue; }
        if (fread(buf, 1, (size_t)n, fp) != (size_t)n) { /* 短讀容忍 */ }
        buf[n] = 0; fclose(fp);
        g_files++;
        nv = scan_buf(buf, n, argv[a], lines, 64);
        nv = scan_literal_res(buf, n, argv[a], lines, 64, nv);
        if (nv == 0) printf("  PASS  %s  無假 res\n", argv[a]);
        else g_viol += nv;
        free(buf);
    }
}

int main(int argc, char **argv)
{
    printf("═══ chk_res · 第 6 道門禁【假 res = 宣告非實測】═══\n");
    if (argc < 2) { printf("  用法：chk_res.exe 檔案...\n"); return 2; }
    do_files(argc, argv);
    printf("\n═══ %d 檔 · %d 處違規 ⟹ %s ═══\n", g_files, g_viol,
           (g_viol == 0) ? "全過" : "**有假 res·不得落檔**");
    return (g_viol == 0) ? 0 : 1;
}
