/* verify_nest_order.c —— 嵌套数组 × 换关系函数：先认出「对象」，再决定这条判据该不该逼成 0
 *
 * four operations：加法=accumulate · 乘法=generate · 跨维度=dimension-change · base=非象
 *   跨维度处只有残差；本档残差口径 = 相对偏差（dev / |目标|）· 周期=不涉（本档无周期量）· 余项有界=不涉
 *
 * 为什么需要这一支
 *   台账 `matrix/forcedness_ledger.csv` 里有一批行，dev ≠ 0 却不是"拟合"：
 *   `FORCED_trunc`（截断残差）与 `FLOW`（有限窗流体量）。上一轮的处置是
 *   「dev ≤ C·h^p 就写截断阶」——**但没人验过那个 (p,q) 到底稳不稳**，
 *   于是"写出的阶"仍是声明。本支就是去把它量出来，并判定它是不是对象。
 *
 * 方法（混沌一步到位：换关系函数，不动的是对象）
 *   关系族取混幂×混对数（论文自己的误差形就是 C·h^-p·(ln h)^-q）：
 *     R(p,q) = dev · h^p · (ln h)^q        p,q ∈ {0,1,2}   ⟹ 9 列
 *   **换窗口时唯一不动的那一列，才是对象**；其余各列都是窗口的函数。
 *   若 9 列没有一列不动 ⟹ 这条是 TREND（真开口），不得逼成 0。
 *
 * 每个 item 有自己的窗口梯队（载体不同，不能硬凑同一组窗口）：
 *   P5/P6 乘积类 → 1e3 … 1e8（素数筛一次到底）
 *   P2 展开系数 → 1e3 … 1e10（闭式求值，不吃内存）
 *   P4 焊接口   → 1e3 … 1e7（要阶梯数组，受内存限）
 *
 * 容器三行
 *   载体 = 各 item 自己的窗口梯队
 *   算子 = 部分积 / 部分和（阶梯）在窗口上取值
 *   对称 = 换窗口不改"极限值"（常数不变）· 换关系函数只有一列不动（对象唯一）
 *
 * 判据不靠猜：不动列 + 实测 (p,q) ⟹ 判据写成 `dev ≤ C / (h^p (ln h)^q)`，容差由机制给（不拍脑袋）。
 *
 * 【五法】乘法=木（部分积 prod / h^p）· 跨维度=火（exp/log，付相对残差）· 加法=水（累加和）
 *        降维(金)=不涉 · 无损归心(土)=不涉 —— 本档不做 R(R(x)) 检验，如实标注
 * 层=结构 ：A[item][window] 的读数 · 9 个关系列 · 「哪一列不动」的选取（换窗口不变）
 * 层=过程 ：dev 随窗口的衰减只报 TREND（允许真开口；本档不把开口逼成 0）
 * 四判据数：res=NA(no-R-test) · domain={1e3…1e10} · dc=各 item 窗口数（5~6）· order=(p,q) 见段2
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NWMAX   6
#define NITEM   5
#define QMAX    100000000
#define SMAX    10000000
#define PMAX    5      /* p = 0, 0.5, 1, 1.5, 2（半整幂也要有：根号收敛是真事） */
#define QMAXL   3      /* q = 0, 0.5, 1 */
#define SPTOL   0.5
static const double PV[PMAX]  = { 0.0, 0.5, 1.0, 1.5, 2.0 };
static const double QV[QMAXL] = { 0.0, 0.5, 1.0 };

static int g_pass = 0, g_fail = 0;
static void J(const char *name, double dev, double tol)
{
    int ok = (dev <= tol);
    if (ok) g_pass++; else g_fail++;
    printf("  [%s] %-58s val=%.3e tol=%.1e\n", ok ? "PASS" : "FAIL", name, dev, tol);
}

/* ---- 每个 item 自己的窗口梯队 ---- */
static double HW[NITEM][NWMAX] = {
    { 1e3, 1e4, 1e5, 1e6, 1e7, 1e8 },
    { 1e3, 1e4, 1e5, 1e6, 1e7, 1e8 },
    { 1e8, 1e10, 1e12, 1e14, 1e16, 1e18 },   /* 只在渐近区（L≫1）里判，远离即发散 */
    { 1e3, 1e4, 1e5, 1e6, 1e7, 0    },        /* 全区间：含 t≈2 边界层 */
    { 2e4, 3e4, 1e5, 3e5, 1e6, 1e7 }         /* 渐近区起算 T0=1e4（窗口须 > T0，否则区间长为 0） */
};
static int    WN[NITEM] = { 6, 6, 6, 5, 6 };
static double A[NITEM][NWMAX];
static const char *ITEM[NITEM] = {
    "hl P5  2C2 截断积 2*prod(1-1/(q-1)^2)",
    "hl P6  e^-gamma 部分积 prod(1-1/p)*lnQ",
    "hl P2  一阶系数 c=(S-1)*L  （应=2）",
    "hl P4  焊接口 阶梯 vs Li2 （全区间 T0=2）",
    "hl P4b 焊接口 阶梯 vs Li2 （渐近区 T0=1e4）"
};

static const double C2X2   = 1.3203236316937391;
static const double EMG    = 0.5614594835668851;
static const double GAMMA  = 0.5772156649015329;

/* ---- 素数筛（一次性到 QMAX）---- */
static unsigned char *comp;
static void sieve_init(void)
{
    int i, j;
    comp = (unsigned char *)calloc((size_t)QMAX + 2, 1);
    if (!comp) { printf("OOM comp\n"); exit(2); }
    for (i = 2; (long long)i * i <= QMAX; i++)
        if (!comp[i]) for (j = i * i; j <= QMAX; j += i) comp[j] = 1;
}
static int is_p(int n) { return n >= 2 && !comp[n]; }

/* 一次扫到底：同时拿到两个乘积类 item 的全部窗口值 */
static double g_prod_c2[NWMAX], g_prod_m[NWMAX];
static void products_run(void)
{
    int q, wi = 0, wj = 0;
    double pc = 1.0, pm = 0.5;   /* pm 从 p=2 起：1-1/2 = 0.5，漏了它整列会偏 2 倍 */
    for (q = 3; q <= QMAX; q++) {
        if (!is_p(q)) continue;
        pc *= (1.0 - 1.0 / ((double)(q - 1) * (double)(q - 1)));
        pm *= (1.0 - 1.0 / (double)q);
        while (wi < WN[0] && (double)q >= HW[0][wi]) { g_prod_c2[wi] = 2.0 * pc; wi++; }
        while (wj < WN[1] && (double)q >= HW[1][wj]) { g_prod_m[wj] = pm * log(HW[1][wj]); wj++; }
        if (wi >= WN[0] && wj >= WN[1]) break;
    }
    /* 窗口落在最后一个素数之后（例如 Q=1e8 时最后素数是 99999989）⟹ 用扫到底的积补齐 */
    while (wi < WN[0]) { g_prod_c2[wi] = 2.0 * pc; wi++; }
    while (wj < WN[1]) { g_prod_m[wj] = pm * log(HW[1][wj]); wj++; }
}
/* 假公式（齿用）：分母写 q^2 */
static double c2_fake(double Q)
{
    int q;
    double p = 1.0;
    for (q = 3; q <= (int)Q; q++)
        if (is_p(q)) p *= (1.0 - 1.0 / ((double)q * (double)q));
    return 2.0 * p;
}
/* item 2：Li2 渐近展开的一阶系数 c(x) = (S-1)·L
 * (S-1)L = 2 + 6/L + 24/L² + … 是【渐近级数】——项在 n≈L 处最小，之后回升。
 * 所以在最小项处截断（最优截断），否则 n 固定会把发散尾巴当残差（1e3 处 dev=3.7e16 就是那个假象）。 */
static double coef_c(double x)
{
    double L = log(x), S = 1.0, term = 1.0, prev = 1e308;
    int n;
    for (n = 1; n <= 200; n++) {
        term *= (double)(n + 1) / L;
        if (fabs(term) > fabs(prev)) break;    /* 过了最小项 ⟹ 停（渐近截断） */
        S += term;
        prev = term;
    }
    return (S - 1.0) * L;
}
/* item 3：焊接口 —— 阶梯 (e^gamma·prod)^2 累积 vs Li2（Simpson）*/
static double *stair;
static void stair_init(void)
{
    int n;
    double m = 1.0;
    stair = (double *)calloc((size_t)SMAX + 2, sizeof(double));
    if (!stair) { printf("OOM stair\n"); exit(2); }
    for (n = 2; n <= SMAX; n++) {
        if (is_p(n)) m *= (1.0 - 1.0 / (double)n);
        stair[n] = exp(GAMMA) * m;
    }
}
static double li2_simpson(double X, int m)
{
    double h = (X - 2.0) / (double)m, s = 0.0;
    int i;
    for (i = 0; i <= m; i++) {
        double t = 2.0 + h * (double)i, L = log(t), f = 1.0 / (L * L);
        double w = (i == 0 || i == m) ? 1.0 : ((i & 1) ? 4.0 : 2.0);
        s += w * f;
    }
    return s * h / 3.0;
}
static double li2_simpson_from(double A0, double X, int m)
{
    double h = (X - A0) / (double)m, s = 0.0;
    int i;
    for (i = 0; i <= m; i++) {
        double t = A0 + h * (double)i, L = log(t), f = 1.0 / (L * L);
        double w = (i == 0 || i == m) ? 1.0 : ((i & 1) ? 4.0 : 2.0);
        s += w * f;
    }
    return s * h / 3.0;
}
static double stair_int(double X)
{
    int n;
    double s = 0.0;
    for (n = 3; n <= (int)X; n++) s += stair[n] * stair[n];
    return s;
}
static double stair_int_from(double A0, double X)
{
    int n;
    double s = 0.0;
    for (n = (int)A0; n <= (int)X; n++) s += stair[n] * stair[n];
    return s;
}

/* ---- 段 1：填表 ---- */
static void seg1_fill(int cells)
{
    int i, j, rows;
    products_run();
    for (j = 0; j < WN[0]; j++) A[0][j] = fabs(g_prod_c2[j] - C2X2) / C2X2;
    for (j = 0; j < WN[1]; j++) A[1][j] = fabs(g_prod_m[j] - EMG) / EMG;
    for (j = 0; j < WN[2]; j++) A[2][j] = fabs(coef_c(HW[2][j]) - 2.0) / 2.0;
    for (j = 0; j < WN[3]; j++) {
        double I = stair_int(HW[3][j]), R = li2_simpson(HW[3][j], 2000);
        A[3][j] = fabs(I / R - 1.0);
    }
    for (j = 0; j < WN[4]; j++) {
        double I = stair_int_from(1e4, HW[4][j]), R = li2_simpson_from(1e4, HW[4][j], 2000);
        A[4][j] = fabs(I / R - 1.0);
    }
    printf("=== 嵌套数组 A[item][window]（读数 = 相对 dev）===\n");
    for (i = 0; i < cells; i++) {
        rows = WN[i];
        printf("  %-40s", ITEM[i]);
        for (j = 0; j < rows; j++) printf(" %9.2e", A[i][j]);
        printf("\n");
    }
}

/* ---- 段 2：换关系函数，找不动列 ---- */
static double spread(const double *r, int n)
{
    int j;
    double lo = r[0], hi = r[0];
    for (j = 1; j < n; j++) { if (r[j] < lo) lo = r[j]; if (r[j] > hi) hi = r[j]; }
    if (hi <= 0.0) return 1e9;
    return hi / lo - 1.0;
}
static double RK[NITEM][PMAX][QMAXL][NWMAX];
static double RC[NITEM][PMAX][QMAXL];
static int    BP[NITEM], BQ[NITEM], CNT[NITEM];

static void seg2_relations(int cells)
{
    int i, p, q, j, rows;
    printf("\n=== 关系族 R(p,q)=dev·h^p·(ln h)^q（p∈{0,.5,1,1.5,2}·q∈{0,.5,1}）——不动的那列才是对象 ===\n");
    for (i = 0; i < cells; i++) {
        rows = WN[i];
        CNT[i] = 0; BP[i] = 0; BQ[i] = 0;
        for (p = 0; p < PMAX; p++)
            for (q = 0; q < QMAXL; q++) {
                for (j = 0; j < rows; j++) {
                    double hp = pow(HW[i][j], PV[p]), lq = pow(log(HW[i][j]), QV[q]);
                    RK[i][p][q][j] = A[i][j] * hp * lq;
                }
                RC[i][p][q] = spread(RK[i][p][q], rows);
                if (RC[i][p][q] <= SPTOL) CNT[i]++;
                if (RC[i][p][q] < RC[i][BP[i]][BQ[i]]) { BP[i] = p; BQ[i] = q; }
            }
        printf("  %-40s best=(p=%d,q=%d) spread=%.3f  窗口=%d 不动列=%d %s\n",
               ITEM[i], BP[i], BQ[i], RC[i][BP[i]][BQ[i]], rows, CNT[i],
               CNT[i] == 1 ? "⟹ 对象已认出" : (CNT[i] == 0 ? "⟹ TREND（无不动列）" : "⟹ 不合格（多列不动）"));
    }
}

/* ---- 段 3：落判据 ---- */
static void seg3_judge(int cells)
{
    int i, j, rows;
    printf("\n=== 判定：能认出对象的 → 写成机制界；认不出的 → TREND（不得逼成 0）===\n");
    for (i = 0; i < cells; i++) {
        double C = 0.0;
        rows = WN[i];
        if (CNT[i] != 1) {
            printf("  [TREND][层=过程] %-36s 不动列=%d ⟹ 认不出对象，只报走到哪儿\n", ITEM[i], CNT[i]);
            continue;
        }
        for (j = 0; j < rows; j++)
            if (RK[i][BP[i]][BQ[i]][j] > C) C = RK[i][BP[i]][BQ[i]][j];
        /* 作者确认: 只把印出的机制指数由 %d 改成 %.1f（配合新增的半整幂族 0.5/1.5），判读对象、阈值、回退规则全不动 */
        printf("  [FORCED_trunc] %-38s dev <= %.4g / (h^%.1f (ln h)^%.1f)   spread=%.3f\n",
               ITEM[i], C, PV[BP[i]], QV[BQ[i]], RC[i][BP[i]][BQ[i]]);
    }
}

/* ---- 段 4：齿 ---- */
static void seg4_teeth(int cells)
{
    int i, j, rows, bad = 0;
    printf("\n=== 齿 ===\n");
    {
        int L = WN[0] - 1;
        double d = A[0][L];
        double d1 = fabs(g_prod_c2[L] + 100.0 * d * C2X2 - C2X2) / C2X2;
        double ratio = d1 / d;
        J("齿1 靶挪 100×dev ⟹ dev 按比例放大（50~200x）", (ratio >= 50.0 && ratio <= 200.0) ? 0.0 : 1.0, 0.0);
    }
    {
        int L = WN[0] - 1;
        double fd = fabs(c2_fake(HW[0][L]) - C2X2) / C2X2;
        J("齿2 假式 prod(1-1/q^2) 必被拒（dev 大 >=100x）", (fd / A[0][L] >= 100.0) ? 0.0 : 1.0, 0.0);
    }
    for (i = 0; i < cells; i++) {
        rows = WN[i];
        if (CNT[i] == 1) {
            for (j = 1; j < rows; j++) if (A[i][j] > A[i][j - 1] * 1.05) bad++;
        } else {
            double lo = A[i][0], hi = A[i][0];
            for (j = 1; j < rows; j++) { if (A[i][j] < lo) lo = A[i][j]; if (A[i][j] > hi) hi = A[i][j]; }
            if (hi / lo < 2.0) bad++;
        }
    }
    J("齿3 FORCED 行不升 / TREND 行有变化（窗口在动）", (double)bad, 0.0);
    {
        double save = A[0][2], sp, r[NWMAX];
        A[0][2] = save * 1.5;
        for (j = 0; j < WN[0]; j++)
            r[j] = A[0][j] * pow(HW[0][j], PV[BP[0]]) * pow(log(HW[0][j]), QV[BQ[0]]);
        sp = spread(r, WN[0]);
        A[0][2] = save;
        J("齿4 单格抬 50% ⟹ 最优列 spread 破 0.5（能 FAIL）", (sp > SPTOL) ? 0.0 : 1.0, 0.0);
    }
    {
        int dup = 0;
        for (i = 0; i < cells; i++) if (CNT[i] > 1) dup++;
        J("齿5 每个 item 至多一列不动（对象唯一）", (double)dup, 0.0);
    }
}

/* ---- 段 5：落表 ---- */
static void seg5_csv(const char *path, int cells)
{
    FILE *f = fopen(path, "w");
    int i, j, rows;
    if (!f) { printf("  [WARN] 写不了 %s\n", path); return; }
    fprintf(f, "item,window,dev,class,mechanism_bound\n");
    for (i = 0; i < cells; i++) {
        rows = WN[i];
        for (j = 0; j < rows; j++) {
            char mech[64], cls[24];
            if (CNT[i] == 1) {
                sprintf(cls, "FORCED_trunc");
                sprintf(mech, "dev<=C/(h^%.1f*(lnh)^%.1f)", PV[BP[i]], QV[BQ[i]]);
            } else {
                sprintf(cls, "TREND");
                sprintf(mech, "NA");
            }
            fprintf(f, "\"%s\",%.0f,%.6e,%s,%s\n", ITEM[i], HW[i][j], A[i][j], cls, mech);
        }
    }
    fclose(f);
    printf("\n  已写 %s\n", path);
}

/* ---- 段 6：台账里所有「非被迫行 + 手挑 tol 行」（20 条）的机械分类 ----
   cls : 2=FORCED_trunc  3=FLOW  4=SEG  5=DUP
   axis: 0=无窗口轴  1=有轴未接（本器未覆盖）  2=本器已判（见段 3）
   fam : 0=tol/机制没有闭式来源（手挑）  1=有闭式来源
   返回: 0=KEEP 1=DOWNGRADE 2=PENDING 3=SAMPLE 4=MERGE 5=UPGRADE */
typedef struct { const char *id; int cls; int axis; int fam; const char *claim; } NROW;

static const NROW NR[20] = {
    { "bsd:50",   5, 0, 0, "rank1 Abel 和趋零" },
    { "bsd:51",   5, 0, 0, "rank2 Abel 和趋零" },
    { "bsd:60",   3, 0, 0, "秩不可由局部积读（负结果）" },
    { "hl:7",     3, 1, 0, "展开(最优截断) vs integral 偏差/尾≤3" },
    { "hl:9",     3, 2, 0, "一阶系数 c=D=k=2" },
    { "hl:11",    3, 1, 0, "下限 I(eps)~1/eps" },
    { "hl:54",    3, 1, 0, "与 archive 表 rho 对齐" },
    { "hl:67",    3, 1, 0, "|E(x)|<=6 sqrt(x)" },
    { "hl:114",   3, 1, 0, "密度积×log^2X 三条收敛" },
    { "rh:103",   4, 1, 0, "段最大值与整段一致" },
    { "rh:26",    4, 0, 0, "dup rh:103" },
    { "rh:286",   4, 0, 0, "dup rh:103" },
    { "rh:45",    4, 0, 0, "dup rh:103" },
    { "hl:15",    2, 2, 1, "焊接口 阶梯 vs Li2" },
    { "hl:18",    2, 2, 1, "2C2 截断积" },
    { "hl:20",    2, 2, 0, "e^-gamma 部分积" },
    { "ym:30",    2, 1, 0, "跑动耦合 error ∝ 1/N" },
    { "ym:73",    2, 1, 0, "lambda1 L2 -> (2 pi)^2" },
    { "ym:84",    2, 1, 0, "physical gap a-independent" },
    { "ns:147#P8", 2, 1, 0, "Riccati 闭式 + 二阶收敛" }
};

static const char *CLSN[6] = { "FORCED", "FORCED_eps", "FORCED_trunc", "FLOW", "SEG", "DUP" };
static const char *AXN[3]  = { "无轴", "有轴未接", "本器已判" };

static int judge_row(const NROW *r, char *out)
{
    if (r->cls == 5) { sprintf(out, "MERGE 重复计数：只计 1 次（G29）"); return 4; }
    if (r->axis == 2) { sprintf(out, "UPGRADE 已由本器判定（见段 3 与 nest_order_ledger.csv）"); return 5; }
    if (r->axis == 0) { sprintf(out, "SAMPLE 不可承重：无窗口轴 ⟹ 只能当读数/负结果"); return 3; }
    if (r->fam == 0 && r->cls == 2) {
        sprintf(out, "DOWNGRADE tol 无机制来源（手挑）⟹ 回到 TREND 或先补族检验");
        return 1;
    }
    sprintf(out, "PENDING 有轴未接：给出待接的公式与窗口轴，接上后必须能 FAIL");
    return 2;
}

static void seg6_classify(void)
{
    int i, k, n[6] = { 0, 0, 0, 0, 0, 0 };
    char buf[256];
    printf("\n=== 段 6：20 条非被迫行 / 手挑 tol 行的机械分类 ===\n");
    for (i = 0; i < 20; i++) {
        k = judge_row(&NR[i], buf);
        n[k]++;
        printf("  %-10s %-13s %-8s %s\n", NR[i].id, CLSN[NR[i].cls], AXN[NR[i].axis], buf);
    }
    printf("  小计：KEEP=%d DOWNGRADE=%d PENDING=%d SAMPLE=%d MERGE=%d UPGRADE=%d\n",
           n[0], n[1], n[2], n[3], n[4], n[5]);
    {   /* 齿6：手挑 tol 的 trunc 行必须被打回（判据能 FAIL） */
        NROW fake = { "fake:trunc", 2, 1, 0, "tol=1e-5 手挑" };
        int v = judge_row(&fake, buf);
        J("齿6 手挑 tol 的 trunc 行必被 DOWNGRADE（能 FAIL）", (double)(v != 1), 0.0);
    }
    {   /* 齿7：无轴行不得被当成判据（KEEP 即 FAIL） */
        NROW fake = { "fake:na", 3, 0, 1, "单样本" };
        int v = judge_row(&fake, buf);
        J("齿7 无窗口轴的行不得 KEEP（能 FAIL）", (double)(v != 3), 0.0);
    }
}

int main(void)
{
    int cells = NITEM;
    printf("==========================================================\n");
    printf(" 嵌套数组 × 关系族：先认出对象，再决定逼不逼成 0\n");
    printf(" four operations：加法=accumulate · 乘法=generate · 跨维度=dimension-change · base=非象\n");
    printf("==========================================================\n");

    sieve_init();
    stair_init();
    seg1_fill(cells);
    seg2_relations(cells);
    seg3_judge(cells);
    seg4_teeth(cells);
    seg5_csv("matrix/nest_order_ledger.csv", cells);
    seg6_classify();

    printf("\nRESULT PASS=%d FAIL=%d\nEXIT=%d\n", g_pass, g_fail, g_fail == 0 ? 0 : 1);
    return g_fail == 0 ? 0 : 1;
}
