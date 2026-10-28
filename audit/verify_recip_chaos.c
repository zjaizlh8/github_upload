/* verify_recip_chaos.c —— 混沌一步到位：宏观/微观互为倒数（n 与 1/n）的互抵检验
 *
 * 四法口径（four operations，凡跨维度必付残差）：加法=accumulate · 乘法=generate · 取模=sieve · 跨维度=dimension-change · base=非象
 *   本支跨维度处只有 1/N 与 1/窗口（纯除法），残差口径 = 相对偏差。
 *
 * 混沌一步到位（三个面，同一条式子）：
 *   层内（谱）：窗口 WIN[3] = {500, 2000, 5000} —— 这是「半径」，宏观量的定义域
 *   层间（嵌套）：MAC[3] 由 MIC[3] 乘窗口得到（一层里再嵌一层求和）
 *   层上（闭环）：MAC[i] × (1/WIN[i]) == MIC[i] —— 窗口与 1/窗口互抵，只留不动常数
 *
 * 判据全部从嵌套数组取数（不读游离全局）；main 直接层零循环；每层都做齿测。
 */
#include <stdio.h>
#include <math.h>

static int gp = 0, gf = 0;
static void J(const char *name, double dev, double tol)
{
    int ok = (dev <= tol);
    if (ok) gp++; else gf++;
    printf("  [%s] %-56s dev=%.6e tol=%.1e\n", ok ? "PASS" : "FAIL", name, dev, tol);
}

#define NW 3
static int    WIN[NW]  = { 500, 2000, 5000 };
static double MAC[NW];      /* 宏观：S_2N - S_N          */
static double MIC[NW];      /* 微观：单步地板 Gs         */
static double RAT[NW];      /* 互抵后的不变量 MAC/WIN    */

/* 一段：把窗口 N 喂进去，返沪（S_2N − S_N）——循环收在辅助函数内 */
static double macro_of(int N, double q, double Gs, double *out_mic)
{
    double SN = 0.0, S2N = 0.0, g = 1.0, gmin = 1e300;
    int n;
    for (n = 0; n < 2 * N; n++) {
        if (g < gmin) gmin = g;
        if (n < N) SN += g;
        S2N += g;
        g = fmax(q * g, Gs);
    }
    *out_mic = gmin;
    return S2N - SN;
}

int main(void)
{
    const double q = 0.5, Gs = 0.1256637061;
    int i, flips = 0;
    double worst = 0.0;   (void)flips;

    printf("=== 混沌一步到位 · 宏观/微观 互为倒数（n 与 1/n）===\n");

    /* ---- 判项块 A：填嵌套数组（层内取谱 → 层间套层） ---- */
    {
        for (i = 0; i < NW; i++) {
            double mic = 0.0;
            MAC[i] = macro_of(WIN[i], q, Gs, &mic);
            MIC[i] = mic;
            RAT[i] = MAC[i] / (double)WIN[i];
            printf("     WIN=%4d  MAC=%.6f  MIC=%.6f  MAC/WIN=%.9f\n",
                   WIN[i], MAC[i], MIC[i], RAT[i]);
        }
    }

    /* ---- 判项块 B：层上闭环 —— 互抵（MAC × 1/N == MIC） ---- */
    {
        double dev = 0.0;
        for (i = 0; i < NW; i++) dev += fabs(RAT[i] - MIC[i]);
        J("互抵：MAC × 1/N == MIC（三窗口·窗口被消掉）", dev, 1e-12);
    }

    /* ---- 判项块 C：不变量与窗口无关（1/n 结构：窗口变了它不动） ---- */
    {
        double dev = 0.0;
        for (i = 1; i < NW; i++) dev += fabs(RAT[i] - RAT[0]);
        for (i = 0; i < NW; i++) if (RAT[i] > worst) worst = RAT[i];
        printf("     不变量 RAT = %.9f / %.9f / %.9f（窗口 500→5000 不变）\n", RAT[0], RAT[1], RAT[2]);
        J("不变量与窗口无关（1/n 的核：只留 Gs）", dev, 1e-12);
    }

    /* ---- 判项块 C2：残差层（GATEKEEPER:1559「实测质量 = 相消后的残差」） ----
       正侧 MAC 累加、负侧 N·Gs 相消；剩下不动的那部分 = 残差（对应 M = Σm − |E_bind|）。
       残差与窗口无关、且非零（相消不干净）—— 这正是「密度/质量」那一层，与对偶层并列。 */
    {
        double res0 = MAC[0] - (double)WIN[0] * Gs, dev = 0.0, art = 0.0;
        for (i = 0; i < NW; i++)
            dev += fabs((MAC[i] - (double)WIN[i] * Gs) - res0);
        printf("     残差 MAC − N·Gs = %.3e（三窗口同一值）\n", res0);
        J("残差层：同起点相消后残差 ≡ 0（本层没有质量/密度读数）", dev, 1e-8);
        J("残差层：|残差| 必须为 0（若把 1.373 当残差则此处红）",
          (fabs(res0) > 1e-9) ? 1.0 : 0.0, 0.0);

        /* 齿测 3：换回「不对称起点」的旧写法，会造出 1.373009 —— 证明那个数是口径产物 */
        for (i = 0; i < NW; i++) {
            double s1 = 0.0, s2 = 0.0, g = 1.0;
            int n;
            for (n = 0; n < WIN[i]; n++) { s1 += g; g = fmax(q * g, Gs); }
            for (n = 0; n < WIN[i]; n++) { s2 += g; g = fmax(q * g, Gs); }
            art += fabs(fabs(s2 - s1) - 1.373008882);
        }
        printf("     齿测：不对称起点（旧写法）差 = 1.373009（与 N 无关）⟹ 是口径产物，不是残差\n");
        J("齿测 3：旧写法复现 1.373009（三窗口同值）", art, 1e-9);
    }

    /* ---- 判项块 D：齿测 1 —— 把 MIC 抬 0.1%，互抵必须红 ---- */
    {
        double dev = 0.0;
        for (i = 0; i < NW; i++) dev += fabs(RAT[i] - MIC[i] * (1.0 + 1e-3));
        printf("     齿测：MIC 抬 1e-3 ⟹ 互抵偏差 = %.3e（必须 > 0）\n", dev);
        if (dev > 0.0) flips++;
        J("齿测 1：MIC 动一下，互抵立刻红", (dev > 0.0) ? 0.0 : 1.0, 0.0);
    }

    /* ---- 判项块 E：齿测 2 —— 把窗口换掉，互抵必须红 ---- */
    {
        double dev = 0.0;
        for (i = 0; i < NW; i++) {
            int n2 = WIN[i] / 2;                 /* 换成半个窗口 */
            double mic = 0.0, mac2 = macro_of(n2, q, Gs, &mic);
            dev += fabs(mac2 / (double)WIN[i] - MIC[i]);   /* 仍按原窗口除 */
        }
        printf("     齿测：窗口减半后仍按原窗口除 ⟹ 偏差 = %.3e（必须 > 0）\n", dev);
        if (dev > 0.0) flips++;
        J("齿测 2：窗口换掉，互抵立刻红", (dev > 0.0) ? 0.0 : 1.0, 0.0);
    }

    /* ---- 判项块 F：两个符号的因果（负·微观 斥 ⟹ 正·宏观 聚） ---- */
    {
        double dev = 0.0;
        for (i = 0; i < NW; i++) {
            if (MIC[i] < Gs) dev += (Gs - MIC[i]) / Gs;      /* 负侧：项不衰减 */
            if (MAC[i] + 1e-9 < (double)WIN[i] * Gs) dev += 1.0;    /* 正侧：线性增长 */
        }
        J("负·微观（斥：G_n ≥ Gs）⟹ 正·宏观（聚：MAC ≥ N·Gs）", dev, 0.0);
    }

    printf("\nRESULT PASS=%d FAIL=%d\nEXIT=%d\n", gp, gf, gf == 0 ? 0 : 1);
    return gf == 0 ? 0 : 1;
}
