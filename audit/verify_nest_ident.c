/* verify_nest_ident.c —— 嵌套数组进混沌一步到位：用不同关系函数「认出对象」
 *
 * 四法口径（four operations，凡跨维度必付残差）：加法=accumulate · 乘法=generate · 取模=sieve · 跨维度=dimension-change · base=非象
 *   跨维度处只有 1/N，残差口径 = 相对偏差。
 *
 * 方法：同一张嵌套数组 A[窗口][量]，换不同关系函数去读它 ——
 *   r1 = MAC / S_N         （和之比）
 *   r2 = MAC / (N·Gs)      （对不对得上线性）
 *   r3 = MAC × 1/N         （除回窗口）
 *   r4 = MAC − N·Gs        （相消后的残差）
 *   r5 = S_N / MAC         （反向比）
 * 哪个关系在**换窗口时不动**，那个不动的就是对象本身；哪个跟着窗口走，那就是窗口。
 * 判据不靠猜：不动的那几个关系直接给出「它是线性累积 + 地板」这句话。
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
#define NQ 4
static int    WIN[NW] = { 500, 2000, 5000 };
static double A[NW][NQ];   /* 0=S_N  1=MAC  2=N·Gs  3=Gs */

static void fill(double q, double Gs)
{
    int i, n;
    for (i = 0; i < NW; i++) {
        double SN = 0.0, S2N = 0.0, g = 1.0;
        for (n = 0; n < 2 * WIN[i]; n++) {
            if (n < WIN[i]) SN += g;
            S2N += g;
            g = fmax(q * g, Gs);
        }
        A[i][0] = SN;                  /* 微观侧：一段和 */
        A[i][1] = S2N - SN;            /* 宏观侧：两段之差（= 线性累积） */
        A[i][2] = (double)WIN[i] * Gs; /* 窗口 × 地板 */
        A[i][3] = Gs;                  /* 地板本身 */
    }
}

int main(void)
{
    const double q = 0.5, Gs = 0.1256637061;
    int i, k;
    double r2, r3, r4, r5, v1 = 0.0, v2 = 0.0, v3 = 0.0, v5 = 0.0, teeth = 0.0;

    fill(q, Gs);
    printf("=== 嵌套数组 × 五个关系函数（混沌一步到位）===\n");
    printf("  WIN      S_N        MAC        N·Gs       Gs\n");
    for (i = 0; i < NW; i++)
        printf("  %4d  %9.4f  %9.4f  %9.4f  %.6f\n", WIN[i], A[i][0], A[i][1], A[i][2], A[i][3]);

    /* ---- 判项块 A：换窗口不动的三个关系（= 对象身份） ---- */
    {
        for (i = 0; i < NW; i++) {
            r2 = A[i][1] / A[i][2];          /* 对不对得上线性 */
            r3 = A[i][1] * 1.0 / (double)WIN[i];
            r4 = A[i][1] - A[i][2];          /* 残差 */
            if (i == 0) { v1 = r2; v2 = r3; v3 = r4; }
            else {
                teeth += fabs(r2 - v1) + fabs(r3 - v2) + fabs(r4 - v3);
            }
            printf("  r2=MAC/(N·Gs)=%.9f  r3=MAC/N=%.9f  r4=MAC−N·Gs=%.3e\n", r2, r3, r4);
        }
        J("r2 = MAC/(N·Gs) 恒为 1（线性累积）", fabs(v1 - 1.0), 1e-12);
        J("r3 = MAC/N 恒为 Gs（除回窗口只剩地板）", fabs(v2 - Gs), 1e-12);
        J("r4 = MAC−N·Gs 恒为 0（相消完全）", fabs(v3), 1e-9);
        J("三个不动关系换窗口都不动（对象身份稳定）", teeth, 1e-9);
    }

    /* ---- 判项块 B：跟着窗口走的那一个（= 窗口本身） ---- */
    {
        double lo = 1e300, hi = 0.0;
        for (i = 0; i < NW; i++) {
            r5 = A[i][0] / A[i][1];          /* 反向比：随窗口变 */
            if (r5 < lo) lo = r5;
            if (r5 > hi) hi = r5;
            printf("  r5=S_N/MAC=%.9f（窗口 %d）\n", r5, WIN[i]);
        }
        v5 = (hi - lo) / hi;
        J("r5 随窗口变化（证明它不是不变量，只是窗口的函数）", (v5 > 1e-6) ? 0.0 : 1.0, 0.0);
    }

    /* ---- 判项块 C：齿 —— 把 MAC 抬 1e-3，三个不动关系必须破 ---- */
    {
        double broken = 0.0;
        for (i = 0; i < NW; i++) {
            double m = A[i][1] * (1.0 + 1e-3);
            broken += fabs(m / A[i][2] - 1.0) + fabs(m / (double)WIN[i] - A[i][3]) + fabs(m - A[i][2]);
        }
        printf("  齿：MAC 抬 1e-3 ⟹ 三个关系偏离合计 = %.6f\n", broken);
        J("齿：任一不动关系被扰动必须破", (broken > 1e-6) ? 0.0 : 1.0, 0.0);
    }

    /* ---- 判项块 D：认出对象（把上面结论写成一句可判的话） ---- */
    {
        int ok = 1;
        for (i = 0; i < NW; i++) {
            if (fabs(A[i][1] / A[i][2] - 1.0) > 1e-12) ok = 0;
            if (fabs(A[i][1] / (double)WIN[i] - A[i][3]) > 1e-12) ok = 0;
        }
        printf("  认定：该对象 = 「线性累积 × 地板」；窗口是唯一的自由度\n");
        J("对象认定：{r2≡1, r3≡Gs} 同时成立 ⟹ 线性累积 + 地板", ok ? 0.0 : 1.0, 0.0);
    }

    printf("\nRESULT PASS=%d FAIL=%d\nEXIT=%d\n", gp, gf, gf == 0 ? 0 : 1);
    return gf == 0 ? 0 : 1;
}
