/* verify_dual_teeth.c —— 对偶层的三个齿（T1 反演 / T2 换角色 / T3 不动点唯一）
 *
 * 四法口径（four operations，凡跨维度必付残差）：加法=accumulate · 乘法=generate · 取模=sieve · 跨维度=dimension-change · base=非象
 *   跨维度处只有 sinh/asinh/sqrt（用 log 与 sqrt 展开），残差口径 = 相对偏差。
 *
 * 学术依据：
 *   Kramers–Wannier： sinh(2K)·sinh(2K*) = 1，自对偶点 sinh(2Kc)=1 ⟹ Kc = ½·ln(1+√2) = 0.4406868…
 *   T-对偶： M²(R) = n²/R² + m²R²/α'²（α'=1），R ↔ 1/R 同时 (n,m) ↔ (m,n) ⟹ 谱不变
 *
 * 这三个齿正是 onestep 矩阵里 L5(verify_dual, 8/17) 与 L13(verify_dual_exact, 21/41) 缺的 29 格的判法。
 */
#include <stdio.h>
#include <math.h>

static int gp = 0, gf = 0;
static void J(const char *name, double dev, double tol)
{
    int ok = (dev <= tol);
    if (ok) gp++; else gf++;
    printf("  [%s] %-58s dev=%.6e tol=%.1e\n", ok ? "PASS" : "FAIL", name, dev, tol);
}

/* sinh / asinh 用 log 与 sqrt 展开（不引 math.h 的 sinh/asinh） */
static double mysinh(double x)
{
    double e = 1.0, t = 1.0; int i;
    for (i = 1; i <= 40; i++) { t *= x / (double)i; if (i % 2 == 1) e += t; else e += 0.0; }
    /* e 只累到奇数项：sinh(x) = Σ x^(2k+1)/(2k+1)!，用另一条更稳的路：exp 级数 */
    {
        double ep = 1.0, s = 1.0, em = 1.0, s2 = 1.0, tp = 1.0, tm = 1.0;
        for (i = 1; i <= 40; i++) {
            tp *= x / (double)i; s += tp;         /* e^x */
            tm *= -x / (double)i; s2 += tm;       /* e^{-x} */
        }
        ep = s; em = s2;
        (void)e;
        return 0.5 * (ep - em);
    }
}
static double myasinh(double y) { return log(y + sqrt(y * y + 1.0)); }
static double kw_map(double K) { return 0.5 * myasinh(1.0 / mysinh(2.0 * K)); }

int main(void)
{
    const double Kc = 0.5 * log(1.0 + sqrt(2.0));      /* 0.4406867935… */
    double Ks[4] = { 0.20, 0.35, 0.4406867935, 0.60 };
    int i;

    printf("=== 对偶层三个齿（Kramers–Wannier + T-对偶）===\n");
    printf("  自对偶点 Kc = ½ln(1+√2) = %.10f\n", Kc);

    /* ---- T1 反演：sinh(2K)·sinh(2K*) == 1；把 K* 换成 K 必须破 ==== */
    {
        double dev = 0.0, teeth = 0.0;
        for (i = 0; i < 4; i++) {
            double K = Ks[i], Ks2 = kw_map(K);
            dev += fabs(mysinh(2.0 * K) * mysinh(2.0 * Ks2) - 1.0);
            teeth += fabs(mysinh(2.0 * K) * mysinh(2.0 * K) - 1.0);   /* 错映射：K* := K */
        }
        J("T1 反演：sinh(2K)·sinh(2K*) = 1（四点）", dev, 1e-10);
        printf("     齿：把 K* 换成 K，偏离量 = %.6f（必须 > 0）\n", teeth);
        J("T1 齿：错映射必须破（K*:=K 时 ≠ 1）", (teeth > 1e-6) ? 0.0 : 1.0, 0.0);
    }

    /* ---- T2 不动点唯一：sinh(2Kc) == 1，且 Kc 处 K* == K ---- */
    {
        double dev = fabs(mysinh(2.0 * Kc) - 1.0);
        double fix = fabs(kw_map(Kc) - Kc);
        double off = 0.0;
        for (i = 0; i < 4; i++) off += fabs(mysinh(2.0 * (Ks[i] + 0.01)) - 1.0);
        J("T2 自对偶点：sinh(2Kc) = 1", dev, 1e-12);
        J("T2 不动点：K*(Kc) = Kc", fix, 1e-12);
        printf("     齿：Kc±0.01 处 sinh(2K)−1 偏离 = %.6f（必须 > 0）\n", off);
        J("T2 齿：偏离 Kc 必须破", (off > 1e-6) ? 0.0 : 1.0, 0.0);
    }

    /* ---- T3 T-对偶：M²(R) 在 R↔1/α'R 且 (n,m)↔(m,n) 下不变 ---- */
    {
        double dev = 0.0, teeth = 0.0;
        int nm[3][2] = { {1, 0}, {0, 1}, {2, 3} };
        double Rs[3] = { 0.5, 1.0, 2.0 };
        for (i = 0; i < 3; i++) {
            int k;
            for (k = 0; k < 3; k++) {
                double R = Rs[k], n = (double)nm[i][0], m = (double)nm[i][1];
                double M1 = n * n / (R * R) + m * m * R * R;
                double R2 = 1.0 / R;
                double M2 = m * m / (R2 * R2) + n * n * R2 * R2;     /* 换角色 + 换半径 */
                double M3 = m * m / (R * R) + n * n * R * R;         /* 只换角色，不换半径 */
                dev += fabs(M1 - M2);
                teeth += fabs(M1 - M3);
            }
        }
        J("T3 T-对偶：R↔1/R 且 (n,m)↔(m,n) 时 M² 不变（九点）", dev, 1e-12);
        printf("     齿：只换角色不换半径，偏离量 = %.6f（必须 > 0）\n", teeth);
        J("T3 齿：漏掉半径反演必须破", (teeth > 1e-6) ? 0.0 : 1.0, 0.0);
    }

    /* ---- T4 边界：本仓的 Gs 与 Kc / √α' 数值上无映射（照实标，不硬凑） ---- */
    {
        double Gs = 0.1256637061;
        printf("  边界：Gs=%.10f  vs  Kc=%.10f  vs  √α'=1.0000000000 ⟹ 无映射（不宣称同一个不动点）\n", Gs, Kc);
        J("T4 边界：不得把 Gs 说成 Kc 或 √α'（本条只作声明，不作等式）",
          (fabs(Gs - Kc) > 1e-3 && fabs(Gs - 1.0) > 1e-3) ? 0.0 : 1.0, 0.0);
    }

    printf("\nRESULT PASS=%d FAIL=%d\nEXIT=%d\n", gp, gf, gf == 0 ? 0 : 1);
    return gf == 0 ? 0 : 1;
}
