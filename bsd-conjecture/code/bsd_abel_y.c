/* bsd_abel_y.c —— Abel 和的 (M, y) 极限检验：TREND 还是 trunc，由这一次测量决定
 *
 * 层=过程（本支只报读数，不作结构判定）
 *
 * 【要判什么】
 *   `verify_bsd.c` 里那三条「Abel 和趋零」用的是**绝对阈值** `|S| < 0.02`，
 *   而代码自己写明「绝对值受核污染，只读零/非零」⟹ 阈值不是机制给的。
 *   本支不改 `verify_bsd.c`（判据不动），只把 `S1_abel_L1` 在四组 (M, y) 上量出来，
 *   看 rank≥1 的残差是**随核尾衰减**（⟹ trunc，有界可算）还是**不动**（⟹ 另有来源）。
 *
 * 【五法口径（前四法沿用，第五法新增）】加法=维度内对消 · 乘法=维度外张成 · 取模=降维约束 · 指数对数=跨维度付残差 · 无损归心=可逆回 1
 *   跨维度（指数对数）：Abel 核 `exp(−2πny/√N)`、`sqrt` 取 √N ⟹ **必付残差**，残差口径 = |ΔS| 读数
 *   取模(金)=降维约束：**本支不涉**（不打包、不切周期）⟹ 周期闭式 = 不涉 · 余项有界 = 0（恒等）
 *   乘法：`a_n/n` 逐项乘核 · 加法：Kahan 补偿求和（沿用 `verify_bsd.c` 的实现）
 *   无损归心：本支无有损映射 ⟹ R(R(x)) − x = 0（平凡成立，无取模可逆）
 *
 * 【机制（若成立则为截断）】
 *   核尾 ~ e^{−2πMy/√N} ⟹ 令 τ = 2πMy/√N，残差应随 τ 增大而指数下降。
 *
 * 【停点两套（S8）】
 *   连续侧（求和）：n 上限 = M；低于 → 计入 · 高于 → 不计（±0）
 *   离散侧（读数）：四组 (M, y) = (8000,1e-3) (8000,4e-3) (32000,1e-3) (32000,4e-3)
 *   对照：rank 0（11a1）必须**不随 y 掉到零**（L(11a1,1)≠0）
 *
 * 【G32 三问】
 *   容器 = 曲线 a_n（n ≤ M，M ≤ 32000）+ 四个 (M, y) 档
 *   空间维度 = 2（rank 通道，核档）
 *   dc = 值域基数：S 的读数 [1e−6, 1] 按 1e−6 量化 ⟹ 1e6
 *   对称 = y ↦ 0 的极限（Abel 定理）；无平移对称
 *   闭式 = 无；本支只报读数
 *
 * 【G34 三关系】
 *   关系1：S 随 y 变（对象：核污染）
 *   关系2：S 随 M 变（对象：截断）
 *   关系3：rank 0 与 rank≥1 的 S 随 y 是否分道（对象：秩判据）
 *
 * 【F / R 双向】
 *   F: 正向 —— 四档 (M,y) 各算 S，看衰减。
 *   R: 反向 —— 齿：① M 加倍而 y 不变 ⟹ S 必须变；② y 放大 4 倍 ⟹ rank≥1 的 |S| 必须更小；③ rank 0 不得跟着掉到 0。
 *
 * 【四判据数】res = |S(8000,1e-3) − S(32000,1e-3)|（同核换 M，实算）· domain = 4 档 · dc = 1000000 · order = 2
 * 【五行六步】土S 锚定 → 火F 虚化 → 金M 坍缩 → 水W 渗透 → 木D 扩展 → 土S 闭环
 *
 * 契约九行：SUBPROBE/D1/D3/D3_TOTAL/S2/ROUTE/CLOSED/DOMAIN/GATE
 */
#define main verify_bsd_main_unused
#include "verify_bsd.c"
#undef main

#include <stdio.h>
#include <math.h>

static const int    MS[2] = { 8000, 32000 };
static const double YS[2] = { 1e-3, 4e-3 };

static void abel_sweep(void)
{
    double S[3][2][2];
    int c, i, j, teeth = 0;
    double r1, r2, r3;

    for (c = 0; c < 3; c++)
        for (i = 0; i < 2; i++)
            for (j = 0; j < 2; j++)
                S[c][i][j] = S1_abel_L1(&S1_CURVES[c], MS[i], YS[j]);

    printf("SUBPROBE=bsd_abel_y\n");
    printf("  曲线         M=8000,y=1e-3   M=8000,y=4e-3   M=32000,y=1e-3  M=32000,y=4e-3\n");
    for (c = 0; c < 3; c++)
        printf("  %-10s  %+.9f     %+.9f     %+.9f     %+.9f\n",
               S1_CURVES[c].name, S[c][0][0], S[c][0][1], S[c][1][0], S[c][1][1]);

    r1 = fabs(S[1][0][0] - S[1][1][0]) + fabs(S[2][0][0] - S[2][1][0]);
    r2 = fabs(S[1][0][0]) / fmax(fabs(S[1][0][1]), 1e-300);
    r3 = fabs(S[0][0][0]) / fmax(fabs(S[1][0][0]) + fabs(S[2][0][0]), 1e-300);

    printf("齿1 M 加倍（y 不变）⟹ |ΔS| = %.3e（须 > 0）\n", r1);
    printf("齿2 y 放 4 倍 ⟹ |S(1e-3)|/|S(4e-3)| = %.3f（rank≥1 须 > 1）\n", r2);
    printf("齿3 rank0 / rank≥1 = %.1f 倍（须 > 10 ⟹ 两支分道）\n", r3);
    printf("核尾指标 τ = 2πMy/√N（M=8000,y=1e-3）：37a1 τ=%.2f · 389a1 τ=%.2f\n",
           2 * M_PI * 8000 * 1e-3 / sqrt((double)S1_CURVES[1].cond),
           2 * M_PI * 8000 * 1e-3 / sqrt((double)S1_CURVES[2].cond));

    printf("D1=%.0f\n", (r1 > 0) ? 0.0 : 1.0);
    printf("四判据数：res=%.3e · domain=4 · dc=1000000 · order=2\n", r1);
    printf("RESIDUAL=%.3e\n", r1);
    if (r1 > 0) teeth++;
    if (r2 > 1.0) teeth++;
    if (r3 > 10.0) teeth++;
    printf("D3=%d\n", teeth);
    printf("D3_TOTAL=3\n");
    printf("S2=%d\n", (teeth == 3) ? 1 : 0);
    printf("ROUTE=0\n");
    printf("CLOSED=1\n");
    printf("DOMAIN=读数层：Abel 和的 (M,y) 依赖；不含判读改动\n");
    printf("GATE=%s\n", (teeth == 3) ? "PASS" : "FAIL");
}

int main(void) { abel_sweep(); return 0; }
