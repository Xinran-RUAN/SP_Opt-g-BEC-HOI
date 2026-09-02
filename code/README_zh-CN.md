# HOI 密度表述的凸优化实现

[English](README.md) | 中文说明 | [中文代码结构说明](docs/ARCHITECTURE.md) | [运行入口说明](runs/README.md)

本目录是当前论文使用的、自包含的 MATLAB 包化实现。主要数学模型和数值实验是一维问题；两个代表性的二维算例复用同一套向量优化流程。

- 本文件介绍数学模型、求解器设计和主要运行方式。
- `docs/ARCHITECTURE.md` 说明代码分层、数据结构、数值不变量和修改风险边界。
- `runs/README.md` 按用途列出全部运行入口。

`legacy/` 中的日期快照只用于追溯，不加入 MATLAB 路径。`results/` 和 `figs/` 保存数值数据与图件，它们是输出目录，不是另一套源代码。

## 数学模型

在周期配置点网格

```matlab
h = 2*L/N;
x = -L + (0:N-1)'*h;
```

上，代码求解以下离散能量的最小化问题：

\[
E_h(\rho)=h\sum_j\left[
\frac{(D\rho)_j^2}{8s_\varepsilon(\rho_j)}+V_j\rho_j
+\frac{\beta}{2}\rho_j^2+\frac{\delta}{2}(D\rho)_j^2
\right],
\]

约束为配置点上的非负性和固定质量：

\[
\rho_j\ge 0,\qquad h\sum_j\rho_j=1.
\]

主线算法是带正则化的凸密度表述，并通过 Fourier 伪谱离散及保持非负性和质量的优化方法求解。线性势能模型仍是基线；其中没有双调和惩罚项，也没有谱滤波。

与论文直接对应的核心接口为：

```text
s_epsilon, ds_epsilon, d2s_epsilon
p_sigma,   dp_sigma,   d2p_sigma
```

`Energy`、`Gradient`、矩阵自由 Hessian、有限差分预条件器以及势能 prox 都直接调用这些函数句柄。简单公式集中写在正式运行脚本的顶部。

当前仍保留下列保持凸性的 Fisher 内置构造：

- `shift_smooth`：\(s=\rho+\varepsilon\)；
- `piecewise_c1`、`piecewise_c2`、`piecewise_c3`：由 `src.regularization.EvaluateDenominator` 实现的凹过渡公式。

旧名称只会在求解开始前转换一次，统一成相同的函数句柄接口；核心求解器内部不会把名称当作数学分支开关。

## 目录结构

```text
+model/                    默认参数、1D/2D 网格、外势和初值
+src/+regularization/      Fisher 分母定义、适配与校验
+src/+potential/           线性势能和凸平滑势能映射
+src/+discretization/+ps/  Fourier plan、D、D^T、能量、梯度与网格传递
+src/+constraints/         质量/可行性和两种独立投影实现
+src/+entropy/             熵值、导数、simplex prox 与 KKT
+src/+diagnostics/         活跃集、支撑集、远场和 Fourier 尾部诊断
+src/+solvers/             ISTA、SPG、FISTA-CD、残差与两类精修器
+experiments/              配置、问题装配、状态传递、计时与保存格式
runs/                      论文实验、比较、验证与诊断入口
post/                      只读取已保存结果的后处理脚本
tests/                     低成本数值验证
results/, figs/            已保存数据和生成图件
docs/                      架构与维护说明
legacy/                    日期化历史快照，不加入 MATLAB path
```

更详细的调用关系和各层职责见 `docs/ARCHITECTURE.md`。

## 运行方式

在本目录中启动 MATLAB，先初始化一次路径：

```matlab
startup_HOI
```

最短的常用流程为：

```matlab
run('runs/run_single_ground_state.m')
run('tests/run_all_tests.m')
```

每个运行脚本都把可编辑参数集中放在顶部。只修改当前实验对应的值，例如 `N`、`regularization_name`、`epsilon`、`solver_name` 或 `projection_name`。

当前最终势能平滑实验的直接入口是：

```matlab
run('runs/run_potential_regularization_effects_L32.m')
```

旧名称 `run_potential_sigma_sweep.m` 只作为兼容转发器保留。大规模论文实验、二维算例、熵分支和诊断脚本的用途与前置条件见 `runs/README.md`。

## 优化算法结构

- `ISTA` 是参考投影梯度方法，用于正确性验证和小规模比较。
- `SPG` 是默认的一阶段凸优化器，采用带保护的 Barzilai–Borwein 步长和非单调投影线搜索。
- `FISTA-CD` 是可选的一阶加速比较方法，不是默认生产求解器。由于算法包含可行性处理和单调重启，这里不声称完整的、无重启 FISTA 收敛率。
- `PDAS / semismooth Newton` 是可选的高精度矩阵自由 KKT 精修阶段，在一阶方法进入局部区域后使用。

SPG 名称中的 “spectral” 指 Barzilai–Borwein 谱步长，与 Fourier 伪谱空间离散无关。

所有方法都使用同一个固定步长投影梯度残差进行停止和比较，即 `residual_step = 1`。能量平台或状态平台只是诊断量，不是驻点证书。默认流程为：

```text
SPG（PG 约 1e-8） -> 可选 PDAS 精修（PG 约 1e-12）
```

## 势能正则化实验

基线势能项 \(V\rho\) 在固定 \(\varepsilon\) 的障碍问题中可能产生人为真空活跃区，因此势能实验将它与

\[
V p_\sigma(\rho),\qquad
p_\sigma(\rho)=\sqrt{\rho^2+\sigma^2}-\sigma
\]

比较。正式运行脚本可以直接定义：

```matlab
epsilon = 1e-3;
s_epsilon   = @(rho) rho + epsilon;
ds_epsilon  = @(rho) ones(size(rho));
d2s_epsilon = @(rho) zeros(size(rho));

sigma = epsilon^4;
p_sigma = @(rho) rho.^2 ./ (hypot(rho,sigma) + sigma);
dp_sigma = @(rho) rho ./ hypot(rho,sigma);
d2p_sigma = @(rho) ...
    (sigma./hypot(rho,sigma)).^2 ./ hypot(rho,sigma);
```

固定网格尺度实验支持 1–4 次幂。旧名称 `sqrt_same_scale` 和 `sqrt_squared_scale` 分别保留为 1 次幂和 2 次幂的别名。

由于 \(p_\sigma'(0)=0\)，并且 \(p_\sigma\) 为凸函数，当 \(V\ge0\) 时，修改后的势能项保持凸性；求解器会检查这一条件。能量、梯度、Hessian 和 prox 使用完全相同的函数句柄。代码中的稳定数值公式避免了相消误差，但标签仍可保留论文中的平方根表达式。函数句柄和标签都会写入 MAT 结果。

`src.potential.Expression` 和 `src.regularization.KineticExpression` 只是说明/兼容辅助函数，不是正式求解中数学公式的来源。

全项目采用统一的 KKT 符号约定：

\[
G(\rho)+\lambda\mathbf 1-\mu=0.
\]

对于严格正的状态，投影后的活跃点数量只用作诊断。interior Newton-PCG 的适用性由严格正性、真空斜率 \(p_\sigma'(0)\) 和强凸性保护条件共同判断；不能仅因某个投影点出现零值，就把具有极小正尾部的状态送入 PDAS。

固定网格 A/B/C 实验对三个变体使用相同的正初值、FISTA-CD 设置和 Bo Lin 半光滑正质量投影。该实验关闭熵项且不做精修，用于检验真空斜率平滑能否消除人为活跃集并改善 Fourier 精度；它并不预设谱精度一定恢复。目标能量、线性势能基线偏差、远场、Fourier 尾部和一阶求解成本分别报告。

对于平方根变体，正式 FISTA-CD 分裂把势能项与非负性、质量守恒共同放入可分离的正质量 prox。光滑部分的回溯只控制 Fisher、beta 和 delta 项。最终仍用完整目标的 `Energy`、`Gradient` 以及独立的 Bo Lin 全梯度映射进行认证。

`legacy_full_gradient` 仅保留给局部分裂 A/B 检查：

```matlab
run('runs/run_potential_splitting_comparison.m')
```

## 消失熵选项（实验性诊断分支）

可选项

\[
\eta H_h(\rho)=\eta h\sum_j\rho_j(\log\rho_j-1)
\]

是计算正则化，不是新的物理项。设置 `entropy.enabled = false` 或 `entropy.eta = 0` 会恢复原始凸密度问题及其单纯形投影路径。

当 `eta>0` 时，一阶阶段使用带 entropy-simplex 复合 prox 的 FISTA-CD（或参考 ISTA）。可选的内部等式约束 Newton 精修不使用 PDAS 活跃集分类器，问题仍保持凸性。

该分支用于检验消除人为活跃集能否改善 Fourier 空间正则性，但不预设熵一定恢复谱精度。熵偏差与固定 eta 的空间误差必须分开研究；解释 Fourier 收敛前还必须检查有限区间尾部。经典 SPG 和 PDAS 始终是 `eta=0` 基线。

## 网格加密与参考解

默认网格实验使用：

```matlab
N_list = [32 64 128 256 512];
N_ref = 1024;
```

主要谱状态比较为：

```text
rho_N 与 P_N rho_ref 比较
```

其中 \(P_N\) 是把参考解 Fourier 系数正交限制到嵌套粗谱空间的算子。程序分别报告已解析模态误差、由 Parseval 关系直接得到的遗漏参考尾部，以及二者的正交总误差。不会对 \(P_N\rho_{ref}\) 额外施加非负投影或质量投影。把 \(\rho_N\) 谱延拓到最细网格只作为次级可视化和比较诊断。

原生网格上的非线性能量差与受保护的公共网格能量诊断分开保存。如果未经改变的过采样三角插值出现非正值，则公共网格熵能量无效，程序不会通过裁剪把它强行变成有效值。

尾部质量和尾部最大值分别保存，避免把有限区间平台误判成 Fourier 误差。只有当最细状态的优化残差和 Fourier 尾部都通过指定充分性阈值时，才称其为参考解。

完整的通用网格实验需显式运行：

```matlab
run('runs/run_mesh_refinement.m')
```

在选择生产用熵参数前，先运行固定网格权衡实验：

```matlab
run('runs/run_entropy_bias_resolution_tradeoff.m')
```

它寻找同时满足物理能量偏差要求和相对密度过渡层网格分辨率要求的 eta。只有找到重叠区并选定 eta 后，才应运行 `run_entropy_mesh_refinement.m`。

Fisher 正则化实验固定空间网格，并从较大 epsilon 到较小 epsilon 使用简单 warm start：

```matlab
run('runs/run_regularization_comparison.m')
```

网格结果和正则化结果分别保存。`post/` 中的文件只加载已有 MAT 文件并作图，不会调用求解器。光滑正则化和分段正则化表现出的谱行为在这里属于数值研究，不作为定理陈述。

## 测试

运行全部低成本测试：

```matlab
run('tests/run_all_tests.m')
```

测试覆盖 Fourier 微分及其斜伴随性、零填充延拓、参考谱精确限制、Parseval 已解析/尾部误差分解、投影等价性、中心差分能量–梯度一致性、矩阵自由 Hessian 一致性、随机离散凸性、SPG 可行性与下降性、FISTA 到精修的交接、interior Newton-PCG 与 PDAS 一致性、熵 prox，以及二维能量–梯度–Hessian 一致性。

## 维护原则

不能只因为两个实现看起来相似，就合并核心数值代码。独立投影路径和兼容适配器也是验证体系的一部分。

以下修改必须给出数学等价理由并运行全部测试：

- 能量、梯度与 Hessian 作用；
- Fourier 归一化和导数/伴随算子；
- 投影与 prox；
- 求解器迭代次序和一阶到精修的交接；
- 停止准则、容差和默认参数。

详细风险分级、核心数值不变量和兼容入口清单见 `docs/ARCHITECTURE.md`。
