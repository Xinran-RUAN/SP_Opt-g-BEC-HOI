# 运行入口说明

所有命令均从项目根目录运行。先执行：

```matlab
startup_HOI
```

脚本名称保持不变，以免破坏论文、MAT 文件和外部笔记中的路径。这里按用途分类，而不是物理移动文件。

## 常用入口

| 目的 | 命令 | 说明 |
| --- | --- | --- |
| 单个 1D 解 | `run('runs/run_single_ground_state.m')` | 最适合改一个参数并快速检查 |
| 全部低成本测试 | `run('tests/run_all_tests.m')` | 修改代码后的首选验证 |
| Fisher 正则化收敛 | `run('runs/run_regularization_comparison.m')` | 固定空间网格的正式实验 |
| 最终势能平滑实验 | `run('runs/run_potential_regularization_effects_L32.m')` | L=32 的 baseline 与固定 sigma 比较 |
| 正式 Fourier 网格实验 | `run('runs/run_mesh_refinement_spectral_harmonic.m')` | 谐振子核心加边界 C-infinity 延拓 |
| 两阶段求解器性能 | `run('runs/RUN_TWO_STAGE_SOLVER_PERFORMANCE.m')` | 支持脚本顶部的 smoke 开关 |
| 共享前缀性能图 | `run('runs/RUN_TWO_STAGE_SOLVER_PERFORMANCE_SHARED_PREFIX.m')` | FISTA 公共前缀与 Newton 分支 |
| 两个 2D 示例 | `run('runs/RUN_TWO_DIMENSIONAL_EXAMPLES.m')` | 支持脚本顶部的 smoke 开关 |

## 论文作图与主实验

- `run_plot_fourier_spatial_convergence.m`：生成 Fourier 空间收敛图；
- `run_plot_sigma_spectral_accuracy.m`：比较固定问题下 sigma 对空间精度的影响；
- `run_potential_regularization_effects_L32.m`：最终势能正则化尺度实验；
- `run_mesh_refinement_spectral_harmonic.m`：正式光滑周期外势网格实验；
- `RUN_TWO_STAGE_SOLVER_PERFORMANCE.m`：求解器成本表和主性能数据；
- `RUN_TWO_STAGE_SOLVER_PERFORMANCE_SHARED_PREFIX.m`：共享 FISTA 前缀的图；
- `RUN_TWO_DIMENSIONAL_EXAMPLES.m`：各向异性谐振子与光晶格示例。

这些脚本通常网格较大、会复用检查点或写入 `results/`/`figs/`。运行前先检查脚本顶部参数和复用开关。

## 基础比较与回归入口

- `run_solver_comparison.m`：一阶求解器和公共精修比较；
- `run_mesh_refinement.m`：通用 eta=0 谱空间网格诊断；
- `run_potential_regularization_comparison.m`：固定网格 P0/P1/P2 比较；
- `run_potential_mesh_refinement.m`：指定势能变体的网格诊断；
- `run_potential_pipeline_single.m`：低成本 FISTA 到 KKT 精修检查；
- `run_potential_splitting_comparison.m`：旧/新 FISTA 分裂 A/B；
- `run_polish_linear_solver_comparison.m`：同一交接点下 GMRES 与 PCG-Schur；
- `run_gradient_consistency.m`：梯度测试的便捷入口；
- `run_projection_equivalence.m`：两种投影实现的便捷入口。

## 可选研究分支

- `run_entropy_active_set_sweep.m`：熵、活跃集和 Fourier 尾部；
- `run_entropy_bias_resolution_tradeoff.m`：固定网格的偏差/分辨率重叠区；
- `run_entropy_mesh_refinement.m`：只应在选定 eta 后运行；
- `run_potential_epsilon_sweep.m`：手动的小 epsilon 刚性诊断；
- `run_spatial_operator_consistency.m`：逐项 Fourier 一致性，不调用优化器。

## 明确标记为诊断的入口

- `run_DIAG_sigma_layer_resolution.m`；
- `run_DIAG_periodic_potential_mesh.m`；
- `run_DIAG_fixed_sigma_periodic_mesh.m`；
- `run_DIAG_harmonic_boundary_sensitivity.m`。

这些脚本用于解释数值现象，不改变正式求解器。结果不能在没有检查问题定义、边界外势和参考状态资格的情况下直接当作论文主结论。

## 兼容入口

`run_potential_sigma_sweep.m` 是旧名称，只转发到 `run_potential_regularization_effects_L32.m`。新说明和新调用应直接使用后者；保留前者是为了不破坏旧笔记和调用命令。

## 后处理

`post/` 中的脚本只读取已有 MAT 文件并作图，不调用求解器：

- `post_mesh_refinement.m`；
- `post_regularization_comparison.m`；
- `post_solver_comparison.m`；
- `diagnose_figure55_pre_switch.m`。

