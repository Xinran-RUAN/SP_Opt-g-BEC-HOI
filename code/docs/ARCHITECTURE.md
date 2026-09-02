# 代码结构与维护说明

本文档说明当前实现的分层、主要数据结构和修改边界。论文对应的数学模型与实验解释见根目录 `README.md`，全部运行入口见 `runs/README.md`。

## 1. 总体调用关系

```text
runs/*.m
  -> experiments.DefaultConfig / experiments.SolveGroundState
       -> model.*                         网格、初值、外势
       -> src.regularization.*            Fisher 分母及其导数
       -> src.potential.*                 势能正则化及复合 prox
       -> src.SolveGroundState1D          求解流程调度
            -> src.solvers.*              一阶迭代与 Newton/PDAS 精修
            -> src.discretization.ps.*    能量、梯度、Hessian 作用、FFT 算子
            -> src.constraints.*          非负与质量约束
            -> src.entropy.*              可选熵正则分支
       -> src.diagnostics.*               误差与谱诊断
  -> experiments.SaveResult               统一保存结果

post/*.m                                   只读取结果并作图
tests/*.m                                  小规模数值一致性检查
```

运行脚本负责声明论文中的具体问题；包目录负责实现可复用算法。正式实验不应在 `runs/` 中复制能量、梯度、投影或 Newton 公式。

## 2. 各层职责

| 路径 | 职责 | 是否属于数值核心 |
| --- | --- | --- |
| `+model/` | 默认参数、1D/2D 网格、初值和外势 | 部分属于 |
| `+experiments/` | 把配置装配为 `problem`，调用核心求解器，保存统一结果 | 装配层 |
| `+src/+discretization/+ps/` | Fourier 伪谱算子、离散能量和梯度 | 是 |
| `+src/+solvers/` | ISTA、SPG、FISTA-CD、Newton-PCG、PDAS 和残差 | 是 |
| `+src/+constraints/` | 单纯形/正质量投影及可行性检查 | 是 |
| `+src/+regularization/` | Fisher 分母接口、内置族和校验 | 是 |
| `+src/+potential/` | 势能正则化接口、校验和复合 prox | 是 |
| `+src/+entropy/` | 可选熵项及 entropy-simplex prox | 是（实验分支） |
| `+src/+diagnostics/` | 活跃集、远场、Fourier 尾部和网格比较 | 否，不应改变求解轨迹 |
| `runs/` | 单次运行、论文实验、验证和诊断入口 | 否 |
| `post/` | 读取已保存数据并作图 | 否 |
| `tests/` | 低成本数值回归测试 | 否 |
| `results/`, `figs/` | 已保存数据和论文图，不参与算法调用 | 否 |
| `legacy/` | 历史快照和迁移说明，不加入 MATLAB path | 否 |

## 3. 主要数据结构

### `config`

由 `experiments.DefaultConfig` 生成，运行脚本只覆盖本实验需要的字段。主要分组为：

- `parameters`：`beta`、`delta`、质量、区域与网格规模；
- `trapping_potential`：外势函数、标签及边界延拓元数据；
- `fisher_regularization` / `regularization`：Fisher 分母及其一、二阶导数；
- `potential_regularization`：`p_sigma` 及其一、二阶导数；
- `entropy`：可选熵项；
- `solver`：主求解器、分裂、停止准则与精修设置；
- `output`：绘图、保存和覆盖策略。

### `problem`

由 `experiments.SolveGroundState` 装配。它包含离散网格、Fourier plan、外势向量、物理参数和已经校验过的函数句柄。核心求解器只读取该结构，不应在迭代中重新解释实验名称。

### `result`

核心状态为 `rho` 和 `energy`；同时保存 `problem`/`solver` 元数据、迭代 `history`、停止与 KKT `diagnostics`、目标能量和基线能量。`experiments.SaveResult` 把常用字段写成稳定的 MAT schema，供 `post/` 使用。

## 4. 核心数值不变量

修改核心代码时必须同时保持以下约束：

1. 离散质量满足 `h*sum(rho)=mass`（2D 为面积权重）。
2. 节点密度保持非负；不能用绘图时裁剪代替求解器可行性。
3. `Energy`、`Gradient` 和 Hessian 作用使用同一组正则化函数句柄。
4. KKT 号约定统一为 `G(rho) + lambda*1 - mu = 0`。
5. 统一用固定步长的 projected-gradient residual 比较求解器。
6. `potential_prox` 分裂只把约定的势能项放入 prox；完整目标仍用于最终认证。
7. 熵关闭或 `eta=0` 时保留原投影路径和原问题含义。
8. Fourier 网格传递与参考谱限制不能额外施加非负投影或裁剪。

## 5. 修改风险分级

### 低风险

- README、注释、运行入口目录说明；
- `post/` 绘图布局和不参与求解的标签；
- `src.diagnostics` 中新增只读诊断；
- 删除 `.asv`、`.DS_Store` 等自动生成文件。

至少检查相关脚本可解析，并确认保存字段名称未被误改。

### 中风险

- `experiments.DefaultConfig` 的默认值；
- `experiments.SolveGroundState` 的装配和兼容转换；
- `experiments.SaveResult` 的 MAT schema；
- `model` 中的网格、外势与初值。

除运行全部测试外，还要检查旧结果/检查点兼容性和默认单次运行。

### 高风险（核心）

- `src.discretization.ps.Energy`、`Gradient`、导数/伴随算子；
- `src.solvers.FullHessianAction`、`HessianAction` 和预条件器；
- `src.constraints`、`src.potential.PositiveConservativeProx`、`src.entropy.SimplexProx`；
- `FISTACD`、`SPG`、`PolishKKT`、`PolishInteriorNewtonPCG`；
- `src.SolveGroundState1D` 中的求解器调度、停止和精修交接。

这类修改不能只凭代码“更简洁”合并。必须给出数学等价理由，并运行全部测试；涉及迭代次序、FFT 归一化、容差或默认值时，还应保存修改前后的同一小规模算例，比较能量、状态、质量误差、PG/KKT 残差和迭代历史。

## 6. 兼容入口

以下文件虽短或内部引用较少，但承担公共/历史兼容职责，不能按静态引用次数直接删除：

- `src.solvers.PolishPDAS`：显式选择 PDAS-GMRES 的兼容入口；
- `experiments.CompareStates`：旧的延拓到参考网格的次级诊断；
- `src.constraints.IsFeasible`：供外部脚本调用的可行性检查；
- `runs/run_potential_sigma_sweep.m`：转发到当前 L=32 最终势能平滑实验；
- 正则化和势能中的 legacy name adapters：用于读取旧配置和旧结果。

若将来确定不再需要兼容性，应先搜索论文脚本、外部笔记和已保存配置，再分一个独立变更删除，不与算法修改混在一起。

## 7. 验证命令

在 MATLAB 中：

```matlab
startup_HOI
run('tests/run_all_tests.m')
```

测试覆盖 Fourier 导数与伴随、一致性/凸性、投影和 prox、完整 Hessian、小问题求解、FISTA 到精修交接、PDAS 与 interior Newton，以及 2D 能量/梯度/Hessian。

