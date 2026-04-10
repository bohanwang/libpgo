# Phase 5：feasibility-preserving line search

## 1. 这阶段要回答什么

这一阶段专门回答下面四个问题：

1. 当前 repo 里的 feasible alpha callback 究竟从哪里安装、在哪里被消费。
2. `ipc-alpha-safety` 在当前实现里到底影响什么。
3. external / self 两条上界为什么能在 runtime 合成一个 merged alpha upper bound。
4. Newton solver 当前实际尊重的是哪条 line search 语义。

如果只看 parser 里的 `ipc-enable-feasible-line-search`，会误以为这只是一个预留开关。  
但当前仓库里，这条 callback 已经真正接进了 dynamic Newton 主路径。

## 2. 这阶段在整条 IPC 主线里解决什么

Phase 4 已经让 barrier 进入了总势能。  
但仅有 barrier 还不够，因为 solver 仍然可能沿着一个让接触距离继续恶化的方向走 full step。

Phase 5 当前真正解决的是：

> 在现有 Newton 框架不重写的前提下，给每一步搜索方向加上一个 contact-aware feasible step upper bound。

因此这一阶段讲的不是 “line search 的一般理论”，而是当前 repo 里：

- contact 层怎样给出 `alphaUpper`
- integrator 怎样透传这个 callback
- solver 怎样从这个上界开始做 step acceptance

## 3. 当前 repo 已经完成了什么

当前仓库里，这一阶段已经具备下面这些事实：

- `ipcEnableFeasibleLineSearch` 已进入 `RunSimContactConfig`，默认值为 `true`
- `runSimCore.cpp` 已会在 barrier 模式下按当前 active set 安装 alpha callback
- external handler 与 self handler 都已提供：
  - `computeEmbeddedAlphaUpperBound(...)`
- runtime 已把 external / self 上界合成为一个 merged callback
- `TimeIntegrator` 已公开：
  - `setAlphaTestFunc(...)`
  - `clearAlphaTestFunc()`
- `TimeIntegratorSolver` 已把 callback 透传给 `NewtonRaphsonSolver`
- `NewtonRaphsonSolver::solve(...)` 已会：
  - 先取 `alphaUpper`
  - `alphaUpper <= 0` 时返回 `SR_FEASIBLE_STEP_ZERO`
  - 否则再进入 line search / step acceptance

这说明 feasible alpha 在当前 repo 里已经不是纸面接口，而是真正参与 dynamic IPC runtime 的 solver 逻辑。

## 4. callback 当前在什么条件下安装

### 4.1 parser 层的总开关

当前是否允许安装 callback，首先取决于：

- `contact-model == "ipc-barrier"`
- `ipc-enable-feasible-line-search == true`

runtime 中对应的布尔量是：

```cpp
const bool shouldUseIpcAlphaFilter =
    contact.contactModel == "ipc-barrier" && contact.ipcEnableFeasibleLineSearch;
```

因此 penalty 模式完全不会走这条路径。

### 4.2 frame 级别还要看 active set 是否存在

当前实现还会进一步检查：

- `hasExternalAlphaUpperBound`
- `hasSelfAlphaUpperBound`

只有至少一条接触路径当前帧真的有 active set 时，runtime 才会安装 callback。  
也就是说：

- barrier 模式但这一帧没有 active sample / active pair
  不会无意义安装 callback

这让 feasible alpha 保持了典型的 frame-local 语义。

## 5. merged callback 当前怎样组织

### 5.1 callback 的几何语义不写在 solver 里

当前 solver 本身并不知道：

- external active sample 是什么
- self active pair 是什么
- `dhat` 和 `dSafe` 分别是什么

这些几何语义都保留在 contact/runtime 层，由 `runSimCore.cpp` 组装出一个 merged callback 再传给 solver。  
这让职责划分很清楚：

- contact 层负责几何可行性
- solver 只负责消费一个 `alphaUpper`

### 5.2 callback 输入的不是绝对位移，而是 Newton 状态

在 `ImplicitBackwardEuler` 里，Newton 优化变量是增量 `z`。  
因此 runtime callback 当前先做：

```cpp
ES::VXd currentU = u + z;
```

然后再用：

- `currentU`
- `dz`

去查询 external / self 的上界。

这一步非常关键，因为它说明当前 feasible alpha 检查的对象是：

```text
当前 frame 的绝对位移状态 + 当前 Newton 搜索方向
```

而不是把 solver 内部变量直接当世界坐标使用。

### 5.3 callback 先把 `dhat` 转成 `dSafe`

当前 runtime 安装 callback 前，会分别计算：

- `externalDSafe = PointPenetrationBarrierEnergy::normalizePositiveZero(externalIpcDhat, -1.0)`
- `selfDSafe = PointTrianglePairBarrierEnergy::normalizePositiveZero(selfIpcDhat, -1.0)`

因此 callback 真正消费的不是裸 `dhat`，而是：

- 与 barrier 数值安全域对齐后的 `dSafe`

这保证了 feasible alpha 与 barrier 在 safe margin 上使用同一套数值语义。

### 5.4 external / self 两条上界在 runtime 合并

当前 callback 内部会分别计算：

- `externalAlphaUpper`
- `selfAlphaUpper`

然后做：

```cpp
alphaUpper = std::min(alphaUpper, externalAlphaUpper);
alphaUpper = std::min(alphaUpper, selfAlphaUpper);
```

这意味着 merged IPC 场景在 solver 看来只是一个最小上界问题。  
solver 不需要理解 “external 优先还是 self 优先”，runtime 已经把它压缩成一个最终上界。

## 6. `ipc-alpha-safety` 当前到底影响什么

当前代码里，`ipc-alpha-safety`：

- 不参与 active set
- 不参与 barrier energy
- 不参与 `kappa`

它只在 external/self handler 计算 alpha upper bound 时起作用。  
也就是说，它当前是：

> solver-side feasibility margin

而不是 contact stiffness 或 activation distance 的别名。

从语义上看：

- `alphaSafety = 1.0`
  使用当前 conservative 几何上界
- `alphaSafety < 1.0`
  在这个上界基础上进一步收紧一步

这也是相关单测里检查 “alpha 随 safety 线性缩放” 的原因。

## 7. callback 当前怎样进入 solver

### 7.1 `TimeIntegrator` 只是保存 callback

`TimeIntegrator` 当前公开了：

- `setAlphaTestFunc(...)`
- `clearAlphaTestFunc()`

它并不在这里解释几何语义，而只是把 callback 存起来。

### 7.2 `TimeIntegratorSolver` 负责往下透传

真正求解时，`TimeIntegratorSolver::solve(...)` 会把 `alphaTestFunc` 一路传给：

- `NewtonRaphsonSolver`

这一步的意义是：IPC 不需要旁路一个专用 Newton；它复用了现有 integrator solver 接缝。

### 7.3 `NewtonRaphsonSolver` 当前怎样消费这个上界

在 `NewtonRaphsonSolver::solve(...)` 里，当前关键逻辑是：

```cpp
alpha = alphaTestFunc ? alphaTestFunc(x, deltax) : 1.0;
```

之后分成两步：

1. 如果 `alpha <= 0.0`
   直接返回 `SR_FEASIBLE_STEP_ZERO`
2. 否则从这个 `alpha` 开始继续 line search

因此当前 solver 已经明确把 feasible alpha upper bound 当成 step acceptance 的前置约束。

## 8. 当前默认 Newton line search 语义是什么

### 8.1 当前默认主路径仍然是 `LSM_SIMPLE`

`NewtonRaphsonSolver::SolverParam` 当前默认 line search method 是：

```cpp
LineSearchMethod lsm = LSM_SIMPLE;
```

因此 walkthrough 中这条 IPC runtime 主线，默认是在 `LSM_SIMPLE` 语义下运行。

### 8.2 `LSM_SIMPLE` 当前真的会从 `alphaUpper` 开始试步

在 `LSM_SIMPLE` 分支下，solver 会：

1. 先用 callback 给出的 `alpha` 作为初始步长
2. 如果能量不下降，再执行 `alpha *= 0.75` 的回缩

因此当前 repo 的 feasible line search 语义可以概括成：

> 先用 contact callback 给出几何上界，再在这个上界之内做默认 Newton step acceptance。

### 8.3 `alphaUpper <= 0` 现在是可观测 solver 状态

如果 callback 返回零，solver 会返回：

- `SR_FEASIBLE_STEP_ZERO`

并且 `statusToString(...)` 会把它转成：

- `feasible_step_zero`

这说明当前 feasible alpha 不只是影响一个局部变量，而是已经成为 solver 状态机的一部分。

## 9. runtime 日志当前怎样暴露 callback 行为

当前 merged callback 会打印下面几类日志：

- `IPC feasible alpha callback active.`
- `IPC feasible alpha merged callback active.`
- `IPC feasible alpha upper bound: ... (external: ..., self: ...)`

这些日志有两个作用：

1. 对本地调试很直接，可以看出 callback 是否安装、是哪条路径在收紧步长。
2. API smoke 会直接依赖这些日志来断言当前 runtime 的 IPC callback 语义已经生效。

## 10. 当前 repo 的验证证据

`tests/core/scene/contact_embedding_test.cpp` 已经固定 external alpha 上界的关键语义：

- `FeasibleStepUpperBoundMatchesAnalyticPlaneBound`
- `FeasibleStepUpperBoundScalesWithAlphaSafety`
- `FeasibleStepUpperBoundIsOneForMotionAwayFromContact`
- `FeasibleStepUpperBoundReturnsOneWithNoActiveSamples`
- `FeasibleStepUpperBoundReturnsZeroWhenAlreadyInsideSafeMargin`
- `FeasibleStepUpperBoundAllowsRecoveryWhenInsideSafeMarginButMovingAway`

`tests/core/scene/self_contact_handler_test.cpp` 则固定 self alpha 上界的关键语义：

- `AlphaUpperBoundShrinksForInwardDirections`
- `AlphaUpperBoundScalesWithAlphaSafety`
- `AlphaUpperBoundAllowsOutwardRecovery`
- `AlphaUpperBoundBlocksFurtherInwardMotionInsideSafeBand`
- `AlphaUpperBoundReturnsOneWithNoActivePairs`

再往上，`tests/api/runSimConfig_parse_test.cpp` 的 smoke 还会验证：

- runtime 确实安装了 callback
- merged 场景确实安装了 merged callback
- 日志中确实存在 `< 1` 的 feasible alpha upper bound

## 11. 这一阶段不展开什么

当前 walkthrough 还不展开下面这些内容：

- inversion-free feasibility filter
- full CCD-based exact global upper bound
- primitive-level PT/EE self IPC
- `LSM_GOLDEN` / `LSM_BRENTS` / `LSM_BACKTRACK` 的完整重构

当前仓库文档只描述已经落地的默认 Newton 主路径，也就是：

- `SO_NEWTON`
- 默认 `LSM_SIMPLE`
- contact-aware alpha upper bound callback

上一阶段： [Phase 4](phase-4-dynamic-incremental-potential.md)  
下一阶段： [Phase 6](phase-6-tests-and-validation.md)
