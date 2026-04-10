# Phase 4：把 barrier 接进 dynamic incremental potential

## 1. 这阶段要回答什么

这一阶段专门回答下面四个问题：

1. 在当前 repo 里，barrier 究竟是怎样进入 timestep 总势能的。
2. 为什么 external / self barrier 都通过 `addGeneralImplicitForceModel(...)` 接进 solver。
3. `TimeIntegrator` 当前怎样把这些 general model 装配进 `ImplicitBackwardEuler`。
4. 当前 frame-local barrier 对象的生命周期怎样管理。

前面几篇已经把 external / self 两条 active-set 与 barrier 前提讲清楚。  
Phase 4 讲的则是：

> 这些 barrier 怎样从 contact 子模块里的对象，变成每一帧 dynamic solve 真的会消费的总势能组成项。

## 2. 这阶段在整条 IPC 主线里解决什么

如果没有这一阶段，前面的 external/self active set 和 barrier energy 仍然只停留在 contact 层。  
真正让 IPC 进入 shared runtime 的，是下面这套 frame protocol：

- 清掉上一帧的 contact runtime
- 用当前位形重建 active set
- 按 active set 新建 barrier energy
- 把 barrier 作为 general `PotentialEnergy` 加进 integrator
- 再执行一步 `ImplicitBackwardEuler` 求解

因此 Phase 4 的本质不是再设计 barrier 公式，而是：

> 把 barrier 变成 current time step incremental potential 的正式组成项。

## 3. 当前 repo 已经完成了什么

当前仓库里，这一阶段已经具备下面这些事实：

- dynamic loop 已形成固定的 `clear -> detect -> rebuild -> attach -> solve` 协议
- external / self barrier 都走同一条 general `PotentialEnergy` 接缝
- `TimeIntegrator::addGeneralImplicitForceModel(...)` 已经能为 barrier 建立：
  - Hessian 模板
  - K/K1/M 缓存
  - `fint` 缓存
- `assembleImplicitModels()` 已经会把这些模型作为 `IMT_GENERAL` 合入 `implicitModelsAll`
- `ImplicitBackwardEulerTimeIntegrator::doTimestep(...)` 已经把它们当成一步总势能的一部分来求解
- external/self barrier buffer 的分配、绑定、释放已经在 `runSimCore.cpp` 中形成稳定生命周期

这意味着 barrier 当前不再只是 contact 模块中的几何对象，而是已经成为 dynamic one-step solve 的正式输入。

## 4. dynamic frame 当前的运行骨架

### 4.1 每一帧先清理旧 contact runtime

当前 dynamic 帧循环开始时，`runSimCore.cpp` 会先做：

```cpp
intg->clearGeneralImplicitForceModel();
intg->clearAlphaTestFunc();
```

这说明当前 repo 的 IPC 设计不是跨帧持有 barrier object，而是：

- 每帧重新检测 active set
- 每帧重新创建 barrier energy
- 每帧重新决定是否安装 feasible alpha callback

也就是说，当前 barrier 是典型的 frame-local rebuild 资源。

### 4.2 `u`、`usurf` 和 `currentPositionFunc` 一起构成 runtime 的坐标接缝

当前 dynamic path 里，体网格位移 `u` 会先通过 `W` 映到 surface displacement `usurf`。  
然后 contact handler 使用 `usurf` 做检测，barrier energy 再通过：

- `makeCurrentPositionFunction(restPosition)`

把当前 solver 变量恢复成实际位置。

这一步把：

- surface-side detection
- volumetric DOF-side energy evaluation

连成了同一条 runtime 链。

## 5. external / self barrier 当前怎样被挂接

### 5.1 external barrier 的挂接顺序

当前 external 路径在 `contact-model == "ipc-barrier"` 时，会执行：

1. `externalContactHandler->execute(usurf.data(), externalIpcDhat)`
2. 如果 `getNumActiveSamples() > 0`
3. `buildBarrierEnergy(externalIpcDhat)`
4. `allocateBuffer()`
5. `setComputePosFunction(currentPositionFunc)`
6. `setBuffer(extBarrierBuffer)`
7. `setCoeff(externalIpcKappa)`
8. `intg->addGeneralImplicitForceModel(extBarrierEnergy, 0, 0)`

这说明 external barrier 当前不是“总是存在，只是有时值为零”，而是：

- 只有当前帧 active sample 非空时才创建
- 每帧按当前 active set 重建

### 5.2 self barrier 的挂接顺序

当前 self 路径在 `contact-model == "ipc-barrier"` 时，会执行：

1. `selfCD->execute(usurf.data(), selfIpcDhat)`
2. 打印 `# self active pairs: ...`
3. 如果 `getNumActivePairs() > 0`
4. `buildBarrierEnergy(selfIpcDhat)`
5. `allocateBuffer()`
6. `setToPosFunction(currentPositionFunc)`
7. `setBuffer(selfBarrierBuffer)`
8. `setCoeff(selfIpcKappa)`
9. `intg->addGeneralImplicitForceModel(selfBarrierEnergy, 0, 0)`

external 和 self 在这一步的核心差异只有 energy 类型与 position callback 名称不同；它们共享同一条 integrator 接缝。

### 5.3 没有 active set 就不会创建 barrier

当前一步总势能中到底有没有 barrier，取决于本帧真实检测结果：

- external active sample 是否非空
- self active pair 是否非空

因此 runtime 当前不是“配置一开就无脑多加两项势能”，而是根据 frame-local geometry 动态决定总势能结构。

## 6. barrier object 在进入 integrator 前要完成什么绑定

在调用 `addGeneralImplicitForceModel(...)` 之前，external / self barrier 都要完成下面这些 runtime 绑定：

- `allocateBuffer()`
- `setBuffer(...)`
- `setCoeff(...)`
- 绑定当前位置函数

external 使用：

```cpp
extBarrierEnergy->setComputePosFunction(currentPositionFunc);
```

self 使用：

```cpp
selfBarrierEnergy->setToPosFunction(currentPositionFunc);
```

这一层绑定的含义是：

- barrier 类型本身只知道自己是 `PotentialEnergy`
- runtime 还需要为它补齐当前帧的 buffer、系数和位置恢复逻辑

完成这些绑定后，它才真正成为本帧可求值的 energy object。

## 7. `addGeneralImplicitForceModel(...)` 当前到底做了什么

### 7.1 这不是简单把指针塞进数组

`TimeIntegrator::addGeneralImplicitForceModel(...)` 当前会做下面几件事：

1. 把 `PotentialEnergy` 记录到 `generalAdditionalForceModels`
2. 调 `fm->createHessian(...)` 建立 Hessian 模板
3. 创建：
   - `generalAdditionalForceModels_K`
   - `generalAdditionalForceModels_K1`
   - `generalAdditionalForceModels_M`
4. 为内部力分配 `generalAdditionalForceModels_fint`
5. 根据 integrator 的阻尼配置保存 damping 参数
6. 标记 `generalForceModelChanged = true`

因此 barrier 一旦被加入 integrator，后续装配就会把它当成正规的 implicit model，而不是 runtime 特判项。

### 7.2 `M` 缓存会沿系统质量矩阵模板初始化

当前实现还会根据 `MasK` 把 general model 的 `M` 缓存填起来。  
这意味着 general barrier model 在 integrator 看来并不是“裸 Hessian”，而是已经适配了现有系统矩阵组装框架的模型条目。

## 8. `assembleImplicitModels()` 当前怎样把 barrier 合进去

### 8.1 integrator 把模型分成三类

当前 `TimeIntegrator` 会把所有 implicit model 分成：

- `IMT_ELASTIC`
- `IMT_SAME_TOPOLOGY`
- `IMT_GENERAL`

external/self barrier 走的都是第三类 `IMT_GENERAL`。

### 8.2 `generalForceModelChanged` 为真时会重建 mapping

当这一标志为真时，`assembleImplicitModels()` 会：

- 把 general model 逐个压入 `implicitModelsAll`
- 为每个 model 建立 `generalAdditionalForceModels_Kmapping`
- 把所有 general Hessian 的 sparsity pattern 合并进 `hessianAll`
- 再用 `small2Big(...)` 建立到系统主 Hessian 的映射

这一步的意义在于：

- barrier 自己维护自己的 Hessian 模板
- integrator 负责把它们映到系统统一矩阵结构中

### 8.3 barrier 从这一刻开始就是 `implicitModelsAll` 的一员

一旦 `assembleImplicitModels()` 完成，barrier 在后续 `updateA()`、`updateb()` 和内部力更新阶段，就和其它 implicit model 处于同一层级。

这正是 walkthrough 中“接进 dynamic incremental potential”的代码级含义。

## 9. `ImplicitBackwardEuler` 当前在哪里真正消费 barrier

当前 `ImplicitBackwardEulerTimeIntegrator::doTimestep(...)` 的关键顺序是：

1. `assembleImplicitModels()`
2. `updateD()`
3. `updateA()`
4. `updateb()`
5. 调 `TimeIntegratorSolver::solve(...)`

由于 barrier 已经在 `implicitModelsAll` 中，所以从这一步开始：

- barrier 梯度会进入一步总残差
- barrier Hessian 会进入系统矩阵
- Newton 求解面对的是包含 barrier 的一步总能量

这一步是当前 IPC 与 “只做外部几何检测但不进 solver” 之间的真正分界线。

## 10. 当前 barrier 的生命周期为什么说是 frame-local

当前 `runSimCore.cpp` 中，barrier buffer 的使用顺序是：

1. 检测到 active set 后再分配 buffer
2. 本帧求解前把 buffer 绑定到 energy
3. `doTimestep(...)` 完成后立即 `freeBuffer(...)`

因此当前 repo 的 barrier 生命周期非常清晰：

- 配置是全局的
- handler 是跨帧持有的
- barrier object 与 barrier buffer 是每帧局部资源

这与 active set 必须按当前几何重建的语义完全一致。

## 11. 当前 repo 的验证证据

这一阶段最直接的回归不是纯数学测试，而是 runtime smoke。  
`tests/api/runSimConfig_parse_test.cpp` 中下面这些用例共同证明了 barrier 已经进入 timestep：

- `RunSimFromConfigCubicDynamicIpcNearContactActivatesExternalBarrier`
- `RunSimFromConfigCubicDynamicIpcSelfBarrierSmokeTest`
- `RunSimFromConfigCubicDynamicIpcMergedBarrierSmokeTest`
- `RunSimFromConfigCubicDynamicIpcMergedDeterministicSmokeTest`

这些 smoke 共同说明：

- barrier 已经进入一步求解
- external/self 两条 barrier 都能独立挂接
- merged runtime 路径也已存在
- deterministic mode 下这条路径的可观测序列可以稳定复现

## 12. 这一阶段还不展开什么

这一阶段还不展开下面这些问题：

- alpha upper bound 怎样进入 Newton
- static path
- constraints solver 路径
- persistent contact state

这些问题里：

- solver-side feasible alpha callback 留给 [Phase 5](phase-5-feasible-line-search.md)
- 其余 static / constraints / persistent-state 相关内容不在当前 IPC walkthrough 展开

上一阶段： [Phase 3](phase-3-self-near-contact-active-set.md)  
下一阶段： [Phase 5](phase-5-feasible-line-search.md)
