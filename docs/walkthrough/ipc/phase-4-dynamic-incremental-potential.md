# Phase 4：Barrier 进入 Dynamic Incremental Potential

前面几篇已经讲清 active set 和 barrier energy 的构造。本篇描述它们如何从 contact 子模块的对象，变成每帧 dynamic solve 真正消费的总势能组成项。

## 帧协议

每帧 dynamic loop 的固定顺序：

```text
clear → detect → rebuild → attach → solve → free
```

具体：

1. `clearGeneralImplicitForceModel()` + `clearAlphaTestFunc()` — 清掉上一帧
2. External / self handler 分别 `execute(...)` — 重建 active set
3. Active set 非空时构造 barrier energy — frame-local 对象
4. `addGeneralImplicitForceModel(...)` — 挂进 integrator
5. `doTimestep(...)` — 求解
6. `freeBuffer(...)` — 释放 barrier buffer

**配置跨帧持有，handler 跨帧持有，barrier object 和 buffer 每帧重建。**

## Barrier 挂接

External 和 self 走同一条接缝，区别仅在 energy 类型和 position callback：

```text
// External (active samples > 0 时)
extBarrierEnergy = externalContactHandler->buildBarrierEnergy(externalIpcDhat)
extBarrierEnergy->allocateBuffer()
extBarrierEnergy->setComputePosFunction(currentPositionFunc)
extBarrierEnergy->setBuffer(extBarrierBuffer)
extBarrierEnergy->setCoeff(externalIpcKappa)
intg->addGeneralImplicitForceModel(extBarrierEnergy, 0, 0)

// Self (active pairs > 0 时)
selfBarrierEnergy = selfCD->buildBarrierEnergy(selfIpcDhat)
selfBarrierEnergy->allocateBuffer()
selfBarrierEnergy->setToPosFunction(currentPositionFunc)
selfBarrierEnergy->setBuffer(selfBarrierBuffer)
selfBarrierEnergy->setCoeff(selfIpcKappa)
intg->addGeneralImplicitForceModel(selfBarrierEnergy, 0, 0)
```

没有 active set 的帧不会创建 barrier — 总势能结构根据当前几何动态决定。

## `addGeneralImplicitForceModel(...)` 做了什么

不只是把指针塞进数组。它会：

1. 记录 `PotentialEnergy` 到 `generalAdditionalForceModels`
2. `createHessian(...)` 建立 Hessian 模板
3. 创建 `K`/`K1`/`M` 缓存
4. 分配 `fint` 缓存
5. 标记 `generalForceModelChanged = true`

## `assembleImplicitModels()` 中的合并

Integrator 把 implicit model 分三类：`IMT_ELASTIC`、`IMT_SAME_TOPOLOGY`、`IMT_GENERAL`。Barrier 走 `IMT_GENERAL`。

当 `generalForceModelChanged` 时，`assembleImplicitModels()` 会：

- 把 general model 逐个压入 `implicitModelsAll`
- 建立 `Kmapping`（local Hessian → 系统主 Hessian 的映射）
- 合并 sparsity pattern 到 `hessianAll`

完成后，barrier 在 `updateA()`、`updateb()` 阶段与其它 implicit model 处于同一层级。

## `doTimestep(...)` 中的消费

`ImplicitBackwardEulerTimeIntegrator::doTimestep(...)` 的关键顺序：

```text
assembleImplicitModels() → updateD() → updateA() → updateb() → solve(...)
```

Barrier 已在 `implicitModelsAll` 中，因此其梯度进入总残差、Hessian 进入系统矩阵。Newton 求解面对的是包含 barrier 的一步总能量。

## 验证

Runtime smoke（`runSimConfig_parse_test.cpp`）覆盖 external-only、self-only、merged、deterministic 四种路径，证明 barrier 已进入一步求解。

---

上一阶段：[Phase 3](phase-3-self-near-contact-active-set.md)
下一阶段：[Phase 5](phase-5-feasible-line-search.md)
