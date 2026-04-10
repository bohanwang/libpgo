# IPC Walkthrough

## 概述

本组文档描述 `contact-model = "ipc-barrier"` 在 `libpgo` dynamic runtime 中的完整实现。
阅读起点建议接在 [Simulation Tick: Pre-Tick Setup](../simulation-tick-pre-tick.md) 之后。

IPC 不是独立的求解器分支，而是挂接在 shared volumetric dynamic path 上：

- 配置层表达 frictionless barrier contact
- runtime 区分 penalty / IPC 两种 contact model
- external / self near-contact active set 按帧重建
- barrier energy 以 `PotentialEnergy` 身份加入 timestep 总势能
- Newton line search 尊重 contact-aware feasible alpha upper bound

## 数据流

一帧 `ipc-barrier` runtime 的执行链：

```text
parseRunSimConfig(...)
  ↓  读取 contact-model, dhat, kappa, alpha-safety, feasible-line-search
runSimFromConfig(...)
  ↓  isDynamicContactEnabled(...) → 决定是否创建 handler
帧循环开始
  ↓  clearGeneralImplicitForceModel() + clearAlphaTestFunc()
  ↓
  ├─ external: execute(usurf, externalIpcDhat)
  │    → active samples → PointPenetrationBarrierEnergy
  │    → addGeneralImplicitForceModel(...)
  │
  ├─ self: execute(usurf, selfIpcDhat)
  │    → BVH query band → seed refine → active pairs
  │    → PointTrianglePairBarrierEnergy
  │    → addGeneralImplicitForceModel(...)
  │
  └─ feasible line search (if enabled):
       → merged alpha callback = min(external, self)
       → setAlphaTestFunc(...)
  ↓
ImplicitBackwardEuler::doTimestep(...)
  ↓  assembleImplicitModels() → barrier 进入 implicitModelsAll
  ↓  updateA() / updateb() → barrier Hessian + gradient 参与总系统
  ↓
NewtonRaphsonSolver::solve(...)
  ↓  alpha = alphaTestFunc(x, dx)  // contact-aware upper bound
  ↓  alpha ≤ 0 → SR_FEASIBLE_STEP_ZERO
  ↓  否则从 alpha 开始 LSM_SIMPLE step acceptance
```

## 建模边界

- 仅覆盖 dynamic path（不展开 static）
- IPC 模式下 friction 被显式禁止
- self path 是 sample-based point-triangle，不是 primitive PT/EE
- feasible alpha 是 contact-aware upper bound，不含 inversion-free filter
- 跟的是 `SO_NEWTON` + `LSM_SIMPLE` 主路径

## Phase Map

| Phase | 职责 | 主要入口 |
| --- | --- | --- |
| [Phase 1](phase-1-config-and-dispatch.md) | 配置解析、参数校验、runtime 分发 | `parseRunSimConfig`, `isDynamicContactEnabled` |
| [Phase 2](phase-2-external-barrier-energy.md) | external active sample、barrier energy、alpha 上界 | `TriangleMeshExternalContactHandler`, `PointPenetrationBarrierEnergy` |
| [Phase 3](phase-3-self-near-contact-active-set.md) | self active pair、BVH query band、seed refine、self barrier | `TriangleMeshSelfContactDetection`, `TriangleMeshSelfContactHandler`, `PointTrianglePairBarrierEnergy` |
| [Phase 4](phase-4-dynamic-incremental-potential.md) | barrier 挂进 timestep 总势能 | `addGeneralImplicitForceModel`, `assembleImplicitModels` |
| [Phase 5](phase-5-feasible-line-search.md) | merged feasible alpha callback 透传到 Newton | `setAlphaTestFunc`, `NewtonRaphsonSolver` |
| [Phase 6](phase-6-tests-and-validation.md) | 验证闭环：数学单测 → handler 回归 → API smoke | `tests/*`, `examples/pulled-cubic-box-self-ipc` |

## 关键代码入口

**配置与 runtime：**

- `src/api/runSimCore.h` / `.cpp`

**Contact handler：**

- `src/core/solve/contact/triangleMeshExternalContactHandler.{h,cpp}`
- `src/core/solve/contact/triangleMeshSelfContactDetection.{h,cpp}`
- `src/core/solve/contact/triangleMeshSelfContactHandler.{h,cpp}`

**Barrier energy：**

- `src/core/solve/contact/pointPenetrationBarrierEnergy.{h,cpp}`
- `src/core/solve/contact/pointTrianglePairBarrierEnergy.{h,cpp}`
- `src/core/solve/nonlinearOptimization/barrierFunction.{h,cpp}`

**Solver：**

- `src/core/solve/simulation/timeIntegrator.{h,cpp}`
- `src/core/solve/simulation/timeIntegratorSolver.cpp`
- `src/core/solve/simulation/implicitBackwardEulerTimeIntegrator.cpp`
- `src/core/solve/nonlinearOptimization/NewtonRaphsonSolver.{h,cpp}`
