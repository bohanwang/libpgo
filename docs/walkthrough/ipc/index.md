# IPC Walkthrough

## 1. 这组 walkthrough 要回答什么

这一组文档专门回答下面几个问题：

1. `contact-model = "ipc-barrier"` 在当前 `libpgo` 仓库里到底已经接到了哪一层。
2. `parseRunSimConfig(...)` 之后，external / self IPC 参数怎样进入 dynamic runtime。
3. active set、barrier energy、feasible alpha callback 分别由哪些类负责。
4. 当前 repo 里的 IPC 能力边界在哪里，哪些行为已经被测试、smoke 和 example 固定下来。

阅读起点建议接在 [Simulation Tick: Pre-Tick Setup](../simulation-tick-pre-tick.md) 之后。  
前一篇已经把 mesh、embedding、deformation model、assembler 和 timestep 前置状态整理清楚；这一组文档从 contact runtime 开始继续，把 “一帧 dynamic solve 里 IPC 是怎么真正参与进去的” 讲完整。

## 2. 这组文档在整条 shared runtime 主线里解决什么

当前仓库的 IPC 不是单独的一条求解器分支，而是挂接在已有 shared volumetric dynamic path 上。  
从 `runSimCore.cpp` 的角度看，当前 `ipc-barrier` 主线解决的是：

- 配置层怎样表达 frictionless barrier contact
- runtime 怎样区分 penalty 和 IPC 两种 contact model
- external / self near-contact active set 怎样按帧重建
- barrier energy 怎样以 general `PotentialEnergy` 身份加入 timestep
- Newton line search 怎样额外尊重 contact-aware feasible alpha upper bound

因此 IPC walkthrough 讨论的不是一个独立 demo，而是 shared runtime 里 contact 这一层的实现形态。

## 3. 当前 repo 已经完成了什么

当前 `libpgo` 仓库里，`contact-model = "ipc-barrier"` 对应的 dynamic 路径已经具备下面这些能力：

- external near-contact active sample 检测
- self near-contact active pair 检测
- external barrier energy：`PointPenetrationBarrierEnergy`
- self barrier energy：`PointTrianglePairBarrierEnergy`
- external / self / merged feasible alpha upper bound callback
- external-only、self-only、ext+self merged 三类 runtime smoke
- deterministic merged smoke

从仓库真相出发，当前 IPC 已经不是“有 barrier 工具函数但 runtime 还没接”的状态，而是：

> `parseConfig -> detect -> barrier build -> addGeneralImplicitForceModel -> alpha callback -> Newton`

这条动态链路已经完整存在。

## 4. 当前 repo 的 phase map

| Phase | 当前阶段在 repo 里的职责 | 主要入口 | 主要验证 |
| --- | --- | --- | --- |
| Phase 1 | 解析 IPC 配置、校验参数、决定 runtime 分发 | `parseRunSimConfig(...)`、`isDynamicContactEnabled(...)` | `RunSimConfigParseTest` 的 parser 用例 |
| Phase 2 | external near-contact sample 检测、external barrier、external alpha 上界 | `TriangleMeshExternalContactHandler`、`PointPenetrationBarrierEnergy` | `contact_embedding_test`、`pointPenetrationBarrierEnergy_test` |
| Phase 3 | self near-contact active pair、self barrier、self alpha 上界 | `TriangleMeshSelfContactDetection`、`TriangleMeshSelfContactHandler`、`PointTrianglePairBarrierEnergy` | `self_contact_handler_test` |
| Phase 4 | 把 external/self barrier 按帧挂进 dynamic incremental potential | `runSimCore.cpp`、`TimeIntegrator::addGeneralImplicitForceModel(...)` | API smoke |
| Phase 5 | 把 feasible alpha callback 透传到 Newton | `setAlphaTestFunc(...)`、`TimeIntegratorSolver`、`NewtonRaphsonSolver` | `contact_embedding_test`、`self_contact_handler_test`、API smoke |
| Phase 6 | 把数学单测、handler 行为回归、runtime smoke 和 example README 组织成验证闭环 | `tests/*`、`examples/pulled-cubic-box-self-ipc` | `ctest` / focused smoke |

这张表有两个关键信息：

1. 每个 phase 在当前仓库里都有明确代码入口。
2. 每个 phase 也都有明确测试或 example 作为外部证据。

## 5. 当前 IPC 主线的数据流

如果把当前 `ipc-barrier` runtime 压成一条执行链，可以按下面顺序理解：

1. `parseRunSimConfig(...)` 读取 `contact-model`、external/self `dhat`、external/self `kappa`、`ipc-alpha-safety` 和 `ipc-enable-feasible-line-search`。
2. `runSimFromConfig(...)` 根据 `isDynamicContactEnabled(...)` 决定是否需要创建 external / self contact handler。
3. dynamic 帧循环开始时，先执行 `clearGeneralImplicitForceModel()` 和 `clearAlphaTestFunc()`，清掉上一帧 contact runtime。
4. external path 调用 `TriangleMeshExternalContactHandler::execute(..., externalIpcDhat)`，把当前 surface sample 映到最近外部面，产生 active samples。
5. self path 调用 `TriangleMeshSelfContactHandler::execute(..., selfIpcDhat)`，通过 BVH query band、sample seed refine 和局部搜索产生 active point-triangle pairs。
6. external active set 非空时，runtime 构造 `PointPenetrationBarrierEnergy`，配置 `setComputePosFunction(...)`、`setBuffer(...)`、`setCoeff(externalIpcKappa)`，并通过 `addGeneralImplicitForceModel(...)` 挂进 timestep。
7. self active set 非空时，runtime 构造 `PointTrianglePairBarrierEnergy`，配置 `setToPosFunction(...)`、`setBuffer(...)`、`setCoeff(selfIpcKappa)`，再通过同一条 general model 接缝挂进 timestep。
8. 如果启用了 feasible line search，runtime 会安装一个 merged callback，在当前 Newton 状态下分别查询 external / self 的 alpha upper bound，再取最小值。
9. `ImplicitBackwardEulerTimeIntegrator` 在 `assembleImplicitModels()` 之后，把 barrier 当成 `IMT_GENERAL` 模型一起参与 `updateA()`、`updateb()` 和一步求解。
10. `NewtonRaphsonSolver` 在 line search 前先消费 `alphaTestFunc` 返回的上界；若上界为零，直接返回 `SR_FEASIBLE_STEP_ZERO`。

这个链条说明当前 repo 的 barrier 不是 solver 外挂，也不是 post-process，而是每帧一步总势能的一部分。

## 6. 关键代码入口

当前 walkthrough 会反复回到下面这些文件：

- `src/api/runSimCore.h`
- `src/api/runSimCore.cpp`
- `src/core/solve/contact/triangleMeshExternalContactHandler.h`
- `src/core/solve/contact/triangleMeshExternalContactHandler.cpp`
- `src/core/solve/contact/triangleMeshSelfContactDetection.h`
- `src/core/solve/contact/triangleMeshSelfContactDetection.cpp`
- `src/core/solve/contact/triangleMeshSelfContactHandler.h`
- `src/core/solve/contact/triangleMeshSelfContactHandler.cpp`
- `src/core/solve/contact/pointPenetrationBarrierEnergy.h`
- `src/core/solve/contact/pointPenetrationBarrierEnergy.cpp`
- `src/core/solve/contact/pointTrianglePairBarrierEnergy.h`
- `src/core/solve/contact/pointTrianglePairBarrierEnergy.cpp`
- `src/core/solve/simulation/timeIntegrator.h`
- `src/core/solve/simulation/timeIntegrator.cpp`
- `src/core/solve/simulation/timeIntegratorSolver.cpp`
- `src/core/solve/simulation/implicitBackwardEulerTimeIntegrator.cpp`
- `src/core/solve/nonlinearOptimization/NewtonRaphsonSolver.h`
- `src/core/solve/nonlinearOptimization/NewtonRaphsonSolver.cpp`

这组文件已经覆盖了配置、检测、能量、装配、solver callback 和验证主入口。

## 7. 当前 walkthrough 依赖的验证证据

这组 walkthrough 对应的验证入口已经全部在仓库中：

- `tests/core/energy/pointPenetrationBarrierEnergy_test.cpp`
- `tests/core/energy/pointTrianglePairBarrierEnergy_test.cpp`
- `tests/core/scene/contact_embedding_test.cpp`
- `tests/core/scene/self_contact_handler_test.cpp`
- `tests/api/runSimConfig_parse_test.cpp`
- `examples/cubic-box-ipc/cubic-box-ipc-ls.json`
- `examples/pulled-cubic-box-self-ipc/pulled-cubic-box-self-ipc.json`
- `examples/pulled-cubic-box-self-ipc/README.md`

因此下面六个 phase 文档描述的都不是“计划中的做法”，而是已经有代码、测试和 example 支撑的当前行为。

## 8. 当前建模边界

同时也要把当前系统的边界看清楚：

- walkthrough 只覆盖 dynamic path，不展开 static path
- IPC 模式下 friction 被显式禁止
- self path 仍然是 sample-based point-triangle active set
- feasible alpha 是 contact-aware upper bound，不包含 inversion-free filter
- walkthrough 跟的是当前 `SO_NEWTON` 主路径
- current self/external geometry 仍然优先复用 repo 现有 sampled surface 路线，而不是 full primitive PT/EE IPC library

因此这组文档描述的不是 paper-equivalent full IPC，而是当前 repo-aligned、dynamic-path-first 的 IPC runtime。

## 9. 建议阅读顺序

- [Phase 1: Config and Dispatch](phase-1-config-and-dispatch.md)
  先看 `ipc-barrier` 怎样作为 contact model 进入 runtime。
- [Phase 2: External Barrier Energy](phase-2-external-barrier-energy.md)
  再看 external near-contact sample 怎样变成 barrier energy 和 external feasible alpha。
- [Phase 3: Self Near-Contact Active Set](phase-3-self-near-contact-active-set.md)
  先把 self path 怎样从 BVH query band 细化到 sampled point-triangle active pair 讲清楚。
- [Phase 4: Dynamic Incremental Potential](phase-4-dynamic-incremental-potential.md)
  然后看 external/self 两条 barrier 怎样真正进入 timestep 的 general `PotentialEnergy` 组装。
- [Phase 5: Feasible Line Search](phase-5-feasible-line-search.md)
  再看 merged alpha callback 怎样接到 Newton。
- [Phase 6: Tests and Validation](phase-6-tests-and-validation.md)
  收尾时看哪些测试、smoke 和 example 已经把这条链路钉住。
