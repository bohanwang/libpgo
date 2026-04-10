# Phase 6：测试与验证

## 1. 这阶段要回答什么

这一阶段专门回答下面四个问题：

1. 当前 repo 里哪些测试真正证明了 IPC runtime 已经接通。
2. external / self / merged 三条路径分别由哪些回归固定下来。
3. 文档里说的 “sample-based self IPC” 和 “merged feasible alpha callback” 有什么直接证据。
4. 当前验证闭环已经覆盖到哪里，哪些边界仍然是已知保留项。

IPC walkthrough 到这里不再讨论“应该怎样实现”，而是回答：

> 当前仓库里，哪些东西已经被自动化验证写死了。

## 2. 当前 repo 的验证闭环已经是什么形状

当前仓库对 IPC walkthrough 的验证不再依赖“看一下动画是否像样”，而是已经形成：

> 数学单测 -> handler / active-set 行为回归 -> API smoke -> example README

这条分层闭环。

换句话说，文档里每个重要断言都能在仓库中找到对应证据。

## 3. 关键验证入口

当前 walkthrough 对应的核心验证入口是：

- `tests/core/energy/pointPenetrationBarrierEnergy_test.cpp`
- `tests/core/energy/pointTrianglePairBarrierEnergy_test.cpp`
- `tests/core/scene/contact_embedding_test.cpp`
- `tests/core/scene/self_contact_handler_test.cpp`
- `tests/api/runSimConfig_parse_test.cpp`
- `examples/pulled-cubic-box-self-ipc/README.md`

这几处分别覆盖：

- barrier 数学
- handler 行为
- runtime 接线
- example 侧的当前能力范围说明

## 4. Workstream A：barrier energy 数学回归

### 4.1 external barrier 的数学测试

external barrier 的数学测试在：

- `tests/core/energy/pointPenetrationBarrierEnergy_test.cpp`

当前覆盖的核心语义是：

- `InactiveRegionReturnsZero`
- `ActiveRegionMatchesFormulaAndFiniteDifference`
- `ClampAvoidsNaNForNonPositiveDistances`

也就是说，external barrier 的：

- 支撑区间
- 一阶导数
- 二阶导数
- 非正距离附近的数值稳定性

都已经有 committed test 固定。

### 4.2 self barrier 的数学测试

self barrier 的数学测试在：

- `tests/core/energy/pointTrianglePairBarrierEnergy_test.cpp`

当前覆盖的核心语义同样是：

- `InactiveRegionReturnsZero`
- `ActiveRegionMatchesFormulaAndFiniteDifference`
- `ClampAvoidsNaNForNonPositiveDistances`

这说明 self barrier 虽然几何上仍是 sample-based，但它作为 `PotentialEnergy` 的数值行为已经被单测锁住。

## 5. Workstream B：handler 与 active-set 行为回归

### 5.1 external handler 的几何与 alpha 上界回归

external handler 的几何与 alpha 上界回归在：

- `tests/core/scene/contact_embedding_test.cpp`

这里不仅检查 sample embedding 与 barycentric embedding 的一致性，还直接覆盖：

- `ExternalBarrierActivationIncludesNearContactSamples`
- `FeasibleStepUpperBoundMatchesAnalyticPlaneBound`
- `FeasibleStepUpperBoundScalesWithAlphaSafety`
- `FeasibleStepUpperBoundIsOneForMotionAwayFromContact`
- `FeasibleStepUpperBoundReturnsOneWithNoActiveSamples`
- `FeasibleStepUpperBoundReturnsZeroWhenAlreadyInsideSafeMargin`
- `FeasibleStepUpperBoundAllowsRecoveryWhenInsideSafeMarginButMovingAway`

这意味着 external IPC 的两件核心事情都已经有直接证据：

- near-contact sample 确实会激活
- feasible alpha upper bound 的 safe-band 语义确实如文档所述

### 5.2 self handler 的几何与 alpha 回归

self handler 的几何与 alpha 回归在：

- `tests/core/scene/self_contact_handler_test.cpp`

它当前已经覆盖：

- near-contact activation
- activation threshold suppression
- zero-band fallback
- legacy penalty path compatibility
- closest target selection
- barrier finite energy
- inward / outward / safe-band / no-active-pair 的 alpha 行为

因此 walkthrough 中关于 self path 的下面这些描述，都已经有 committed handler 级测试支撑：

- sample seed refine
- one target per sample
- zero-band fallback
- sample-based feasible alpha

## 6. Workstream C：API smoke

`tests/api/runSimConfig_parse_test.cpp` 把 `runSimCore` 层的 IPC runtime 固定成多组 smoke：

1. `RunSimFromConfigCubicDynamicIpcSmokeTest`
   基础 IPC dynamic smoke
2. `RunSimFromConfigCubicDynamicIpcNearContactActivatesExternalBarrier`
   external-only IPC near-contact smoke
3. `RunSimFromConfigCubicDynamicIpcSelfBarrierSmokeTest`
   self-only IPC smoke
4. `RunSimFromConfigCubicDynamicIpcMergedBarrierSmokeTest`
   ext+self merged smoke
5. `RunSimFromConfigCubicDynamicIpcMergedDeterministicSmokeTest`
   merged deterministic smoke

这些 smoke 不只是检查程序返回值，还直接检查日志里的关键行为：

- `# external active samples:`
- `# self active pairs:`
- `IPC feasible alpha callback active.`
- `IPC feasible alpha merged callback active.`
- `IPC feasible alpha upper bound:`

并且会显式断言旧 warning 已经消失：

- `IPC barrier energy not yet implemented for external contact.`
- `IPC barrier energy not yet implemented for self contact.`

这意味着当前 runtime walkthrough 的 “已接通” 不是主观判断，而是 API 级断言。

## 7. Workstream D：example 与 README

当前 walkthrough 对应的主 example 是：

- `examples/pulled-cubic-box-self-ipc/pulled-cubic-box-self-ipc.json`
- `examples/pulled-cubic-box-self-ipc/README.md`

这个 example 在当前 repo 里的职责非常明确：

- 回归 self IPC runtime
- 回归 merged ext+self 场景
- 给 sample-based self path 提供一个可重复的 integration scene

README 里已经显式写清当前实现范围：

- external and self `ipc-barrier` runtime branches
- self near-contact active-set construction
- self barrier energy
- external / self / merged feasible-alpha runtime filtering
- sample-based self path

因此 example 文档本身也是 walkthrough 的一部分证据，而不只是一个演示目录。

## 8. 当前验证口径已经能证明什么

把这些测试和 example 放在一起看，当前 repo 已经能够验证下面几件事：

1. external barrier 的数学与数值稳定性正确
2. self barrier 的数学与数值稳定性正确
3. external near-contact active sample 会被激活
4. self near-contact active pair 会被激活
5. external / self feasible alpha upper bound 都会生效
6. merged callback 确实在 runtime 中工作
7. merged IPC 路径在 deterministic mode 下可重复

这正是 walkthrough 中每个 phase 都能直接落到代码和测试上的原因。

## 9. 当前保留边界

同时，Phase 6 也把当前验证边界固定了下来：

- self 仍然是 sample-based，不是 primitive PT/EE guarantee
- walkthrough 不把当前实现宣称成 full primitive IPC benchmark
- 复杂场景下 solver 仍可能出现软终止状态，例如 `max_iter_reached`
- 当前验证不覆盖 friction、adaptive `kappa`、primitive-level 演化和 inversion-free filter

因此这一阶段的作用不是把系统“说成完整”，而是把：

- 当前已经成立的能力
- 当前仍然保留的边界

同时写死。

## 10. 建议的验证命令

如果要本地复核 walkthrough 相关事实，当前最直接的入口是：

```bash
cd libpgo/build/core-debug/bin
./core_scene_self_contact_handler_test
./api_runSim_config_parse_test --gtest_filter='RunSimConfigParseTest.RunSimFromConfigCubicDynamicIpcNearContactActivatesExternalBarrier:RunSimConfigParseTest.RunSimFromConfigCubicDynamicIpcSelfBarrierSmokeTest:RunSimConfigParseTest.RunSimFromConfigCubicDynamicIpcMergedBarrierSmokeTest:RunSimConfigParseTest.RunSimFromConfigCubicDynamicIpcMergedDeterministicSmokeTest'
```

以及更偏数学的：

```bash
ctest --test-dir ../ -R "core\\.energy\\.(point_penetration_barrier|point_triangle_pair_barrier)|core\\.scene\\.(contact_embedding|self_contact_handler)|api\\.runSim_config\\.parse" --output-on-failure
```

这两组命令基本覆盖了 walkthrough 中最关键的数学、handler 和 runtime 断言。

上一阶段： [Phase 5](phase-5-feasible-line-search.md)  
返回总览： [IPC Overview](index.md)
