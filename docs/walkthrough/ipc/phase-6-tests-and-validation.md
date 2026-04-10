# Phase 6：测试与验证

验证闭环分四层：数学单测 → handler 行为回归 → API smoke → example。文档中每个重要断言都能在仓库中找到对应证据。

## A. Barrier Energy 数学回归

### `pointPenetrationBarrierEnergy_test.cpp`（external）

| 测试 | 验证内容 |
| --- | --- |
| `InactiveRegionReturnsZero` | $d \geq \hat{d}$ 时 energy/gradient/Hessian 为零 |
| `ActiveRegionMatchesFormulaAndFiniteDifference` | active region 的一阶/二阶导数与有限差分一致 |
| `ClampAvoidsNaNForNonPositiveDistances` | $d \leq 0$ 时不产生 NaN |

### `pointTrianglePairBarrierEnergy_test.cpp`（self）

同上三项。Self barrier 虽然几何上是 sample-based，但 `PotentialEnergy` 层的数值行为已被锁住。

## B. Handler 行为回归

### `contact_embedding_test.cpp`（external handler）

| 测试 | 验证内容 |
| --- | --- |
| `ExternalBarrierActivationIncludesNearContactSamples` | near-contact sample 被正确激活 |
| `FeasibleStepUpperBoundMatchesAnalyticPlaneBound` | alpha 上界匹配解析 plane bound |
| `FeasibleStepUpperBoundScalesWithAlphaSafety` | alpha 随 safety 线性缩放 |
| `FeasibleStepUpperBoundIsOneForMotionAwayFromContact` | outward 方向不收缩 |
| `FeasibleStepUpperBoundReturnsOneWithNoActiveSamples` | 无 active sample 时返回 1 |
| `FeasibleStepUpperBoundReturnsZeroWhenAlreadyInsideSafeMargin` | safe-band 内 inward 返回 0 |
| `FeasibleStepUpperBoundAllowsRecoveryWhenInsideSafeMarginButMovingAway` | safe-band 内 outward 允许 recovery |

### `self_contact_handler_test.cpp`（self handler）

| 测试 | 验证内容 |
| --- | --- |
| `NearContactActivationUsesUnsignedDistance` | 使用 unsigned distance 作为激活判据 |
| `ActivationDistanceThresholdSuppressesPairs` | 缩小 band 时 active pair 被压掉 |
| `ZeroBandFallsBackToCollisionOnlyWhenNotColliding` | 零带宽回退到旧碰撞路径 |
| `ZeroBandFallbackHandlesCollidingPairs` | 零带宽 fallback 正确处理碰撞 |
| `ClosestTargetSelectionKeepsOneTargetPerSample` | 每个 sample 只保留一个 target |
| `LegacyPenaltyPathStillWorks` | 旧 penalty 路径兼容 |
| `BarrierBuilderProducesFiniteEnergyForNearContactPairs` | barrier 产生有限正能量 |
| `AlphaUpperBoundShrinksForInwardDirections` | inward 方向收缩 alpha |
| `AlphaUpperBoundScalesWithAlphaSafety` | safety scaling |
| `AlphaUpperBoundAllowsOutwardRecovery` | outward recovery |
| `AlphaUpperBoundBlocksFurtherInwardMotionInsideSafeBand` | safe-band 内阻止 inward |
| `AlphaUpperBoundReturnsOneWithNoActivePairs` | 无 active pair 时返回 1 |

## C. API Smoke

`runSimConfig_parse_test.cpp` 中的 runtime smoke：

| 测试 | 覆盖路径 |
| --- | --- |
| `RunSimFromConfigCubicDynamicIpcSmokeTest` | 基础 IPC dynamic |
| `...IpcNearContactActivatesExternalBarrier` | external-only near-contact |
| `...IpcSelfBarrierSmokeTest` | self-only |
| `...IpcMergedBarrierSmokeTest` | ext+self merged |
| `...IpcMergedDeterministicSmokeTest` | merged deterministic |

这些 smoke 检查日志中的关键行为：

- `# external active samples: N`
- `# self active pairs: N`
- `IPC feasible alpha callback active.` / `...merged callback active.`
- `IPC feasible alpha upper bound: ...`
- 旧 warning（`...not yet implemented...`）已消失

## D. Example

`examples/pulled-cubic-box-self-ipc/`：ext+self merged IPC 场景，附 README 说明实现范围。

## 本地验证命令

```bash
# 聚焦 IPC 相关测试
ctest --test-dir libpgo/build/core-debug \
  -R "core\\.energy\\.(point_penetration_barrier|point_triangle_pair_barrier)|core\\.scene\\.(contact_embedding|self_contact_handler)|api\\.runSim_config\\.parse" \
  --output-on-failure
```

## 保留边界

- Self 是 sample-based，不是 primitive PT/EE
- 不覆盖 friction、adaptive $\kappa$、inversion-free filter
- 复杂场景下 solver 可能出现 `max_iter_reached`
- 不等同于 paper-equivalent full IPC

---

上一阶段：[Phase 5](phase-5-feasible-line-search.md)
返回总览：[IPC Overview](index.md)
