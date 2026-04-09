# Phase 2：external contact 的 barrier energy

## 1. 这阶段要回答什么

这一阶段专门回答下面四个问题：

1. external IPC 为什么在当前 repo 里不是“穿透以后再补一个力”，而是“近接触 sample 提前激活”。
2. `TriangleMeshExternalContactHandler` 当前到底缓存了哪些 active sample 数据。
3. `PointPenetrationBarrierEnergy` 怎样把这批 sample 数据变成能量、梯度和 Hessian。
4. external feasible alpha upper bound 现在建立在什么几何近似上。

如果只看到 `buildBarrierEnergy(...)` 这个接口，会误以为这是一个纯数学 energy 类。  
但当前 external IPC 真正完成的是一条完整 vertical slice：

> sampled surface point -> near-contact active sample -> fixed plane snapshot -> barrier energy -> feasible alpha upper bound

## 2. 这阶段在整条 IPC 主线里解决什么

Phase 1 只是把 `ipc-barrier` 作为 contact model 接进 runtime。  
真正让 external IPC “开始工作” 的，是当前这一阶段：

- external handler 不再只关心 penetration
- active set 的判据从 `d <= 0` 扩展到 `d < dhat`
- barrier energy 不再依赖 penalty energy 的二次型
- solver 侧可行步长过滤开始有 external 几何依据

所以 Phase 2 的本质不是“再讲一遍接触检测”，而是：

> 把 external contact 从 penetration-only 的 penalty runtime，推进成 near-contact barrier runtime。

## 3. 当前 repo 已经完成了什么

当前仓库里，external 路径已经具备下面这些落地能力：

- `TriangleMeshExternalContactHandler::execute(..., activationDistance)` 已经存在
- active set 的筛选语义已经是“保留所有 `signedDistance < activationDistance` 的 sample”
- handler 已经缓存 `constraintNormals`、`constraintTargetPositions`、`constraintSignedDistances`
- handler 已经保存 `barycentricIdx`、`barycentricWeights` 和 `contactedSamples`
- `buildBarrierEnergy(...)` 已经能构造 `PointPenetrationBarrierEnergy`
- external handler 已经能直接给出：
  - `computeSurfaceAlphaUpperBound(...)`
  - `computeEmbeddedAlphaUpperBound(...)`

这说明 external IPC 在当前 repo 里已经不只是“检测近接触”，而是完整的 runtime 接口层。

## 4. `TriangleMeshExternalContactHandler` 当前是什么形状

### 4.1 external path 仍然建立在 sampled surface 上

handler 构造阶段已经完成：

- surface triangle sampling
- sample 去重
- `sampleInfoAndIDs`
- `sampleTriangleIDs`
- sample 到真实 DOF 的插值矩阵 `interpolationMatrix`
- `sampleWeights`

因此 runtime 每帧真正变化的是：

- surface 顶点位置
- sample 当前位置
- 最近的外部目标快照

而不是 sample 集合本身。

这也解释了当前 repo 的建模选择：

> external 接触基元是 sampled surface point，而不是 primitive-level point/edge/face 对。

### 4.2 handler 还维护了 external target 侧的加速结构

除了 sample 数据，构造函数还会建立：

- `surfaceMeshRuntime`
- `surfaceMeshNormals`
- `externalSurfaceBVTrees`
- `externalSurfaceNormals`

这意味着 external handler 当前不仅持有源 surface 的采样数据，也持有外部障碍物侧的：

- BVH
- pseudo normal
- runtime mesh snapshot

因此 `execute(...)` 时可以直接做最近三角形查询，而不需要重新初始化外部几何结构。

## 5. active sample 当前怎样被生成

### 5.1 runtime 先把 surface 位移映回 sample 位移

当前 external path 在 dynamic loop 中会调用：

```cpp
externalContactHandler->execute(usurf.data(), externalIpcDhat);
```

handler 内部的顺序是：

1. 用 `restP + u` 得到当前 surface 顶点位置 `curP`
2. 更新 `surfaceMeshRuntime`
3. 用 `computeSamplePosition(...)` 把顶点位置映到 `sampleCurP`
4. 再进入真正的 external 查询

所以 external near-contact 检测发生的真正空间是：

- 当前 sampled surface

而不是原始 surface 顶点数组本身。

### 5.2 每个 sample 只保留一个 best external candidate

当前 `execute(double activationDistance)` 对每个 sample 的逻辑可以概括成：

1. 遍历所有 external surface
2. 对每个 object 做 `closestTriangleQuery(...)`
3. 根据 closest triangle feature 取 pseudo normal
4. 做法向朝向过滤：`tgtNormal.dot(srcNormal) > 0` 的候选被排除
5. 计算有符号距离 `signedDistance = (srcPos - tgtPos) · tgtNormal`
6. 只保留 `signedDistance < activationDistance` 的候选
7. 在这些候选中选 `signedDistance` 最小的那个

这里的选择标准很明确：

- 已穿透时，选更负的 signed distance
- 近接触但未穿透时，选更小的正距离

因此当前 external active set 的单位不是 “一个 sample 对多个候选面”，而是：

```text
一个 sample -> 一个最近 external target snapshot
```

### 5.3 `# external active samples` 日志就是当前 active set 规模

`execute(...)` 完成后，handler 会打印：

```text
# external active samples: N
```

这条日志并不是调试残留，而是 API smoke 会直接依赖的可观测行为。  
当前 repo 把 external active set 是否真正建立起来，显式暴露成了 runtime 日志。

## 6. 当前 active sample 快照里到底保存了什么

当一个 sample 被保留下来后，handler 会把下面这些量写入 frame-local 缓存：

- `constraintCoeffs[count] = sampleWeights[info.sId]`
- `constraintNormals`
- `constraintTargetPositions`
- `constraintSignedDistances`
- `contactedSamples`
- `contactedTriangles`
- `barycentricIdx`
- `barycentricWeights`

这组缓存有两个重要含义：

### 6.1 barrier energy 不需要重新做最近点查询

`PointPenetrationBarrierEnergy` 直接消费：

- 当前 sample 与真实 DOF 的插值关系
- 当前 sample 对应的 target point
- 当前 sample 对应的 pseudo normal

因此 energy 层不需要再去 external BVH 查询 closest point。

### 6.2 feasible alpha 也复用同一份快照

`computeSurfaceAlphaUpperBound(...)` / `computeEmbeddedAlphaUpperBound(...)` 读取的也是这组 active sample 快照。  
这意味着当前 repo 的 external barrier 与 external alpha upper bound 使用的是同一帧检测结果，而不是两套彼此独立的几何查询。

## 7. `PointPenetrationBarrierEnergy` 当前怎样工作

### 7.1 energy 的局部距离模型是什么

对每个 active sample，当前 external barrier 使用的局部距离就是：

```text
d = (p - p0) · n
```

其中：

- `p`
  由 `setComputePosFunction(...)` 和 `barycentricWeights` 从当前 DOF 恢复出来
- `p0`
  来自 handler 在当前帧缓存的 `constraintTargetPositions`
- `n`
  来自 handler 在当前帧缓存的 `constraintNormals`

因此 external barrier 当前是典型的 sampled point vs frozen target plane snapshot 语义。

### 7.2 energy 的支撑区间已经是 barrier 语义

`PointPenetrationBarrierEnergy` 当前行为明确分成两段：

- `distance >= dhat`
  返回零能量、零梯度、零 Hessian
- `distance < dhat`
  使用 `BarrierFunctions::logBarrierEnergy / Gradient / Hessian`

这让 external barrier 在当前 repo 里成为真正的支撑区间 barrier，而不是简单把 penalty 系数调大。

### 7.3 数值安全域通过 `normalizePositiveZero(...)` 固定下来

external barrier 当前提供：

```cpp
PointPenetrationBarrierEnergy::normalizePositiveZero(dhat, positiveZero)
```

其实现会把 `positiveZero` 约束在一个安全上界以内，并在未显式给值时使用相对 `dhat` 的推荐值。  
这一步的作用不是改几何语义，而是避免 `logBarrier*` 在非正距离附近出现 NaN 或不稳定值。

### 7.4 Hessian 模板仍然沿用 sample-to-DOF 插值关系

energy 构造时会根据每个 active sample 的 `barycentricIdx` 建立 Hessian 模板。  
在真正求 Hessian 时，再把：

- `n n^T`
- barrier 二阶导
- sample 权重
- barycentric 权重

共同装配回真实 DOF。

因此当前 external barrier 的数学结构已经和 runtime 的 sample embedding 严格对齐。

## 8. external feasible alpha upper bound 当前怎样计算

### 8.1 两个 public 入口只是输入空间不同

当前 handler 提供：

- `computeSurfaceAlphaUpperBound(...)`
- `computeEmbeddedAlphaUpperBound(...)`

两者最终都会转到：

- `computeSampleAlphaUpperBound(...)`

差别只在于：

- 一个输入的是 surface-sized displacement
- 一个输入的是 full embedded displacement

但最终都会映到 sample 空间再计算。

### 8.2 upper bound 使用的是同一份 active plane snapshot

对每个 active sample，当前上界计算会使用：

- `sampleID`
- `constraintNormals`
- `constraintTargetPositions`
- `sampleRestP`
- 当前 sample displacement
- 当前 sample 方向 `sampleDu`

然后计算：

- 当前距离 `d0`
- 法向闭合速度 `dDot`

并基于 `dSafe` 构造 conservative `alphaUpper`。

### 8.3 safe band 语义当前已经被显式编码

当前实现有两条关键分支：

- `d0 <= dSafe`
  如果方向继续 inward，则直接返回 `0`
  如果方向 outward，则允许继续恢复
- `dDot < 0`
  才会真正收缩 `alphaUpper`

这意味着当前 external feasible alpha 不是“只要在 safe band 内就冻结”，而是：

> 阻止继续 inward，但允许 outward recovery

### 8.4 `ipc-alpha-safety` 当前只作用于 upper bound

当前 external handler 会在最终结果上乘：

- `std::clamp(alphaSafety, 0.0, 1.0)`

因此 `ipc-alpha-safety` 当前影响的是 solver-side feasible step margin，而不是 barrier energy 本身。

## 9. 当前 repo 的验证证据

`tests/core/energy/pointPenetrationBarrierEnergy_test.cpp` 已经固定了 external barrier 的数学语义：

- `InactiveRegionReturnsZero`
- `ActiveRegionMatchesFormulaAndFiniteDifference`
- `ClampAvoidsNaNForNonPositiveDistances`

`tests/core/scene/contact_embedding_test.cpp` 则固定了 handler 和 alpha 上界的几何语义：

- `ExternalBarrierActivationIncludesNearContactSamples`
- `FeasibleStepUpperBoundMatchesAnalyticPlaneBound`
- `FeasibleStepUpperBoundScalesWithAlphaSafety`
- `FeasibleStepUpperBoundIsOneForMotionAwayFromContact`
- `FeasibleStepUpperBoundReturnsOneWithNoActiveSamples`
- `FeasibleStepUpperBoundReturnsZeroWhenAlreadyInsideSafeMargin`
- `FeasibleStepUpperBoundAllowsRecoveryWhenInsideSafeMarginButMovingAway`

再往上，`tests/api/runSimConfig_parse_test.cpp` 中的：

- `RunSimFromConfigCubicDynamicIpcNearContactActivatesExternalBarrier`

会直接验证 runtime 日志里存在：

- `# external active samples: ...`
- `IPC feasible alpha callback active.`

因此 Phase 2 的 external path 已经有数学、handler 和 API 三层证据。

## 10. 这一阶段不展开什么

这一阶段还不展开下面这些内容：

- self near-contact active set
- barrier 怎样进入 timestep 的总势能
- merged feasible alpha callback 的 runtime 组装
- static path
- primitive-level exact IPC

当前 external IPC 仍然是 repo-aligned 的 sample-based 路线，而且 barrier energy 内部不会重新做 closest-point 查询。

上一阶段： [Phase 1](phase-1-config-and-dispatch.md)  
下一阶段： [Phase 3](phase-3-dynamic-incremental-potential.md)
