# Phase 5：self-contact 的 near-contact active set

## 1. 这阶段要回答什么

这一阶段专门回答下面四个问题：

1. self-contact 在当前 repo 里怎样从 triangle-pair 候选细化成最终 active pair。
2. 为什么 current self IPC 仍然是 sample-based，而不是 primitive PT/EE。
3. `TriangleMeshSelfContactDetection` 和 `TriangleMeshSelfContactHandler` 当前分别负责什么。
4. self barrier energy 和 self feasible alpha upper bound 共享的是哪一批缓存。

如果只说“self 也有 barrier 了”，会掩盖当前 self 路径真正复杂的地方。  
Phase 5 的核心不在 barrier 标量公式，而在：

> broad phase 候选、sample seed refine、局部搜索和最终 point-triangle active pair 怎样被组织成一条稳定 runtime 链。

## 2. 这阶段在整条 IPC 主线里解决什么

external 路径相对直接，因为目标外部面是独立的 external object。  
self 路径更难的地方在于：

- broad phase 不能只看“已经相交”
- active set 不能直接等于 triangle pair
- barrier 和 feasible alpha 都需要同一份稳定的 point-triangle snapshot

当前仓库已经把这一层整理成三段：

1. `TriangleMeshSelfContactDetection`
   先给出进入 `dhat` 带宽的 triangle-pair 候选
2. `TriangleMeshSelfContactHandler`
   把 triangle-pair 候选细化成 sampled point-triangle active pair
3. `PointTrianglePairBarrierEnergy`
   与 `computeEmbeddedAlphaUpperBound(...)`
   共同消费这份 active pair 快照

所以 Phase 5 实际上是 self IPC 的 geometry-and-runtime composition phase。

## 3. 当前 repo 已经完成了什么

当前仓库里，self 路径已经具备下面这些能力：

- `TriangleMeshSelfContactDetection` 已支持带 `activationDistance` 的查询接口
- self handler 已支持：
  - `execute(u0, activationDistance)`
  - 零带宽 fallback 到旧碰撞路径
- near-contact 路径已形成：
  - candidate triangle pair
  - sample seed refine
  - final active pair
  这三层结构
- self handler 当前已缓存：
  - `contactedTrianglePairs`
  - `contactedTriangleIDs`
  - `activeClosestPoints`
  - `activeNormals`
  - `activeDistances`
  - `activeBarycentricWeights`
- `buildBarrierEnergy(...)` 已能构造 `PointTrianglePairBarrierEnergy`
- `computeEmbeddedAlphaUpperBound(...)` 已能基于同一份 active pair 快照给出 self feasible alpha

也就是说，当前 self IPC 已经形成了一条完整的 sample-based vertical slice：

> BVH query band -> sampled point-triangle active pair -> self barrier energy -> self feasible alpha

## 4. 当前 self path 的基础形状是什么

### 4.1 self 几何单位仍然是 sampled surface point

`TriangleMeshSelfContactHandler` 构造时已经完成：

- surface triangle sampling
- sample 去重
- `sampleInfoAndIDs`
- `sampleTriangleIDs`
- `vertexID2SampleIDsLinear`
- `interpolationMatrix`
- `sampleWeights`

因此当前 self IPC 的基础几何单位并不是原始 mesh 顶点，也不是 primitive-level PT/EE，而是 sampled surface point。

### 4.2 detection 与 handler 已经明确拆层

当前 self 路径不是一个类做完所有事情，而是明确分成两层：

- `TriangleMeshSelfContactDetection`
  负责 BVH traversal、AABB query band、candidate triangle pair
- `TriangleMeshSelfContactHandler`
  负责 sample seed refine、局部目标选择、barrier 构造和 alpha upper bound

这让 broad phase 与 active-set 细化职责保持清晰分开。

## 5. broad phase 当前怎样生成 candidate triangle pairs

### 5.1 detection 已支持 near-contact 查询模式

当前 `TriangleMeshSelfContactDetection` 支持：

- `execute(const double* positions0, double activationDistance)`
- `execute(const double* positions0, const double* positions1, double activationDistance)`

当 `activationDistance > 0` 时，它不再只做 collision-only 判据，而是改用 query band 语义。

### 5.2 node 层和 triangle 层都使用 query band

当前 detection 内部会分别使用：

- `nodesWithinQueryBand(...)`
- `trianglesWithinQueryBand(...)`

其中 band 取：

- `std::max(activationDistance, 0.0)`

因此当 `activationDistance > 0` 时，broad phase 当前回答的问题不再是：

```text
哪些 pair 已经相交
```

而是：

```text
哪些 pair 可能进入了 selfIpcDhat 带宽
```

### 5.3 DCD / CCD 与 query band 共用同一条接口

当前 detection 的内部逻辑分成四种情况：

- DCD + `band <= 0`
  继续使用传统 `DCDIntersect`
- CCD + `band <= 0`
  使用 `CCDIntersect`
- DCD + `band > 0`
  用 `AABBWithinBand(lastAABBs[a], lastAABBs[b], band)`
- CCD + `band > 0`
  用 `SweptAABB(...)` 再做 `AABBWithinBand(...)`

这说明当前 self near-contact broad phase 不是另起一套系统，而是在现有 BVH traversal 上增加了 query-band 语义。

### 5.4 candidate pair 还会经过共享顶点过滤和去重

在 leaf 层 `assignTriangles(...)` 中，当前 detection 仍然保留：

- `triIdA != triIdB`
- `!shareVertex(triIdA, triIdB)`

然后把结果写入 `candidateTrianglePairs`，最后统一：

- 规范 pair 顺序
- 排序
- 去重

因此 `getCandidateTrianglePairs()` 返回的已经是一份去重后的 broad-phase 候选集，而不是未整理的中间结果。

## 6. `execute(u0, activationDistance)` 当前怎样分叉

### 6.1 零带宽会回退到旧碰撞路径

当前 `TriangleMeshSelfContactHandler::execute(u0, activationDistance)` 在 `activationDistance <= 0.0` 时会：

1. 调 `executeDCD()`
2. 如果没有 `collidingTrianglePairs`
   直接清空 active pair
3. 如果有碰撞
   调 `handleContactDCD(0.0, 100)`

因此零带宽路径当前保留了 legacy penalty / collision 语义。

### 6.2 正带宽会进入 near-contact DCD

当 `activationDistance > 0.0` 时，handler 会进入：

```cpp
executeNearContactDCD(activationDistance, 100);
```

也就是说，当前 self IPC near-contact 路径并不是对旧 `handleContactDCD(...)` 的微调，而是有自己专门的入口。

## 7. `executeNearContactDCD(...)` 当前怎样做 sample seed refine

### 7.1 当前的三个中间缓冲

这一层最重要的三个中间数据是：

- `sampleSeedRecords`
- `sampleVisited`
- `activeSeedRecords`

其中：

- `sampleSeedRecords`
  记录每个 sample 当前最好的 seed triangle 与 `bestDist2`
- `sampleVisited`
  标记该 sample 是否已有 seed
- `activeSeedRecords`
  只保留最终进入下一轮细化的 sample seed

### 7.2 seed 的判据是 unsigned point-triangle distance

当前 `tryUpdateSeed(...)` 使用的不是 plane-signed depth，而是：

- `getSquaredDistanceToTriangle(sampleP, va, vb, vc)`

并要求：

- `dist2 < activationDistance^2`

这意味着当前 self near-contact active set 的 canonical 距离语义是：

> sampled point 到 target triangle 的 unsigned distance

这也是对应单测 `NearContactActivationUsesUnsignedDistance` 的直接来源。

### 7.3 每个 sample 只保留一个 best seed triangle

一旦某个 sample 在多个 triangle pair 候选里都落入 band，当前实现会只保留：

- `dist2` 最小的那个 seed triangle

因此这一层的结果已经不是 triangle-pair list，而是：

```text
sample -> best seed triangle
```

这一步对规模控制很重要，因为它避免了一对多 target 导致 active set 爆炸。

## 8. `finalizeActivePairsFromSeeds(...)` 当前怎样生成最终 active pair

### 8.1 active seed 之后还有一轮局部搜索

seed refine 完成后，handler 不会直接把 seed triangle 当最终 target。  
它还会进入：

```cpp
finalizeActivePairsFromSeeds(maxSearchingNumTriangles);
```

这一层会为每个 active seed 建立一个 `LocalSearchBuffer`，其中包含：

- `sampleID`
- `priority_queue`
- `searchFilter`

然后从 seed triangle 开始，沿 triangle 邻接关系做局部搜索。

### 8.2 局部搜索的目标是找到最终最近 triangle

当前每个候选 triangle 都会被评估：

- 最近点 `closestPt`
- barycentric weight `triBaryWeight`
- `dist2`

如果某个邻接 triangle 更近，就替换当前 best result。  
当队列中的候选已经明显远于当前 best result 时，搜索就停止。

这一步解释了为什么当前 repo 不直接把 seed triangle 当最终 active pair：  
seed 只是进入局部搜索的起点，真正的 active pair 仍然要经过最近目标选择。

### 8.3 当前还会过滤无效 sample/triangle 组合

在确定最终 closest triangle 之后，handler 还会检查：

- sample 是否本来就属于目标 triangle
- `excludedTriangles` 是否要求排除这组组合

只有通过这些检查的 pair 才会进入最终 active set。

### 8.4 当前最终缓存的不只是 pair ID

最终 active pair 写入的缓存包括：

- `contactedTrianglePairs`
- `contactedTriangleIDs`
- `activeClosestPoints`
- `activeNormals`
- `activeDistances`
- `activeBarycentricWeights`

其中 `contactedTrianglePairs` 当前是 4 个 sample id：

- `pair[0]`
  source sampled point
- `pair[1]`
  target triangle 第一个 corner sample
- `pair[2]`
  target triangle 第二个 corner sample
- `pair[3]`
  target triangle 第三个 corner sample

因此最终 active set 已经不是抽象 triangle id，而是直接可供 barrier 与 alpha upper bound 消费的 sampled point-triangle 连接关系。

### 8.5 normal 与 barycentric weight 也在这一层固定

当前实现中：

- `activeBarycentricWeights`
  来自 closest point 的 barycentric 坐标
- `activeNormals`
  优先取 `sample - closestPoint` 的单位方向
  若法向长度过小，则回退到 triangle 几何法向

这意味着 self barrier 与 self alpha 当前共享的是：

- 同一份 closest-point 几何快照
- 同一份 normal
- 同一份 barycentric weight

## 9. self barrier energy 当前怎样消费这份 active pair 快照

### 9.1 `buildBarrierEnergy(...)` 直接读取 handler 缓存

当前 runtime 中一旦 `getNumActivePairs() > 0`，就会执行：

```cpp
selfBarrierEnergy = selfCD->buildBarrierEnergy(selfIpcDhat);
```

`PointTrianglePairBarrierEnergy` 当前直接消费：

- `contactedTrianglePairs`
- `sampleEmbeddingWeights`
- `sampleWeights`
- `activeNormals`
- `activeBarycentricWeights`

因此 barrier energy 不会再自己去做 self-contact 几何搜索。

### 9.2 当前距离模型是 sampled linearization

对每个 active pair，当前 self barrier 使用的距离是：

```text
d = n · (p - (w0 t0 + w1 t1 + w2 t2))
```

其中：

- `p`
  source sampled point 的当前位置
- `t0/t1/t2`
  target triangle 三个 corner sample 的当前位置
- `n`
  来自 handler 的当前 frame active normal
- `w`
  来自 handler 的当前 frame barycentric weight

这一步把当前 self IPC 的建模边界说得很清楚：

- 已经是 near-contact barrier
- 但仍然是 sample-based point-triangle linearization
- 还不是 primitive-level exact IPC

## 10. self feasible alpha upper bound 当前怎样复用同一份快照

`computeEmbeddedAlphaUpperBound(...)` 与 barrier 一样，也会映到 sample 空间。  
然后对每个 active pair 读取：

- `contactedTrianglePairs`
- `activeNormals`
- `activeBarycentricWeights`

再计算：

- 当前距离 `d0`
- 相对法向速度 `dDot`

其中 triangle 侧的位移增量会先按 barycentric weight 组合成：

- `closestPoint`
- `closestDelta`

因此当前 self alpha upper bound 与 self barrier 使用的是同一份 active-pair snapshot，而不是两次独立的几何搜索。

## 11. 当前 repo 的验证证据

`tests/core/scene/self_contact_handler_test.cpp` 已经把 self path 的关键语义固定成 committed coverage，包括：

- `NearContactActivationUsesUnsignedDistance`
- `ActivationDistanceThresholdSuppressesPairs`
- `ZeroBandFallsBackToCollisionOnlyWhenNotColliding`
- `ZeroBandFallbackHandlesCollidingPairs`
- `ClosestTargetSelectionKeepsOneTargetPerSample`
- `LegacyPenaltyPathStillWorks`
- `BarrierBuilderProducesFiniteEnergyForNearContactPairs`
- `AlphaUpperBoundShrinksForInwardDirections`
- `AlphaUpperBoundScalesWithAlphaSafety`
- `AlphaUpperBoundAllowsOutwardRecovery`
- `AlphaUpperBoundBlocksFurtherInwardMotionInsideSafeBand`
- `AlphaUpperBoundReturnsOneWithNoActivePairs`

这些测试共同证明了下面几件事：

- self near-contact 判据使用 unsigned distance
- band 改小时 active pair 会被压掉
- 零带宽仍然保留旧碰撞 fallback
- 每个 sample 当前只保留一个目标 triangle
- self barrier 能产生有限、正的能量
- self feasible alpha upper bound 已经按 inward/outward/safe-band 语义工作

再往上，`tests/api/runSimConfig_parse_test.cpp` 的下面几个 smoke 会覆盖 self runtime：

- `RunSimFromConfigCubicDynamicIpcSelfBarrierSmokeTest`
- `RunSimFromConfigCubicDynamicIpcMergedBarrierSmokeTest`
- `RunSimFromConfigCubicDynamicIpcMergedDeterministicSmokeTest`

example 侧则有：

- `examples/pulled-cubic-box-self-ipc/pulled-cubic-box-self-ipc.json`
- `examples/pulled-cubic-box-self-ipc/README.md`

## 12. 这一阶段不展开什么

当前 self path 仍然明确保留下面这些边界：

- 不是 primitive PT/EE IPC
- 不包含 exact CCD-based self feasible alpha upper bound
- 不包含 frictional IPC
- 不包含 inversion-free filter
- active pair 当前仍然是 sampled point-triangle 数据结构

因此 walkthrough 对当前 self IPC 的准确描述应当是：

> sample-based self near-contact barrier runtime

而不是 full primitive self IPC。

上一阶段： [Phase 4](phase-4-feasible-line-search.md)  
下一阶段： [Phase 6](phase-6-tests-and-validation.md)
