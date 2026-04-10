# Phase 3：self-contact 的 near-contact active set

## 1. 这阶段要回答什么

这一阶段专门回答下面四个问题：

1. self-contact 在当前 repo 里怎样从 triangle-pair 候选细化成最终 active pair。
2. 为什么 current self IPC 仍然是 sample-based，而不是 primitive PT/EE。
3. `TriangleMeshSelfContactDetection` 和 `TriangleMeshSelfContactHandler` 当前分别负责什么。
4. self barrier energy 和 self feasible alpha upper bound 共享的是哪一批缓存。

如果只说“self 也有 barrier 了”，会掩盖当前 self 路径真正复杂的地方。  
Phase 3 的核心不在 barrier 标量公式，而在：

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

所以 Phase 3 实际上是 self IPC 的 geometry-and-runtime composition phase。

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

### 4.3 BVH 在当前 self IPC 路径里到底改了什么

如果只看 `TriangleMeshSelfContactDetection` 这个类名，再看到文档里频繁出现 “BVH query band”，很容易形成一个过强的印象：

> 为了接 self IPC，repo 重新实现了一棵新的 BVH tree。

这其实不是当前代码的真实形状。  
当前 self IPC 对 BVH 的改动，更准确的说法是：

> repo 复用了现有 `Mesh::TriMeshBVTree`，但把它在 self-contact broad phase 中回答的问题，从 collision-only 扩展成了 near-contact candidate generation。

也就是说，这里的核心变化不在 “树怎么建”，而在 “树怎样被查询、查询结果怎样被组织并继续往下游传”。

### 4.4 当前复用的仍然是已有 `TriMeshBVTree`

在 `TriangleMeshSelfContactDetection` 的构造阶段，当前实现做的仍然是：

- 持有一个 `Mesh::TriMeshBVTree bvTree`
- 对 rest mesh 调 `buildByInertiaPartition(triangleMeshRef)`
- 初始化 node AABB 缓存
- 初始化 triangle AABB 缓存

这说明当前 self IPC 并没有引入一个新的、IPC 专属的树类型。  
tree 本体仍然是 repo 现有几何层里的通用 BVH；IPC 相关的逻辑并没有被塞回 `boundingVolumeTree.*` 去污染通用数据结构。

这点很重要，因为它决定了当前设计的职责边界：

- `boundingVolumeTree`
  继续负责通用 mesh hierarchy 与 bounding volume 更新
- `TriangleMeshSelfContactDetection`
  在 contact 层决定怎样消费这棵树、怎样解释 “命中”

因此如果要描述当前仓库中的改动，最准确的表述不是：

```text
我们把 BVH tree 重写成了 IPC 版本
```

而是：

```text
我们在现有 BVH traversal 之上，增加了 IPC self near-contact 需要的查询语义和候选组织方式
```

### 4.5 当前新增的是 query semantics，而不是建树算法

在 legacy self-contact 语义里，broad phase 主要回答的是：

```text
哪些 triangle pair 已经发生碰撞，或者可能在 CCD 中发生碰撞
```

这对 penalty / collision fallback 足够，但对 IPC self barrier 不够。  
因为 self barrier 需要的是：

```text
哪些 triangle pair 虽然还没有真正相交，但已经进入了 selfIpcDhat 的近接触带宽
```

所以当前 repo 里真正新增的不是另一套树构建策略，而是下面这些 query-side 语义：

- `execute(..., activationDistance)` 这一组 overload
- `nodesWithinQueryBand(...)`
- `trianglesWithinQueryBand(...)`
- `AABBWithinBand(...)`
- `SweptAABB(...)`

这几个接口共同表达的意思是：

- 当 `activationDistance <= 0` 时
  broad phase 仍然保持旧的 collision-only 判据
- 当 `activationDistance > 0` 时
  broad phase 开始回答 near-contact band 查询

因此当前 self IPC 对 BVH 的关键升级，不是让 tree “更复杂”，而是让同一棵 tree 能服务两种不同的问题：

1. 旧的 collision detection 问题
2. 新的 near-contact candidate generation 问题

### 4.6 IPC 语义为什么放在 detection 层，而不是放回 tree 层

当前实现把 IPC 语义保留在 `TriangleMeshSelfContactDetection`，而没有去改通用 `boundingVolumeTree` 接口，本质上是一个很清楚的分层决定。

原因是这里的 “band” 不是 BVH 的普适概念，而是当前 self IPC runtime 的 contact 语义。  
它依赖的是：

- `activationDistance`
- DCD 还是 CCD
- self-contact 当前帧的 `positions0 / positions1`
- broad phase 最终要输出的是 candidate triangle pair，而不是 generic query hit

这些都属于 contact runtime，而不是通用几何容器的职责。

把这层语义放在 detection 层，有两个直接结果：

1. `TriMeshBVTree` 仍然保持通用性  
   其他模块不需要知道什么是 `selfIpcDhat`。
2. self IPC 可以在不改 tree 本体的情况下演化自己的查询语义  
   例如当前的 near-contact band、零带宽 fallback、CCD swept band，都可以只在 detection 层收敛。

这也是为什么当前 repo 的 self IPC 看起来像是在 “改 BVH”，但实际上改动集中在 contact detection 代码，而不是 mesh hierarchy 代码。

### 4.7 当前 detection 对这棵树外挂了哪些 runtime 组织

虽然 tree 本体没有被 IPC 重写，但 `TriangleMeshSelfContactDetection` 的 runtime 组织已经明显是为 self IPC broad phase 服务的。  
当前类内部除了 `bvTree` 本身，还额外管理了：

- `lastAABBs` / `curAABBs`
  node 层 AABB 快照
- `lastTriangleAABBs` / `curTriangleAABBs`
  triangle 层 AABB 快照
- `travseralFrontier0` / `travseralFrontier1`
  BVTT frontier
- `trianglePairBufferTLS`
  线程本地 triangle pair 收集缓冲
- `bvttNodeBufferTLS`
  线程本地下一层 frontier 缓冲
- `parallel_threshold`
  控制串行 / 并行 traversal 的切换点

这意味着当前 repo 的改动不是简单地 “多加一个 band 判据” 就结束了。  
它同时把 broad phase 的运行形态整理成了一个稳定的 runtime pipeline：

```text
updateBoundingVolumes
-> refresh node/triangle AABB snapshots
-> traverse BVTT frontier
-> emit candidate triangle pairs into TLS buffers
-> merge / normalize / sort / unique
```

因此你可以把当前 self IPC 对 BVH 的改进理解成两层：

1. 语义层  
   从 collision-only 查询变成 near-contact 查询
2. 运行层  
   把这类查询组织成一个可并行、可去重、可供 handler 消费的稳定候选生成器

### 4.8 broad phase 的输出也被改成了更适合 IPC 下游的形状

旧的 self collision 视角更关心的是：

- 有没有碰撞
- 哪些 pair 需要进入 CCD test

但当前 self IPC broad phase 对下游真正有价值的不是 “碰撞是否成立” 本身，而是：

- 哪些 triangle pair 值得继续做 sample seed refine

所以 detection 当前返回的是：

- `candidateTrianglePairs`

它已经不再强调 “colliding” 这个命名，而是强调：

```text
这是一份 broad-phase 候选
```

这和后面的 handler 角色正好对应起来：

- detection
  负责把 BVH 输出压成可管理的 triangle-pair 候选集
- handler
  负责把这份候选集继续细化成真正的 sampled point-triangle active pair

换句话说，BVH 在当前 self IPC 链路里的改进，不只是 “查得更宽一点”，而是：

> 它的输出语义被重新定义成了 active-set pipeline 的第一层中间表示。

### 4.9 这一层改造为什么值得单独讲清楚

如果不把这一段单独讲清楚，读者很容易在两种误解之间来回摇摆：

- 误解一：
  以为 repo 只是给旧 self collision 检测多加了一个阈值
- 误解二：
  以为 repo 为 IPC 另起炉灶写了一套新的 BVH library

当前代码真实处在这两者之间：

- 不是只多了一个阈值，因为 query 语义、缓存组织、输出形状都变了
- 也不是重写 BVH library，因为建树与 bounding volume 更新仍然复用现有 `TriMeshBVTree`

因此这一阶段里更准确的总结应该是：

> 当前 self IPC 对 BVH 的改进，本质上是把现有 `TriMeshBVTree` 从 collision-only broad phase，提升成了一个 near-contact-aware 的候选生成层。

有了这个心智模型，下面第 5 节再去看 `nodesWithinQueryBand(...)`、`trianglesWithinQueryBand(...)` 和 `candidateTrianglePairs` 的具体行为时，就不会把代码读成 “重写一棵树”，也不会把它误读成 “只是旧逻辑上加了个 if”。

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

上一阶段： [Phase 2](phase-2-external-barrier-energy.md)  
下一阶段： [Phase 4](phase-4-dynamic-incremental-potential.md)
