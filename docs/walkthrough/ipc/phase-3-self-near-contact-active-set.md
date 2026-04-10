# Phase 3：Self Near-Contact Active Set

Self IPC 的 vertical slice：

```text
BVH query band → candidate triangle pairs → sample seed refine → local search → active point-triangle pairs → self barrier energy → self feasible alpha
```

核心复杂度不在 barrier 公式（与 external 类似），而在如何从 BVH 候选稳定地生成 active pair。

## Self vs External：为什么 self contact 更难

### External Contact 的简单性

External barrier 能做到"一趟查询 + 冻结快照"，根本原因是两侧几何完全独立：

1. **几何独立**：source mesh 和 external obstacle 是不同对象，顶点集合不重叠，不需要排除"自己碰自己"
2. **单向检测**：只需要 source sample → external surface 的最近点查询，external surface 本帧不动
3. **冻结 target**：最近点 $p_0$ 和法向 $n$ 在检测后固定为常量，barrier 中只有 source 侧的 DOF 参与梯度/Hessian
4. **无拓扑过滤**：不存在"sample 和 target 共享顶点"的情况，所有候选都是合法的

### Self Contact 引入的额外复杂度

同一个 mesh 做自碰撞查询时，上述每一条都不再成立：

**1. 共享顶点问题**

Source sample 和 target triangle 来自同一个 mesh。如果 sample 本身就是 target triangle 的一个顶点（或来自邻接三角形），距离天然为零但不代表穿透。必须显式过滤掉拓扑上相邻的 pair，否则 active set 里全是虚假接触。

实现中通过 `excludedTriangles`（每个 sample 预计算的邻接排除集）和 `contactedTrianglePairs` 的组合检查来处理。

**2. Broad Phase 更贵**

External 只需要一棵 BVH（external surface），做 point → tree 查询。Self contact 需要同一棵 BVH 的 **self-query**：遍历 BVTT（BV Testing Tree）的所有非自身子树组合，候选数量量级更高。Detection 层引入 `BVTT frontier` + TLS buffer + normalize/sort/unique 来控制开销。

**3. 三级细化 pipeline**

External 的检测是一步到位：最近三角形查询 → 法向过滤 → 距离过滤。Self 需要三级：

| 级别 | 操作 | 目的 |
| --- | --- | --- |
| BVH query band | `trianglesWithinQueryBand` | AABB 级快速排除，产出候选 triangle pair |
| Seed refine | unsigned distance < $\hat{d}$ | 对每个 sample 从候选中选唯一最近 seed triangle |
| Local search | 从 seed 沿邻接扩展 | seed 未必是全局最近，局部搜索修正最终 target |

多出的 seed → local search 两步是 self contact 特有的，因为 BVH band 返回的候选太粗，直接用会导致 active set 不稳定。

**4. Target 不再冻结**

External 的 target 是外部障碍物，本帧固定。Self 的 target triangle 三个顶点也是优化变量 — barrier energy 对 **source 和 target 双方的 DOF** 都有梯度。这意味着：

- Gradient 需要同时散回 source sample 和 target triangle 三顶点的 DOF
- Hessian 从 external 的 rank-1 $nn^T$（只涉及 source barycentric 组）变成 4 组 sample 之间的耦合块
- 稀疏模板更密：每条 active pair 关联 $4 \times 4 = 16$ 个 sample-pair 的 $3 \times 3$ block（而 external 只有 source barycentric 组的 $k \times k$ block）

**5. 法向稳定性**

External 的法向来自外部 mesh 的 pseudo normal，不受优化影响。Self 的法向取自 `sample - closestPoint` 方向，当距离很小时这个方向不稳定，需要 fallback 到 triangle 几何法向。

### 总结对比

| 维度 | External | Self |
| --- | --- | --- |
| 几何关系 | 不同对象，无共享顶点 | 同一 mesh，需拓扑过滤 |
| Broad phase | point → external BVH | 同一 BVH self-query (BVTT frontier) |
| 检测层数 | 1 级（最近点查询） | 3 级（BVH band → seed → local search） |
| Target 角色 | 冻结常量（不参与优化） | 优化变量（参与 gradient/Hessian） |
| 每条约束涉及的 DOF | 1 个 source sample 的 barycentric 组 | 1 source + 3 target corner = 4 组 |
| Hessian block 结构 | $k \times k$ blocks（$k$ = source 插值顶点数） | $4k \times 4k$ blocks（4 组 sample 交叉） |
| 法向来源 | external pseudo normal（稳定） | sample-closestPoint 方向（需 fallback） |

## 分层架构

Self path 明确拆成两层：

| 层 | 类 | 职责 |
| --- | --- | --- |
| Detection | `TriangleMeshSelfContactDetection` | BVH traversal、query band、candidate triangle pair 生成 |
| Handler | `TriangleMeshSelfContactHandler` | sample seed refine、local search、barrier 构造、alpha upper bound |

基础几何单位仍然是 **sampled surface point**，与 external 一致。

## Broad Phase：BVH Query Band

### BVTT 是什么

这里的 `BVTT` 指 **BV Testing Tree**。在这份实现里，它不是单独新建的一棵树，而是“对同一棵 BVH 做 self-query 时，正在遍历的 node-pair 空间”。

可以把它理解成：

- `BVH` 是单棵层次包围盒树
- `BVTT` 是把这棵树和它自己配对后，从 `(root, root)` 开始逐层展开的遍历过程

它的职责不是直接产生最终 active pair，而是先把不可能接近的 node-pair 剪掉，只把可能进入 query band 的部分继续下钻到 triangle 级别。

一个很小的遍历图：

```text
(root, root)
    |
    v
(child, child) ...
    |
    v
(leaf, leaf)
    |
    v
(triangle, triangle) pairs
```

Detection 复用现有 `Mesh::TriMeshBVTree`，没有引入新的树类型。IPC 改动集中在 **query 语义**：

- `activationDistance ≤ 0`：保持旧的 collision-only 判据（`DCDIntersect` / `CCDIntersect`）
- `activationDistance > 0`：改用 `nodesWithinQueryBand(...)` / `trianglesWithinQueryBand(...)`，回答"哪些 triangle pair 进入了 $\hat{d}$ 带宽"

Detection 额外管理 node/triangle AABB 快照、BVTT frontier、TLS triangle pair buffer，形成可并行的候选生成 pipeline：

```text
updateBoundingVolumes → refresh AABB → traverse BVTT frontier → emit pairs → normalize/sort/unique
```

输出 `candidateTrianglePairs` 已经是去重后的 broad-phase 候选集，供 handler 继续细化。

## Handler：从候选到 Active Pair

### `execute(u0, activationDistance)` 分叉

- `activationDistance ≤ 0`：回退到旧碰撞路径（`executeDCD()` → `handleContactDCD()`）
- `activationDistance > 0`：进入 `executeNearContactDCD(activationDistance, maxSearchTriangles)`

### Sample Seed Refine

对每个 candidate triangle pair，检查该 pair 中的 sample point 到 target triangle 的 **unsigned distance**：

$$\text{dist}^2 = \texttt{getSquaredDistanceToTriangle}(p, v_a, v_b, v_c)$$

- 只保留 $\text{dist}^2 < \hat{d}^2$ 的候选
- 每个 sample 只保留距离最小的一个 seed triangle
- 结果：`sample → best seed triangle`（避免一对多导致 active set 爆炸）

### Local Search

Seed 只是候选起点，不保证已经是最终最近的 target triangle。
- `finalizeActivePairsFromSeeds(...)` 会从 seed 出发，沿 triangle 邻接做一个小范围局部搜索，不断检查相邻 triangle 的最近点距离和 barycentric weight；
- 如果邻居更近，就更新 best result。priority queue 用来优先扩展更近的候选，并在“待搜索候选已经明显差于当前 best”时提前停止。
- 最后还会过滤掉拓扑上不合法的 pair，例如 sample 本身属于 target triangle，或命中 `excludedTriangles` 的情况。

### Active Pair 缓存

最终写入的缓存：

| 字段 | 含义 |
| --- | --- |
| `contactedTrianglePairs[i]` | `[src_sample, tgt_corner0, tgt_corner1, tgt_corner2]` — 4 个 sample id |
| `activeClosestPoints` | 最近点坐标 |
| `activeNormals` | 优先取 `sample - closestPoint` 方向；过短时回退到 triangle 几何法向 |
| `activeBarycentricWeights` | 最近点的 barycentric 坐标 |
| `activeDistances` | unsigned distance |

这份缓存被 barrier energy 和 alpha upper bound **共享消费**。

## `PointTrianglePairBarrierEnergy`

### 距离模型

每条 active pair 包含 4 个 sample：source point（index 0）和 target triangle 三个 corner（index 1, 2, 3）。距离定义为：

$$d = n \cdot \bigl(p_0 - (w_0\, p_1 + w_1\, p_2 + w_2\, p_3)\bigr)$$

- $p_0$：source sample 当前位置
- $p_1, p_2, p_3$：target triangle 三个 corner sample 的当前位置
- $n$：active pair 快照中的法向（`constraintNormals`）
- $w_0, w_1, w_2$：target 最近点的 barycentric 坐标（`constraintBarycentricWeights`）

与 external 的关键区别：**target 侧的 $p_1, p_2, p_3$ 也是优化变量**，不是冻结常量。

### Energy

Barrier 函数与 external 相同：

$$E_i = w_i \cdot \kappa \cdot b(d_i, \hat{d}), \quad b(d, \hat{d}) = -(d - \hat{d})^2 \ln\!\left(\frac{d}{\hat{d}}\right)$$

- $w_i$：source sample 的面积权重（`sampleWeights`）
- $\kappa$：全局系数（`coeffAll`）
- $d \geq \hat{d}$ 时跳过，$d \leq 0$ 时 clamp 到 `positiveZero`

### Gradient

定义 4 个 local scale 系数，表示距离对每个 sample 位置的偏导方向：

$$s_0 = +1, \quad s_1 = -w_0, \quad s_2 = -w_1, \quad s_3 = -w_2$$

则距离对第 $k$ 个 sample 位置的梯度为：

$$\frac{\partial d}{\partial p_k} = s_k \cdot n$$

Barrier 对第 $k$ 个 sample 的 3D 梯度：

$$\nabla_{p_k} E_i = w_i \cdot \kappa \cdot b'(d_i, \hat{d}) \cdot s_k \cdot n$$

每个 sample 不是直接自由度，而是通过 `sampleEmbeddingWeights` 从真实 DOF 插值得到。设第 $k$ 个 sample 的插值展开为 $p_k = \sum_j \alpha_j \, x_j$，则散回全局 DOF：

$$\nabla_{x_j} E_i = \alpha_j \cdot s_k \cdot w_i \cdot \kappa \cdot b'(d_i, \hat{d}) \cdot n$$

代码中对应逻辑：

```cpp
const double barrierGradient = pairWeight * BF::logBarrierGradient(distance, dhat, positiveZero);
const double localScale[4]   = {1.0, -baryW[0], -baryW[1], -baryW[2]};

for (int localIdx = 0; localIdx < 4; ++localIdx) {
    const ES::V3d localGrad = normal * (barrierGradient * localScale[localIdx]);
    // 遍历 sampleEmbeddingWeights 的非零列，散回全局 DOF
    for (InnerIterator it(...); it; ++it) {
        grad.segment<3>(offset) += localGrad * it.value();
    }
}
```

并行写入通过 per-DOF `spin_mutex` 保护。

### Hessian

距离对 $p$ 是线性的，因此没有二阶几何项。Barrier 的二阶导完全来自标量 barrier Hessian：

$$h_i = w_i \cdot b''(d_i, \hat{d})$$

对第 $k$ 个和第 $l$ 个 sample，Hessian 的 $3 \times 3$ block 为：

$$H_{kl}^{(i)} = h_i \cdot s_k \cdot s_l \cdot n n^T$$

这是一个秩 1 矩阵，但因为涉及 4 个 sample（而非 external 的 1 个），总共产生 $4 \times 4 = 16$ 个 sample-pair 的 $3 \times 3$ block。

进一步展开到真实 DOF：设 sample $k$ 的第 $i$ 个插值顶点权重为 $\alpha_i$，sample $l$ 的第 $j$ 个插值顶点权重为 $\beta_j$，则全局 Hessian 的 $3 \times 3$ block：

$$H_{\text{global}}^{(i,j)} = \alpha_i \cdot \beta_j \cdot s_k \cdot s_l \cdot h_i \cdot n n^T$$

代码中对应逻辑：

```cpp
const double barrierHessian = pairWeight * BF::logBarrierHessian(distance, dhat, positiveZero);
const ES::M3d nnT           = ES::tensorProduct(normal, normal) * barrierHessian;
const double localScale[4]  = {1.0, -baryW[0], -baryW[1], -baryW[2]};

for (int vi = 0; vi < 4; ++vi) {
    for (int vj = 0; vj < 4; ++vj) {
        const ES::M3d Kij = nnT * (localScale[vi] * localScale[vj]);
        // 双重遍历两个 sample 的 embedding weights，散回全局稀疏矩阵
        // 通过预建的 entryMap 直接写入 valuePtr()
    }
}
```

### 稀疏模板

构造时预遍历所有 active pair 的 4 个 sample × 4 个 sample 的 embedding weight 非零列组合，为每个 $(row, col)$ 的 $3 \times 3$ block 分配 triplet。运行时通过 `entryMap` 直接定位写入位置。

与 external 对比：external 每条约束只涉及 1 个 source sample 的 barycentric 组（通常 1–3 个顶点），稀疏模板很小。Self 每条约束涉及 4 个 sample 各自的 embedding 展开，交叉组合后 triplet 数量大得多。

### 建模边界

这是 **sample-based point-triangle linearization**，不是 primitive-level exact IPC。法向和 barycentric 坐标在检测时冻结，每次 Newton 迭代重新计算距离但不更新几何投影。

## Self Feasible Alpha Upper Bound

### 入口

`computeEmbeddedAlphaUpperBound(currentU, du, alphaSafety, dSafe)` 先通过 `interpolationMatrix` 把全局 DOF 位移和搜索方向映射到 sample 空间，然后调用 `computeSampleAlphaUpperBound(...)`。

### 每条 active pair 的计算

对第 $i$ 条 active pair（4 个 sample：source point + target triangle 三顶点），读取快照中的法向 $n$ 和 barycentric 权重 $w$：

1. **恢复当前位置和搜索方向**

   $$p = \text{restP}_{s} + u_{s}, \quad dp = du_{s}$$

   $$q = w_0 t_0 + w_1 t_1 + w_2 t_2, \quad dq = w_0\,dt_0 + w_1\,dt_1 + w_2\,dt_2$$

   其中 $t_k = \text{restP}_{k} + u_{k}$ 是 target corner 的当前位置，$dt_k$ 是搜索方向。$q$ 和 $dq$ 分别是 target 最近点的当前位置和搜索方向。

   与 external 的区别：external 的 target 不动所以 $dq = 0$；self 的 target 也在运动，$dq$ 非零。

2. **当前距离和法向闭合速度**

   $$d_0 = (p - q) \cdot n, \quad \dot{d} = (dp - dq) \cdot n$$

   $\dot{d} < 0$ 表示 source 和 target 在沿法向靠近（closing），$\dot{d} \geq 0$ 表示远离或平行。

3. **分支逻辑**

   | 条件 | 行为 | 原因 |
   | --- | --- | --- |
   | $d_0 \leq d_{\text{safe}}$ 且 $\dot{d} \leq 0$ | 立即返回 $\alpha = 0$ | 已在安全带内且继续靠近，必须完全阻止 |
   | $d_0 \leq d_{\text{safe}}$ 且 $\dot{d} > 0$ | 跳过（不收缩） | 已在安全带内但正在远离，允许 recovery |
   | $d_0 > d_{\text{safe}}$ 且 $\dot{d} < 0$ | $\alpha \leq \frac{d_0 - d_{\text{safe}}}{-\dot{d}}$ | 收缩步长使距离不低于 $d_{\text{safe}}$ |
   | $d_0 > d_{\text{safe}}$ 且 $\dot{d} \geq 0$ | 不收缩 | 远离方向，无需干预 |

   第三种情况的推导：线搜索后距离为 $d_0 + \alpha \dot{d}$，要求 $d_0 + \alpha \dot{d} \geq d_{\text{safe}}$，由于 $\dot{d} < 0$，解得 $\alpha \leq (d_0 - d_{\text{safe}}) / (-\dot{d})$。

### 全局聚合

遍历所有 active pair，取最小的 $\alpha$ 上界。如果没有任何 pair 处于 closing 方向，直接返回 $1.0$。否则：

$$\alpha_{\text{final}} = \text{clamp}\bigl(\alpha_{\min} \cdot \text{alphaSafety},\; 0,\; 1\bigr)$$

`alphaSafety` 是 solver 侧的额外安全裕度（通常略小于 1），防止数值误差导致恰好踩到安全带边界。

## 验证

`self_contact_handler_test.cpp` 覆盖：unsigned distance 判据、threshold suppression、zero-band fallback、one target per sample、barrier finite energy、inward/outward/safe-band alpha 行为。

---

上一阶段：[Phase 2](phase-2-external-barrier-energy.md)
下一阶段：[Phase 4](phase-4-dynamic-incremental-potential.md)
