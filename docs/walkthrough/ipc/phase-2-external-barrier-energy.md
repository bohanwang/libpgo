# Phase 2：External Barrier Energy

External IPC 的完整 vertical slice：

```text
sampled surface point → near-contact detection → active sample snapshot → barrier energy → feasible alpha upper bound
```

## Active Sample 检测

### 基础几何单位

External handler 建立在 **sampled surface** 上。构造阶段完成 surface triangle sampling、去重、插值矩阵 `interpolationMatrix`。Runtime 每帧变化的是 surface 顶点位置和最近外部目标快照，而不是 sample 集合本身。

Handler 同时持有外部障碍物侧的 BVH、pseudo normal 和 runtime mesh snapshot，`execute(...)` 时直接做最近三角形查询。

### `execute(usurf, externalIpcDhat)` 流程

1. `restP + u` → 当前 surface 顶点 `curP`
2. `computeSamplePosition(...)` → 当前 `sampleCurP`
3. 对每个 sample，遍历所有 external surface：
     - `closestTriangleQuery(...)` 找最近三角形
     - 法向朝向过滤：`tgtNormal · srcNormal > 0` 的候选被排除
     - 计算 `signedDistance = (srcPos - tgtPos) · tgtNormal`
     - 保留 `signedDistance < activationDistance` 的候选
     - 多候选中选 signed distance 最小的
4. 结果：**每个 sample 至多一个最近 external target snapshot**

### Active sample 快照内容

- `constraintCoeffs`、`constraintNormals`、`constraintTargetPositions`、`constraintSignedDistances`
- `barycentricIdx`、`barycentricWeights`、`contactedSamples`、`contactedTriangles`

Barrier energy 和 feasible alpha 都**直接消费这份快照**，不会重新做最近点查询。

这些字段是在 `TriangleMeshExternalContactHandler::execute(...)` 末尾一次性写好的。可以把它理解成：

```text
active contact #ci
  ↕
sample id = contactedSamples[ci]
  ↕
当前样本点由 barycentricIdx[ci] / barycentricWeights[ci] 从真实 DOF 恢复
  ↕
距离平面由 constraintTargetPositions[ci] / constraintNormals[ci] 定义
  ↕
该条约束的权重是 constraintCoeffs[ci]
```

### 每个字段具体含义

- `constraintCoeffs[ci]`
  这条 active sample 约束的标量权重。当前实现里它来自 `sampleWeights[sId]`，本质上是 sample 在表面采样中的面积权重，随后会乘到 barrier energy / gradient / Hessian 上。

- `constraintNormals.segment<3>(ci * 3)`
  这条约束对应的目标法向 `n`。它不是 source sample 自己的法向，而是最近 external triangle 在最近特征处的 pseudo normal。后续距离都按
  $$
  d = (p - p_0)\cdot n
  $$
  里的 `n` 来算。

- `constraintTargetPositions.segment<3>(ci * 3)`
  这条约束冻结下来的目标点 `p_0`，也就是当前帧最近点查询返回的 `closestPosition`。后面的 barrier energy 不会重新查最近点，而是把它当成固定 target snapshot。

- `constraintSignedDistances[ci]`
  `execute(...)` 当下测得的有符号距离
  $$
  d_{\text{snapshot}} = (p_{\text{sample}} - p_0)\cdot n
  $$
  它主要用于调试、日志、测试和检查 active set 内容；真正求能量时会重新用当前迭代位置计算距离，而不是直接复用这个标量。

- `barycentricIdx[ci]`
  一个整数列表，表示“这条 sample 约束会影响哪些真实顶点 / DOF 块”。它来自 `interpolationMatrix` 对应该 sample 行的非零列索引，代码里是把列号除以 3 得到顶点号。
  如果 sample 就是表面顶点，它通常只含一个顶点；如果 sample 在三角形内部，通常会含 3 个顶点；如果 surface 又嵌到更大的 embedded DOF 中，这里还可能展开成 embedding 后的一组顶点。

- `barycentricWeights[ci]`
  与 `barycentricIdx[ci]` 一一对应的插值权重。后续恢复 sample 位置时，用
  $$
  p = \sum_j w_j\,p_j
  $$
  把真实 DOF 上的位置插值回 sample 点；gradient / Hessian 也沿同一组权重散回全局系统。

- `contactedSamples[ci]`
  这条 active contact 对应的是哪个 sample ID。它把“约束编号 `ci`”连回 handler 内部的 sample 空间，用于 feasible alpha upper bound 阶段重新取该 sample 的当前位移和搜索方向。

- `contactedTriangles[ci]`
  这条 active sample 最初来自哪个 source surface triangle。它主要是追踪信息，便于调试、可视化和定位是哪片三角形激活了外部接触；barrier 核心计算本身更直接依赖的是 `contactedSamples` 和 barycentric 映射。

### 它们如何一起定义一条外部 barrier 约束

对第 `ci` 条 active sample：

1. 用 `barycentricIdx[ci]` 和 `barycentricWeights[ci]` 从当前 DOF 恢复 sample 点位置 `p`
2. 用 `constraintTargetPositions[ci]` 取冻结目标点 `p_0`
3. 用 `constraintNormals[ci]` 取冻结法向 `n`
4. 计算当前距离 `d = (p - p_0) · n`
5. 用 `constraintCoeffs[ci]` 作为该条 barrier 项的权重

所以这组快照本质上定义的是：

```text
一个“当前可动 sample 点” vs 一个“本帧冻结的外部目标平面”
```

而不是“每次评估都重新做 closest point query”的动态接触模型。

## `PointPenetrationBarrierEnergy`

### 距离模型

每个 active sample 的局部距离：

$$d = (p - p_0) \cdot n$$

- $p$：由 `setComputePosFunction(...)` + barycentric weights 从当前 DOF 恢复
- $p_0$：handler 缓存的 `constraintTargetPositions`
- $n$：handler 缓存的 `constraintNormals`

典型的 **sampled point vs frozen target plane snapshot** 语义。

### Barrier 函数

使用 `BarrierFunctions::logBarrierEnergy`（定义在 `barrierFunction.{h,cpp}`）：

$$b(d, \hat{d}) = -(d - \hat{d})^2 \ln\!\left(\frac{d}{\hat{d}}\right)$$

- $d \geq \hat{d}$：返回零（支撑区间外）
- $d < \hat{d}$：barrier 激活
- $d \leq 0$：被 clamp 到 `positiveZero` 避免 NaN

总 external barrier energy = $\sum_i w_i \cdot \kappa \cdot b(d_i, \hat{d})$

其中 $w_i$ 是 sample weight，$\kappa$ 是 `externalIpcKappa`。

### `normalizePositiveZero(dhat, positiveZero)`

把 `positiveZero` 约束在安全上界内（`min(positiveZero, dhat * 0.5)`），未显式给值时用 `max(1e-10, 1e-6 * dhat)` 作为推荐值。

它的角色是一个**数值安全阈值**，不是额外的物理接触参数。

原因是 barrier 函数内部要计算 `log(d / dhat)`。当 sample 已经明显穿透、导致 `d <= 0` 时，直接取对数会产生 `NaN` 或 `-inf`。因此实现不会把非正距离原样送进 barrier，而是先把它钳到一个很小的正数上；这个“很小的正数”就是 `positiveZero`。

可以把它理解成：

- `dhat`：决定 barrier 从哪里开始激活的几何阈值
- `kappa`：决定 barrier 罚得多重
- `positiveZero`：当距离已经坏到 `d <= 0` 时，保证 barrier 仍然可计算的**数值保险丝**

### Gradient 装配

对第 $i$ 条 active sample，barrier energy 对距离 $d$ 的标量导数是：

$$g_i = w_i \cdot \kappa \cdot b'(d_i, \hat{d})$$

其中 $b'$ 是 `logBarrierGradient`：

$$b'(d, \hat{d}) = -2(d - \hat{d})\ln\!\left(\frac{d}{\hat{d}}\right) - \frac{(d - \hat{d})^2}{d}$$

由于 $d = (p - p_0) \cdot n$，距离对 sample 点位置 $p$ 的梯度就是法向本身：

$$\frac{\partial d}{\partial p} = n$$

所以 barrier 对 sample 点的 3D 梯度为：

$$\nabla_p E_i = g_i \cdot n$$

但 sample 点不是自由度 — 它由真实 DOF 通过 barycentric 插值恢复：$p = \sum_j w_j \, p_j$。因此需要把梯度散回每个参与顶点：

$$\nabla_{p_j} E_i = w_j \cdot g_i \cdot n$$

代码中对应的逻辑（`gradient()` 方法）：

```cpp
// 对每个 active sample ci
const double barrierGradient = constraintCoeffs[ci] * BF::logBarrierGradient(distance, dhat, positiveZero);
const ES::V3d gradSample = n * barrierGradient;  // 3D sample-space gradient

for (int vi = 0; vi < barycentricIdx[ci].size(); vi++) {
    const int vid = barycentricIdx[ci][vi];
    const double w = barycentricWeights[ci][vi];
    grad.segment<3>(vid * 3) += gradSample * w;   // scatter 到全局 DOF
}
```

最后整个 `grad` 再乘以全局系数 `coeffAll`（即 $\kappa$）。并行通过 per-vertex `spin_mutex` 保护写入。

### Hessian 装配

Barrier 对距离的标量二阶导是：

$$h_i = w_i \cdot \kappa \cdot b''(d_i, \hat{d})$$

其中 $b''$ 是 `logBarrierHessian`：

$$b''(d, \hat{d}) = -2\ln\!\left(\frac{d}{\hat{d}}\right) - \frac{4(d - \hat{d})}{d} + \frac{(d - \hat{d})^2}{d^2}$$

由于距离模型是 $d = (p - p_0) \cdot n$（对 $p$ 是线性的），Hessian 不含二阶几何项，直接是：

$$\frac{\partial^2 E_i}{\partial p \,\partial p} = h_i \cdot n n^T$$

这是一个 $3 \times 3$ 秩 1 矩阵。散回真实 DOF 后，顶点 $j$ 和顶点 $k$ 之间的 $3 \times 3$ 块为：

$$H_{jk}^{(i)} = w_j \cdot w_k \cdot h_i \cdot n n^T$$

代码中对应的逻辑（`hessian()` 方法）：

```cpp
// 对每个 active sample ci
const double barrierHessian = constraintCoeffs[ci] * BF::logBarrierHessian(distance, dhat, positiveZero);
const ES::M3d nnT = ES::tensorProduct(n, n) * barrierHessian;  // h_i * n n^T

for (int vi = 0; vi < barycentricIdx[ci].size(); vi++) {
    const double wi = barycentricWeights[ci][vi];
    for (int vj = 0; vj < barycentricIdx[ci].size(); vj++) {
        const double wj = barycentricWeights[ci][vj];
        const ES::M3d hLocal = nnT * wi * wj;  // 3x3 block for (vi, vj)

        // scatter 9 个标量到全局稀疏 Hessian
        for (int dofi = 0; dofi < 3; dofi++)
            for (int dofj = 0; dofj < 3; dofj++)
                hess(vi*3+dofi, vj*3+dofj) += hLocal(dofi, dofj);
    }
}
```

### Hessian 稀疏模板

构造时预建稀疏模板：遍历所有 active sample 的 `barycentricIdx`，为每对 $(v_i, v_j)$ 的 $3 \times 3$ 块预分配 triplet。这样运行时只需要按 `hessianEntryMap` 直接写入 `valuePtr()`，不需要动态插入。

## External Feasible Alpha Upper Bound

Handler 提供 `computeSurfaceAlphaUpperBound(...)` 和 `computeEmbeddedAlphaUpperBound(...)`，两者都映到 sample 空间调用 `computeSampleAlphaUpperBound(...)`。

对每个 active sample，计算：

- 当前距离 $d_0$、法向闭合速度 $\dot{d}$
- 安全距离 $d_{\text{safe}}$ = `normalizePositiveZero(dhat, -1.0)`

分支逻辑：

| 条件 | 行为 |
| --- | --- |
| $d_0 \leq d_{\text{safe}}$ 且 $\dot{d} < 0$（继续 inward） | 返回 $\alpha = 0$ |
| $d_0 \leq d_{\text{safe}}$ 且 $\dot{d} \geq 0$（outward） | 允许 recovery |
| $\dot{d} < 0$ | 收缩 $\alpha$ 使距离不低于 $d_{\text{safe}}$ |
| $\dot{d} \geq 0$ | 不收缩 |

最终结果乘以 `ipcAlphaSafety`（solver-side margin）。

## 验证

- 数学：`pointPenetrationBarrierEnergy_test.cpp` — 支撑区间、有限差分、NaN clamp
- Handler：`contact_embedding_test.cpp` — near-contact 激活、alpha 上界的解析验证、safe-band 行为
- Runtime：`runSimConfig_parse_test.cpp` — API smoke 检查日志 `# external active samples: ...`

---

上一阶段：[Phase 1](phase-1-config-and-dispatch.md)
下一阶段：[Phase 3](phase-3-self-near-contact-active-set.md)
