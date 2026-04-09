# Phase 06：Full Runtime 与 Example 场景

## 1. 这阶段要回答什么

这一阶段专门回答下面四个问题：

1. cubic mesh 现在怎样真正走进 `runSimCore`？
2. `BarycentricCoordinates`、`GenerateMassMatrix`、surface embedding `W` 在 cubic 路径下怎样复用？
3. runtime setup 在 cubic 路径里已经推进到了哪些层？
4. repo 里的 cubic example 分别说明了哪类运行能力已经落地，边界又在哪里？

## 2. 这阶段在整条主线里解决什么

Phase 05 已经说明 cubic 路径能进入 global energy assembly。  
Phase 06 要解决的是：

> 这些 assembly 能力怎样真正进入 `runSimFromConfig(...)`，并和 surface mesh、质量矩阵、接触、时间积分一起组成可运行的 runtime 主链。

也就是说，这一步不再只是“理论上 solver 可以消费 cubic”，而是“入口函数已经怎样消费 cubic”。

## 3. 当前 repo 已经完成了什么

当前 `runSimCore` 已经能在 cubic 路径下完成下面这些 runtime setup：

- 读取 `CubicMesh`
- 创建 `SimulationMesh`
- 读取并缩放 surface mesh
- 计算 surface 到 volumetric DOFs 的嵌入矩阵 `W`
- 创建 `DeformationModelManager`
- 创建 `DeformationModelAssembler`
- 创建 `DeformationModelEnergy`
- 初始化 plastic 参数向量
- 计算初始 stiffness `K`
- 生成质量矩阵 `M`
- 进入 dynamic 分支时继续组织接触处理器

这说明当前 repo 的 cubic 已经不再只是：

- element-level capability
- manager-level capability
- assembler smoke capability

而是已经进入 shared volumetric runtime。

## 4. `runSimFromConfig(...)` 当前的结构应该怎样读

### 4.1 先读 config

Phase 02 已经解释过 parser 语义。到 `runSimFromConfig(...)` 真正开始执行时，关键输入已经变成：

- `mesh`
- `scene`
- `simulation`
- `contact`
- `solver`
- `runtime`

这一步意味着 cubic 路径不再由字符串或文件扩展名驱动，而是由已经解析好的结构化配置驱动。

### 4.2 runtime 先设置并行环境

当前 runtime 一开始就会根据：

- `std::thread::hardware_concurrency()`
- `runtime.deterministicMode`

设置 TBB 的：

- `max_allowed_parallelism`
- `thread_stack_size`

这说明 cubic 路径从进入 runtime 的第一步起，就已经和 tet 共享同一套执行环境，而不是另开一套专用入口。

## 5. concrete mesh load 怎样收口成 shared volumetric path

### 5.1 concrete 分叉

当前实现会先按 `mesh.volumetricMeshType` 分叉：

- `TET` -> `TetMesh`
- `CUBIC` -> `CubicMesh`

### 5.2 在 concrete mesh 上先做统一几何变换

无论是 tet 还是 cubic，当前 runtime 都会先在 concrete mesh 上遍历顶点并应用：

```cpp
vertex *= mesh.scale;
```

这个顺序很重要。它意味着：

- `scale` 的语义属于 scene / volumetric 几何层
- `SimulationMesh` 看到的就是已经缩放后的坐标

### 5.3 立刻创建 solver-side `SimulationMesh`

接下来分别调用：

- `SimulationMesh::createFromTetMesh(...)`
- `SimulationMesh::createFromCubicMesh(...)`

### 5.4 收口到通用 `VolumetricMesh*`

最后再把 concrete mesh move 到：

```cpp
std::unique_ptr<VolumetricMesh> volumetricMesh;
```

这一步之后，下游 setup 尽量都只再看：

- `simMesh`
- `volumetricMesh.get()`

因此这一段可以被看成 shared volumetric runtime 的第一个关键拐点。

## 6. surface mesh 和 embedding `W` 怎样复用

runtime 接下来会：

1. 读取 `surface-mesh`
2. 对 surface 顶点应用同样的 `scale`
3. 组装 `surfaceRestPositions`
4. 构造：

```cpp
BarycentricCoordinates(surfaceMesh.numVertices(),
                       surfaceRestPositions.data(),
                       volumetricMesh.get())
```

5. 生成插值矩阵 `W`

这说明 cubic 路径在 surface embedding 这一层并没有特殊旁路。它完全复用了 `VolumetricMesh` 家族统一接口：

- `containsVertex(...)`
- `computeBarycentricWeights(...)`

从架构角度看，这一步意义很大，因为它说明 cubic 在 runtime 中已经不只是“内部弹性模型能算”，而是：

> 边界表面和内部自由度之间的映射也已经进入 shared path。

## 7. `SimulationMesh -> DMM -> Assembler -> Energy` 在 runtime 里怎样落地

前一阶段已经解释过这条装配链的角色分工。放到 `runSimCore` 里，它具体会落成下面这些对象：

1. `DeformationModelManager`
2. `DeformationModelAssembler`
3. `DeformationModelEnergy`
4. `plasticity`
5. `restPosition`
6. 初始 stiffness `K`

其中 `plasticity` 初始化的方式也值得注意：

- runtime 会遍历所有 element
- 从每个 element 的 plastic model 中拿到参数维度
- 用单位矩阵 `I` 初始化 plastic deformation gradient 语义

这说明 cubic 路径到这里时，已经真正进入了“求解前的完整 FEM setup”，而不是停在对象构造演示层。

## 8. 质量矩阵 `M` 怎样进入这条链

当前 runtime 用的是通用接口：

```cpp
GenerateMassMatrix::computeMassMatrix(volumetricMesh.get(), M, true);
```

这里最值得强调的不是 `M` 被算了出来，而是：

- runtime 不需要知道自己拿到的是 tet 还是 cubic
- 质量矩阵层消费的是 `VolumetricMesh*`，不是某个具体子类专用入口

这正是 shared volumetric runtime 的核心设计。

## 9. fixed vertices、外物体和接触在 cubic 路径里怎样组织

当前 runtime 在算完 `K` 之后，会继续处理：

- `fixed-vertices`
- `external-objects`
- 重力与外力
- dynamic branch 下的 contact handler

其中：

- `fixed-vertices` 通过文件读取顶点 id，再构造 `MultipleVertexPulling`
- `external-objects` 通过表面 mesh 和运动轨迹组织成 kinematic objects

一旦进入 dynamic 分支，runtime 还会根据 contact 配置创建：

- `TriangleMeshExternalContactHandler`
- `TriangleMeshSelfContactHandler`

这说明 cubic 路径现在已经能进入：

- fixed vertex pulling
- 外物体接触
- self contact

这些更高层的运行语义，而不只是内部弹性部分。

## 10. contact / IPC 当前已经推进到哪里

当前 `RunSimContactConfig` 已经支持：

- `penalty`
- `ipc-barrier`

`runSimCore` 里也已经按 contact model 组织：

- 是否启用 dynamic contact
- external IPC 参数
- self IPC 参数
- feasible line search 开关

因此对 cubic 路径最准确的表述不是“以后也许能接 IPC”，而是：

> cubic 已经能进入 runtime 的 IPC / penalty contact 入口层。

是否所有长时场景都已经过全面回归，是另一回事；但入口层和例子层已经在仓库里存在。

## 11. 当前 repo 的 example 场景分别说明什么

### 11.1 `examples/cubic-box/`

这是最基础的 cubic 盒子下落场景：

- volumetric mesh: `cubic-box.veg`
- surface mesh: `cubic-box.obj`
- 动力学
- 外部平面障碍物
- penalty contact

它说明 cubic shared volumetric runtime 的最小动态场景已经在 repo 中被组织出来。

### 11.2 `examples/cubic-box-ipc/`

这组例子把 contact model 推进到 IPC barrier。它说明：

- 同一组 cubic 资产可以在 runtime 中切换 contact mode
- cubic 路径已经能进入 barrier-based contact setup

### 11.3 `examples/pulled-cubic-box-self-ipc/`

这组例子进一步把 runtime 场景推进到：

- fixed vertex pulling
- self contact
- IPC feasible line search

它说明 cubic 路径已经能与约束、self contact 和 barrier 参数一起工作。

### 11.4 `examples/dragon-cubic/`

这个例子的重点不是 contact，而是资产来源：

- 它不是规则盒子
- 它来自 `triangle-mesh` 体素化路径

也就是说，Phase 01 的 triangle-mesh asset generation 和 Phase 06 的 runtime consumption 在这里真正闭环。

## 12. 当前 repo 的现实边界应该怎样写

这里需要非常明确，因为 walkthrough 必须区分“仓库里已经存在的能力”和“本次会话重新验证过的能力”。

当前可以稳定声称的是：

- cubic runtime 主链已经在源码中存在
- 多组 cubic example 资产已经 committed
- cubic 已经能进入 penalty / IPC / self-contact 相关 runtime 入口

当前不能在这份文档里声称的是：

- 本次会话重新跑通了所有 example
- 当前环境里已经重新构建了完整 docs 站点
- 所有长时 dynamic/contact 场景都已经在这里重新做过端到端验证

因此本阶段的结论应该写成：

> 仓库里已经存在完整的 cubic runtime 主链和多组 committed example，但这份 walkthrough 记录的是 repo truth，不是一次新的全量运行报告。

## 13. 与现有 walkthrough 的关系

如果你已经读过 [Simulation Tick: Pre-Tick Setup](../simulation-tick-pre-tick.md)，可以把本阶段理解为它的 cubic 视角补充说明：

- 那篇文档按 runtime 时间顺序讲 shared setup
- 本系列则从 cubic 路径视角解释它是怎样被补进这条 shared setup 的

上一阶段： [Phase 05](phase05_energy_assembly.md)  
返回总览： [Cubic Mesh Walkthrough](index.md)
