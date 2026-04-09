# Phase 05：`SimulationMesh -> DMM -> Assembler -> Energy` 装配链

## 1. 这阶段要回答什么

这一阶段专门回答下面三个问题：

1. cubic 单元模型是怎样真正接进全局 FEM 装配链的？
2. `DeformationModelManager`、`DeformationModelAssembler`、`DeformationModelEnergy` 三层各自负责什么？
3. 当前 repo 的 cubic 路径究竟只停在 manager dispatch，还是已经进入全局 energy / gradient / Hessian 装配？

如果这一阶段只写成“DMM 新增了一个 `CUBIC` 分支”，会明显低估当前 repo 的完成度。

## 2. 这阶段在整条主线里解决什么

Phase 04 已经给出了 cubic 单元级 deformation model。  
Phase 05 要解决的则是：

> 单元级结果怎样进入全局装配链，并变成 solver 认识的势能对象。

这一步是从“有局部 element model”到“有全局 FEM energy”的关键跨越。

## 3. 先把 4 个核心对象分清楚

### 3.1 `SimulationMesh`

它是 solver-side 的几何和材料语义容器，负责提供：

- 顶点
- element connectivity
- element type
- per-element material binding

它本身不计算 energy。

### 3.2 `DeformationModelManager`

它不是全局装配器，而是：

> 逐单元 deformation model 工厂 + 注册表

它根据 mesh type、材料绑定和 plastic model 选择，为每个 element 生产一个具体 `DeformationModel` 对象。

### 3.3 `DeformationModelAssembler`

它才是真正遍历所有单元、把局部结果加回全局向量和稀疏矩阵的人。

### 3.4 `DeformationModelEnergy`

它把 assembler 再包装成统一 `PotentialEnergy` 风格的接口，让求解器能直接调：

- `func`
- `gradient`
- `hessian`

## 4. 当前 repo 已经完成了什么

当前仓库里，这条链已经完整存在：

```text
SimulationMesh
  -> DeformationModelManager
  -> DeformationModelAssembler
  -> DeformationModelEnergy
  -> solver / integrator
```

而且 cubic 已经不是只停在 `DeformationModelManager` 的一个 mesh-type dispatch 分支上。它已经能继续被：

- assembler
- energy wrapper
- runtime setup

所消费。

## 5. `DeformationModelManager::init(...)` 现在怎样处理 cubic

### 5.1 它的逻辑分层

当前 `init(...)` 的主体逻辑可以分成三层：

1. 根据 `SimulationMesh` 的 primary / secondary material 选择 elastic material
2. 根据 plastic model 枚举构造 plastic model
3. 根据 `SimulationMeshType` 选择具体 `DeformationModel` 子类

这个分层很重要，因为它说明：

- mesh type dispatch 只是最后一层
- 前面的材料选择逻辑本身已经尽量 mesh-type 无关

### 5.2 cubic 在最后一层怎样接入

对 cubic 路径来说，最后一层现在已经明确存在：

- `SimulationMeshType::CUBIC`
- 收集 8 个局部 rest 顶点
- 组成 `restPosition[24]`
- 构造 `new CubicMeshDeformationModel(restPosition.data(), elementMaterial, pm)`

这一步的关键意义不是“终于不 throw 了”，而是：

- cubic 复用了现有 volumetric material / plastic model 选择逻辑
- 唯一的新增只发生在 mesh-type dispatch 最后一层

这说明 cubic 接入被限制在 mesh-type dispatch 的最小必要位置。

## 6. 为什么说 manager 不是“装配器”

这里非常值得再强调一次，因为 `DeformationModelManager` 这个名字很容易让人误解。

它不做：

- 全局能量求和
- 全局梯度装配
- 全局 Hessian 装配

它做的是：

- 为每个 element 创建一个能独立计算局部 energy / force / stiffness 的对象
- 把这些对象保存起来，供 assembler 之后复用

一句话说：

> `Manager` 解决“为这个 element 创建哪个局部模型”，`Assembler` 才解决“把局部量怎样高效加回全局”。

## 7. `DeformationModelAssembler` 在 cubic 路径里到底做了什么

当前 assembler 的核心职责可以分成两部分：

### 7.1 构造期

构造时它会：

- 读取 mesh 规模
- 保存 element deformation models
- 为每个单元分配 `CacheData`
- 预构造全局 Hessian 稀疏模板
- 为局部到全局的 scatter 建立索引映射

### 7.2 查询期

在 `computeEnergy(...)`、`computeGradient(...)`、`computeHessian(...)` 里，它会：

1. 从全局位置向量里 gather 当前 element 的局部顶点位置
2. gather plastic / elastic 参数
3. 调用 element `prepareData(...)`
4. 调用 element `computeEnergy / compute_dE_dx / compute_d2E_dx2`
5. 把局部结果 scatter 回全局向量或稀疏矩阵

### 7.3 cubic 路径的关键点

对 cubic 路径来说，当前实现已经按 24 个局部自由度组织：

- 局部位置：`V24d localp`
- 局部梯度：`V24d localGradx`
- 局部刚度：`24 x 24` 缓冲

这说明当前 repo 已经明显超出“只是在 DMM 里加一个分支”的最小状态，而是：

> assembler 已经能以 cubic element 的自由度规模做真实全局装配。

## 8. `DeformationModelEnergy` 在 cubic 路径里扮演什么角色

`DeformationModelEnergy` 不区分 tet 或 cubic。它做的事情非常工程化，但很重要：

1. 如果给了 `restPosition`，把优化变量解释成位移 `u`
2. 还原当前绝对坐标 `x = X_rest + u`
3. 调 assembler 的 `computeEnergy / computeGradient / computeHessian`

这一步的意义在于：

- element-level mathematics 到这里被收口成统一的势能接口
- 后面 solver 不需要知道底下是 tet 还是 cubic

因此 cubic 一旦进入 `DeformationModelEnergy`，就说明它已经被真正接到统一求解接口上了。

## 9. 当前 repo 到底完成到哪一步

这里需要非常明确写一句：

当前 repo 的 cubic 路径**不是只接到 manager 层**。

它已经至少接到：

- `SimulationMesh`
- `DeformationModelManager`
- `DeformationModelAssembler`
- `DeformationModelEnergy`

因此如果今天仍把 cubic 描述成“还停在 manager 之前”，就已经不对了。

## 10. 当前 repo 的验证证据

`tests/core/energy/simulationMesh_material_binding_test.cpp` 已经给这条链提供了 committed smoke coverage：

- `DeformationModelManagerCubicStableNeoSmokeTest`
  - 证明 `DeformationModelManager` 已经能为 cubic element 创建 `CubicMeshDeformationModel`
- `CubicAssemblerEnergySmokeTest`
  - 证明 cubic 路径已经能经过：
    - `SimulationMesh::createFromCubicMesh(...)`
    - `DeformationModelManager`
    - `DeformationModelAssembler`
    - `DeformationModelEnergy`
  - 并且能得到：
    - 有限的 energy
    - 接近零的 rest-state gradient
    - 非空 Hessian

这些测试的意义在于：

- 它们验证的不只是某个单元对象能否构造成功
- 而是验证 cubic 已经能穿过整条 global energy assembly 主链

## 11. 这一阶段的边界

这一阶段仍然不深入展开：

- `runSimCore` 里怎样读取 surface mesh
- 怎样构造 embedding `W`
- 怎样构造 mass matrix `M`
- dynamic/contact/IPC 在 runtime 里怎样继续使用这条链

这些属于下一阶段的 full runtime 话题。

上一阶段： [Phase 04](phase04_element_model.md)  
下一阶段： [Phase 06](phase06_full_runtime_examples.md)
