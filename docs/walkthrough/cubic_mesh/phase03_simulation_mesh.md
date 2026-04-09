# Phase 03：`CubicMesh -> SimulationMesh` Handoff

## 1. 这阶段要回答什么

这一阶段专门回答三个问题：

1. 为什么 `CubicMesh` 不能直接交给 deformation model 层使用？
2. `SimulationMesh::createFromCubicMesh(...)` 具体复制了哪些信息，又刻意没有复制哪些信息？
3. 当前 repo 里，cubic 路径在这一层和 tet 路径到底对齐到了什么程度？

这一步如果没讲清楚，后面非常容易把两个层次混在一起：

- scene/asset 层的 `CubicMesh`
- solver 层的 `SimulationMesh`

它们都在描述“单元、顶点和材料”，但语义并不相同。

## 2. 这阶段在整条主线里解决什么

Phase 02 已经让 cubic 资产能从 config 进入 runtime。  
Phase 03 要解决的是：

> 把 scene-side 的 `CubicMesh`，稳定搬运成 solver-side 的 `SimulationMeshType::CUBIC`。

这一步听上去像“搬数据”，但它实际解决了后面很多层都需要的统一语义：

- solver 看到的 element type 是什么
- solver 怎样访问局部顶点编号
- per-element material 怎样绑定
- volumetric 路径中的材料类型怎样从 `VolumetricMesh` 转成 `SimulationMeshMaterial`

## 3. 当前 repo 已经完成了什么

当前仓库里，这一步所需的核心前置条件都已经在位：

- `SimulationMeshType::CUBIC`
- `SimulationMesh::createCubic(...)`
- `SimulationMesh::createFromCubicMesh(...)`
- `ElementMaterialBinding`
- `getPrimaryMaterial(...)`
- `getSecondaryMaterial(...)`
- typed construction `createTyped<N>(...)`

因此 cubic 在进入 solver 时，已经不再只是一个原始 `.veg` 读入对象，而是一个：

- 带明确 element type
- 带明确 element arity
- 带 per-element material binding
- 带统一访问接口

的 solver-side mesh 容器。

## 4. 为什么这里需要单独的 `SimulationMesh`

### 4.1 `CubicMesh` 的职责仍然偏 scene / volumetric 表示层

`CubicMesh` 负责的是：

- 体网格几何
- 元素拓扑
- `VolumetricMesh` 家族内部的材料表示
- 插值和质量矩阵等 volumetric 通用能力

### 4.2 `SimulationMesh` 的职责是 solver-side 统一容器

`SimulationMesh` 负责的是：

- solver 视角下的 element type
- solver 视角下的 typed element arity
- `SimulationMeshMaterial` 语义
- per-element material binding
- 后续 `DeformationModelManager`、assembler、energy 都能复用的访问接口

也就是说，`SimulationMesh` 不是简单的包装层，而是：

> deformation model 层和 runtime 求解层真正消费的统一网格语义容器。

## 5. `SimulationMesh` 当前的结构是什么

当前 `SimulationMesh` 直接持有：

- `vertices_`
- `elements_`
- `elementUVs_`
- `elementMaterialBindings_`
- `materials_`
- `meshType_`

它暴露的关键查询接口则包括：

- `getNumElements()`
- `getNumVertices()`
- `getNumElementVertices()`
- `getVertexIndex(ele, j)`
- `getVertex(ele, j, ...)`
- `getPrimaryMaterial(ele)`
- `getSecondaryMaterial(ele)`

对于 cubic 路径来说，这组接口的意义在于：

- 后面所有使用 `SimulationMesh` 的模块都不必再关心“原始输入是不是 `CubicMesh`”
- 它们只需要知道当前 element type 是 `CUBIC`

## 6. `createFromCubicMesh(...)` 到底做了什么

### 6.1 复制全局顶点

第一段逻辑是遍历 `cubicMesh.getNumVertices()`，把每个顶点位置复制进：

```cpp
std::vector<Vertex> vertices;
```

这一步不做任何材料处理，也不做任何物理计算。它只是把 scene 层的顶点表搬进 solver 统一容器。

### 6.2 复制 8 点单元连通性

第二段逻辑遍历 `cubicMesh.getNumElements()`，为每个单元收集：

```cpp
{
  cubicMesh.getVertexIndex(ei, 0),
  cubicMesh.getVertexIndex(ei, 1),
  cubicMesh.getVertexIndex(ei, 2),
  cubicMesh.getVertexIndex(ei, 3),
  cubicMesh.getVertexIndex(ei, 4),
  cubicMesh.getVertexIndex(ei, 5),
  cubicMesh.getVertexIndex(ei, 6),
  cubicMesh.getVertexIndex(ei, 7)
}
```

这里最关键的不是“它复制了 8 个 index”，而是：

- 它严格沿用 `CubicMesh` 的局部顶点顺序
- 没有在 solver 层重新发明另一套 hexa 节点编号

这件事对后面的 `CubicMeshDeformationModel` 至关重要，因为单元模型里的 shape function 顺序、积分点几何和 `dFdx` 都依赖同一套局部顶点顺序。

### 6.3 做材料 handoff

这一层的材料策略和 tet 路径保持一致：

1. 从 `CubicMesh` 读出每个 element 的 `VolumetricMesh::ENuMaterial`
2. downcast 成 `ENuMaterial`
3. 转成 solver 层的 `SimulationMeshENuMaterial`
4. 把每个 element 绑定到一个 primary material 槽

这意味着当前 repo 在 volumetric 路径里采用的策略是：

- 每个 element 在 handoff 时先有一个 primary material
- cubic 在这一步不引入 secondary material
- `density` 仍然留在 volumetric / mass matrix 语义里消费

最后这一点尤其值得强调，因为它解释了为什么 `SimulationMesh` 主要承接的是：

- 几何
- element type
- 弹性材料参数 handoff

而不是把所有 volumetric mesh 属性都无差别抄一遍。

### 6.4 建立 `ElementMaterialBinding`

当前实现会为每个 cubic element 创建：

```cpp
{ primary = ei, secondary = -1 }
```

这表示：

- 当前每个 element 对应一个 primary material 实例
- secondary material 暂时未使用

这个策略与 tet 路径对齐，目的不是最节省内存，而是让 solver 语义更直接、更稳定。

## 7. 为什么最终要收口到 `createCubic(...) -> createTyped<8>(...)`

`createFromCubicMesh(...)` 最终不是自己 new 一个对象，而是调用：

```cpp
return createCubic(vertices, elements, elementMaterialBindings, materials);
```

再由 `createCubic(...)` 继续落到：

```cpp
createTyped<8>(SimulationMeshType::CUBIC, ...)
```

这样做的好处有三点：

1. element arity 与 mesh type 的一致性由统一入口检查
2. material binding 数量与 element 数量的匹配由统一入口检查
3. material index 的合法性由统一入口检查

因此这一步的准确描述应当是：

> cubic handoff 不是写一段独立旁路逻辑，而是把 cubic 正式并入 `SimulationMesh` 的 typed construction 体系。

## 8. 当前 repo 已经对齐到了什么程度

如果把 tet 与 cubic 放在这一层对照，可以看到它们已经有相当高的一致性：

| 环节 | Tet 路径 | Cubic 路径 |
| --- | --- | --- |
| solver-side 类型 | `SimulationMeshType::TET` | `SimulationMeshType::CUBIC` |
| typed factory | `createTet(...)` | `createCubic(...)` |
| element arity | 4 | 8 |
| material handoff | `ENuMaterial -> SimulationMeshENuMaterial` | `ENuMaterial -> SimulationMeshENuMaterial` |
| binding 模式 | primary-only | primary-only |

也就是说，这一步并不是“cubic 终于勉强能读进 solver”，而是：

> cubic 已经在 solver-side mesh 这一层获得了和 tet 对称的正式表示。

## 9. 当前 repo 的验证证据

`tests/core/energy/simulationMesh_material_binding_test.cpp` 已经把这一步最关键的行为写成 committed coverage：

- `CreateFromCubicMeshBuildsPrimaryOnlyBinding`
- `CreateFromCubicMeshCopiesConnectivityAndMaterialParameters`

这些测试明确确认：

1. 进入 solver 后的 element type 确实是 `SimulationMeshType::CUBIC`
2. `getNumElementVertices()` 确实返回 8
3. 顶点坐标和 element 连通性没有在 handoff 时被破坏
4. `E / nu` 已经从 `CubicMesh` 的 element material 正确转移到 `SimulationMeshENuMaterial`

## 10. 这一阶段的边界

这一阶段仍然不展开：

- cubic 单元内部怎样计算 `F`
- `CubicMeshDeformationModel` 怎样组织 8 个积分点
- `DeformationModelManager` 怎样选择 element model
- runtime 中怎样构造 `W`、`M`、`K`

它解决的是 `CubicMesh -> SimulationMesh` 这一次关键 handoff，而不是后面所有层。

上一阶段： [Phase 02](phase02_config_runtime_entry.md)  
下一阶段： [Phase 04](phase04_element_model.md)
