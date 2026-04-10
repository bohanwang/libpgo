# Phase 03：`CubicMesh → SimulationMesh` Handoff

## 为什么需要 `SimulationMesh`

`CubicMesh`（scene/asset 层）和 `SimulationMesh`（solver 层）虽然都描述"单元、顶点和材料"，但语义不同：

| | `CubicMesh` | `SimulationMesh` |
| --- | --- | --- |
| 职责 | 体网格几何、`VolumetricMesh` 家族内材料、插值 | solver-side 统一容器：typed element、per-element binding、统一访问接口 |
| 消费者 | 资产工具、barycentric 查询 | `DeformationModelManager`、`Assembler`、`Energy`、runtime |

Handoff 之后，下游模块不需要关心原始输入是 `CubicMesh` 还是 `TetMesh`。

## 成员布局

`SimulationMesh` 直接持有五个核心成员，没有间接层：

```cpp
std::vector<EigenSupport::V3d>                      vertices_;
std::vector<std::vector<int>>                       elements_;
std::vector<ElementMaterialBinding>                 elementMaterialBindings_;
std::vector<std::unique_ptr<SimulationMeshMaterial>> materials_;
SimulationMeshType                                  meshType_;
```

### `vertices_`

`std::vector<V3d>`，全局顶点位置表。下标就是全局 vertex id，所有 element 的连通性都引用这张表。

### `elements_`

`std::vector<std::vector<int>>`，每个 element 是一个 vertex index 列表。长度由 `meshType_` 决定：TET=4、CUBIC=8、TRIANGLE=3、EDGE_QUAD=4、SHELL=6。

虽然内层用了 `std::vector<int>` 而不是固定大小的 `std::array`，但 `createTyped<N>(...)` 在构造时把 typed `std::array<int, N>` 通过 `memcpy` 搬入，运行时通过 `getNumElementVertices()` 取第一个 element 的 size 来推断 arity。

### `elementMaterialBindings_`

```cpp
struct ElementMaterialBinding {
    int primary = -1;
    int secondary = -1;
};
```

每个 element 一条 binding，索引到 `materials_` 表。`primary` 必须有效，`secondary` 为 `-1` 表示不使用（当前 cubic/tet 路径都是 primary-only）。Shell 类型可能使用 secondary。

下游通过 `getPrimaryMaterial(ele)` / `getSecondaryMaterial(ele)` 访问，内部做 bounds 校验。

### `materials_`

`std::vector<unique_ptr<SimulationMeshMaterial>>`，拥有所有材料实例的所有权。

构造时对传入的每个 `SimulationMeshMaterial*` 调用 `clone()` 做深拷贝，所以 `SimulationMesh` 不依赖外部材料的生命期。`setMaterial(matID, mat)` 也是 clone 语义。

常见子类是 `SimulationMeshENuMaterial`（存 Young's modulus $E$ 和 Poisson's ratio $\nu$），cubic/tet 路径都用它。

### `meshType_`

枚举值 `{TET, CUBIC, TRIANGLE, EDGE_QUAD, SHELL}`，由 typed factory 在构造时写入。`DeformationModelManager` 用它做 switch 分发到不同的 element model。

### 额外成员：`elementUVs_`

`std::vector<std::vector<V2d>>`，per-element per-vertex 的 UV 坐标。通过 `assignElementUVs(...)` 后置赋值，不是构造时必须的。主要用于 shell/triangle 类型。

### Typed Factory

公共构造走 typed static factory：

```text
createTet(...)    → createTyped<4>(TET, ...)
createCubic(...)  → createTyped<8>(CUBIC, ...)
createTriangle(...) → createTyped<3>(TRIANGLE, ...)
...
```

`createTyped<N>(...)` 做三件事：

1. 校验 `N` 与 `meshType` 匹配、element/binding 数量一致
2. 深拷贝 vertices、elements、materials、bindings
3. 校验所有 binding 的 material index 在合法范围内

构造函数是 `protected`，只有 factory 能创建实例，保证 type-arity 一致性。

## `createFromCubicMesh(...)` 做了什么

三步：

### 1. 复制顶点

遍历 `cubicMesh.getNumVertices()`，把顶点位置搬入 `std::vector<Vertex>`。

### 2. 复制 8 点单元连通性

遍历 `cubicMesh.getNumElements()`，收集 8 个 vertex index。**严格沿用 `CubicMesh` 的局部顶点顺序** — 不重新发明编号，因为 `CubicMeshDeformationModel` 的 shape function 依赖同一套顺序。

### 3. 材料 handoff

与 tet 路径策略一致：

1. 从 `CubicMesh` 读出每个 element 的 `ENuMaterial`
2. 转成 solver 层的 `SimulationMeshENuMaterial`
3. 每个 element 绑定一个 primary material（`secondary = -1`）
4. `density` 留在 volumetric / mass matrix 语义消费，不搬入 `SimulationMesh`

最终调用 `createCubic(...)` → `createTyped<8>(SimulationMeshType::CUBIC, ...)`。

## Tet vs Cubic 对齐程度

| 环节 | Tet | Cubic |
| --- | --- | --- |
| solver-side type | `TET` | `CUBIC` |
| typed factory | `createTet(...)` | `createCubic(...)` |
| element arity | 4 | 8 |
| material handoff | `ENuMaterial → SimulationMeshENuMaterial` | 同左 |
| binding | primary-only | primary-only |

两者在 solver-side mesh 层获得了对称的正式表示。

## 验证

`simulationMesh_material_binding_test.cpp` 覆盖：

- `CreateFromCubicMeshBuildsPrimaryOnlyBinding` — type 为 CUBIC、arity 为 8
- `CreateFromCubicMeshCopiesConnectivityAndMaterialParameters` — 顶点/连通性/E/nu 正确转移

---

上一阶段：[Phase 02](phase02_config_runtime_entry.md)
下一阶段：[Phase 04](phase04_element_model.md)
