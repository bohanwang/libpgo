# Phase 05：全局装配链

从单元级 `CubicMeshDeformationModel` 到 solver 消费的 `PotentialEnergy`，需要经过四个对象：

```text
SimulationMesh → DeformationModelManager → DeformationModelAssembler → DeformationModelEnergy
```

## 四层角色分工

### `SimulationMesh`

Solver-side 几何和材料容器。提供顶点、element connectivity、element type、per-element material binding。**不计算 energy。**

### `DeformationModelManager`

逐单元 deformation model **工厂 + 注册表**（不是装配器）。根据 mesh type、材料绑定和 plastic model，为每个 element 生产一个具体 `DeformationModel` 对象。

Cubic 分支：

```cpp
case SimulationMeshType::CUBIC:
    // 收集 8 个 rest 顶点 → restPosition[24]
    new CubicMeshDeformationModel(restPosition.data(), elementMaterial, pm);
```

Mesh type dispatch 只是最后一层，前面的材料/plastic 选择逻辑已经 mesh-type 无关。

### `DeformationModelAssembler`

真正遍历所有单元、把局部结果 scatter 回全局向量和稀疏矩阵的层。

**构造期：** 分配 per-element `CacheData`、预构造全局 Hessian 稀疏模板、建立 local→global scatter 索引映射。

**查询期（每次 energy/gradient/hessian 求值）：**

1. 从全局位置向量 gather 局部顶点位置
2. Gather plastic/elastic 参数
3. `prepareData(...)`
4. `computeEnergy / compute_dE_dx / compute_d2E_dx2`
5. Scatter 回全局向量或稀疏矩阵

Cubic 路径按 24 个局部自由度组织（`V24d localp`、`V24d localGradx`、`24×24` 刚度缓冲）。

### `DeformationModelEnergy`

把 assembler 包装成统一 `PotentialEnergy` 接口：

- 如果给了 `restPosition`，把优化变量解释成位移 $u$，还原 $x = X_{\text{rest}} + u$
- 调 assembler 的 `computeEnergy / computeGradient / computeHessian`

这一层不区分 tet/cubic — 一旦进入 `DeformationModelEnergy`，solver 不需要知道底下的单元类型。

## 验证

`simulationMesh_material_binding_test.cpp` 覆盖：

- `DeformationModelManagerCubicStableNeoSmokeTest` — DMM 能为 cubic element 创建 `CubicMeshDeformationModel`
- `CubicAssemblerEnergySmokeTest` — cubic 穿过整条装配链，得到有限 energy、接近零的 rest-state gradient、非空 Hessian

---

上一阶段：[Phase 04](phase04_element_model.md)
下一阶段：[Phase 06](phase06_full_runtime_examples.md)
