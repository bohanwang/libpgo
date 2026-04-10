# Phase 06：Full Runtime 与 Example 场景

## `runSimFromConfig(...)` 中的 Cubic 路径

Cubic 进入 shared runtime 后的完整 setup 序列：

```text
1. 读取 CubicMesh (.veg)
2. 统一 scale
3. createFromCubicMesh(...) → SimulationMesh
4. move 到 unique_ptr<VolumetricMesh>
   ─── 以下与 tet 共享 ───
5. 读取 surface mesh (.obj)，统一 scale
6. BarycentricCoordinates → embedding 矩阵 W
7. DeformationModelManager → Assembler → Energy
8. 初始化 plasticity 参数
9. 计算初始 stiffness K
10. GenerateMassMatrix::computeMassMatrix(volumetricMesh.get(), M)
11. fixed vertices / external objects / gravity
12. dynamic branch → contact handler (penalty / ipc-barrier)
13. time integration
```

步骤 5 之后，runtime 只看 `simMesh` + `volumetricMesh.get()`，cubic 和 tet 完全收口。

## Surface Embedding `W`

```cpp
BarycentricCoordinates(surfaceMesh.numVertices(),
                       surfaceRestPositions.data(),
                       volumetricMesh.get())
```

复用 `VolumetricMesh` 家族统一接口（`containsVertex`、`computeBarycentricWeights`），无 cubic 专用旁路。

## 质量矩阵 `M`

```cpp
GenerateMassMatrix::computeMassMatrix(volumetricMesh.get(), M, true);
```

消费 `VolumetricMesh*`，不需要知道是 tet 还是 cubic。

## Contact / IPC

Cubic 路径能进入 `RunSimContactConfig` 支持的所有模式（`penalty` / `ipc-barrier`），包括 external contact、self contact、feasible line search。入口层和 tet 共享。

## Example 场景

| 目录 | 场景 | 资产来源 |
| --- | --- | --- |
| `examples/cubic-box/` | 基础下落 + penalty contact | `uniform` |
| `examples/cubic-box-ipc/` | 下落 + IPC barrier contact | `uniform` |
| `examples/pulled-cubic-box-self-ipc/` | fixed vertex pulling + self IPC + feasible line search | `uniform` |
| `examples/dragon-cubic/` | triangle-mesh 体素化产物运行 | `triangle-mesh --resolution 64` |

`dragon-cubic` 闭合了 Phase 01 的 triangle-mesh 资产生成到 Phase 06 的 runtime 消费。

---

上一阶段：[Phase 05](phase05_energy_assembly.md)
返回总览：[Cubic Mesh Walkthrough](index.md)
