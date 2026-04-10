# Cubic Mesh Walkthrough

## 概述

本组文档描述 cubic/hexa FEM 在 `libpgo` 中从资产生成到 runtime 求解的完整链路。阅读起点建议接在 [Simulation Tick: Pre-Tick Setup](../simulation-tick-pre-tick.md) 之后 — 那篇按 runtime 执行顺序讲 setup，本系列从 cubic 路径视角补充。

Cubic 不是独立的 solver 分支，而是通过最少必要分叉接入 shared volumetric runtime。

## 数据流

从资产到 runtime 求解的完整执行链：

```text
cubicMesher (asset generation)
  ↓  体素化 → flood fill → surface extraction
  ↓  输出 CubicMesh (.veg) + surface triangle mesh (.obj)
  ↓
parseRunSimConfig(...)
  ↓  mesh-type = "cubic" → 加载 CubicMesh
  ↓  读取材料参数、timestep、contact 配置
  ↓
SimulationMesh::createFromCubicMesh(cubicMesh, surfaceMesh)
  ↓  枚举 8 节点 hex element、3D embedding、surface 映射
  ↓  产出 SimulationMesh（type = CUBIC）
  ↓
DeformationModelManager::createForCubicMesh(simulationMesh)
  ↓  为每个 hex element 创建 CubicMeshDeformationModel
  ↓  每个 element：8 节点、24 DOF、2×2×2 Gauss 积分点
  ↓
DeformationModelAssembler
  ↓  element energy/gradient/Hessian → 全局 DOF 装配
  ↓  稀疏模板预建（element DOF mapping → triplet）
  ↓
DeformationModelEnergy (PotentialEnergy 接口)
  ↓  func() / gradient() / hessian() → 调用 Assembler
  ↓
runSimFromConfig(...)
  ↓  addGeneralImplicitForceModel(deformationEnergy)
  ↓  + gravity + contact (if enabled) + boundary conditions
  ↓
ImplicitBackwardEuler::doTimestep(...)
  ↓  Newton iteration：deformation Hessian + barrier + ...
  ↓  收敛 → 更新 position/velocity → 下一帧
```

## Phase Map

| Phase | 主题 | 主要入口 |
| --- | --- | --- |
| [Phase 01](phase01_mesh_assets.md) | 网格表示与资产生成 | `CubicMesh`, `cubicMesher` |
| [Phase 02](phase02_config_runtime_entry.md) | 配置入口与 shared runtime 前置整理 | `RunSimMeshConfig`, `parseRunSimConfig` |
| [Phase 03](phase03_simulation_mesh.md) | `CubicMesh → SimulationMesh` handoff | `SimulationMesh::createFromCubicMesh` |
| [Phase 04](phase04_element_model.md) | 单元模型：energy / gradient / Hessian | `CubicMeshDeformationModel` |
| [Phase 05](phase05_energy_assembly.md) | 全局装配链 | `DeformationModelManager`, `Assembler`, `Energy` |
| [Phase 06](phase06_full_runtime_examples.md) | Full runtime 与 example 场景 | `runSimFromConfig`, example 目录 |

## 关键代码入口

**资产层：**

- `src/core/scene/volumetricMesh/cubicMesh.{h,cpp}`
- `src/tools/cubicMesher/cubicMesher.cpp`

**Solver-side mesh：**

- `src/core/energy/solidDeformationModel/simulationMesh.{h,cpp}`
- `src/core/energy/solidDeformationModel/cubicMeshDeformationModel.{h,cpp}`

**装配链：**

- `src/core/energy/solidDeformationModel/deformationModelManager.cpp`
- `src/core/energy/solidDeformationModel/deformationModelAssembler.cpp`
- `src/core/energy/solidDeformationModel/deformationModelEnergy.cpp`

**Runtime：**

- `src/api/runSimCore.{h,cpp}`
