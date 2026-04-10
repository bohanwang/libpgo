# Cubic Mesh Walkthrough

这组 walkthrough 的目标是把 repo 里分散在源码、测试和 example 资产中的 cubic/hexa 线索，整理成一条**对齐当前 repo 现状**的完整主线。

阅读时需要先有一个总体判断：

- 当前仓库里的 cubic FEM 已经不是停留在表示层的草稿。
- 它已经覆盖资产生成、配置入口、solver handoff、单元模型、装配链和 runtime example。
- 因此这里不按“待实现步骤”组织，而是按“现在仓库已经完成了哪几个阶段”组织。

如果你已经读过 [Simulation Tick: Pre-Tick Setup](../simulation-tick-pre-tick.md)，可以把本系列看成它的补充视角：

- `simulation-tick-pre-tick.md` 站在 runtime 视角，按 `runSimCore` 的执行顺序讲 setup。
- 本系列站在 cubic 路径视角，解释 cubic 是怎样一步一步被补进这条 shared volumetric runtime 的。

## 这组文档试图回答什么

本系列集中回答六个问题：

1. libpgo 里的 cubic/hexa 体网格到底是什么离散对象，它和 tet 的关系是什么？
2. repo 现在怎样生成 cubic 资产，而不是只假设外部已经准备好 `.veg`？
3. cubic mesh 是怎样从 JSON 配置进入 `runSimCore` 的？
4. `CubicMesh -> SimulationMesh -> CubicMeshDeformationModel` 这三层 handoff 分别解决什么问题？
5. cubic 路径现在究竟只到了 manager，还是已经接进 assembler、energy 和 runtime？
6. 现有 example 资产说明了 repo 已经支持哪些 cubic 运行场景，边界又在哪里？

## 为什么按 phase 组织

当前 cubic 路径跨越了资产层、solver-side mesh、单元模型、装配链和 runtime setup。

如果直接顺着代码目录或函数调用临时跳读，很容易把几类语义混在一起：

- asset generation
- solver-side handoff
- element-level mathematics
- full runtime setup

因此这组文档采用 phase 结构，把这些层按依赖顺序拆开，避免把“当前 repo 已完成的事实”和“某一层局部实现细节”混成一团。

## 阶段总览

| 阶段 | 主题 | 当前 repo 状态 | 入口文档 |
| --- | --- | --- | --- |
| Phase 01 | 网格表示与资产生成 | 已完成，覆盖规则体块和 triangle-mesh 两条资产生成路径 | [Phase 01](phase01_mesh_assets.md) |
| Phase 02 | 配置入口与 shared runtime 前置整理 | 已完成 | [Phase 02](phase02_config_runtime_entry.md) |
| Phase 03 | `CubicMesh -> SimulationMesh` handoff | 已完成 | [Phase 03](phase03_simulation_mesh.md) |
| Phase 04 | `CubicMeshDeformationModel` 单元模型 | 已完成 | [Phase 04](phase04_element_model.md) |
| Phase 05 | `SimulationMesh -> DMM -> Assembler -> Energy` 装配链 | 已完成 | [Phase 05](phase05_energy_assembly.md) |
| Phase 06 | `runSimCore` shared volumetric runtime 与例子场景 | 已完成，但运行层边界需要单独说明 | [Phase 06](phase06_full_runtime_examples.md) |

## 当前 repo 状态应该怎样理解

如果只看零散类名或旧印象，很容易以为 cubic 还停在“表示层已经有了，但 handoff 都还没接”的阶段。对当前仓库来说，这个判断已经过时了。

现在仓库里已经有：

- `CubicMesh` 表示层
- `cubicMesher` 工具，并且有 `uniform` 和 `triangle-mesh` 两种资产生成路径
- `RunSimConfig` 的 cubic 配置入口
- `SimulationMesh::createFromCubicMesh(...)`
- `CubicMeshDeformationModel`
- `DeformationModelManager` 的 `CUBIC` 分支
- `DeformationModelAssembler` / `DeformationModelEnergy` 的 cubic smoke 覆盖
- `runSimCore` 里的 cubic shared volumetric runtime 路径
- 多个 cubic example 目录

因此本系列的主线不再是“怎样规划 cubic FEM”，而是：

> 当前 repo 里的 cubic FEM 已经怎样落地，以及每个阶段到底完成了哪些真实工作。

## 关键入口文件

如果只想抓主干，优先看下面这些入口：

- `src/core/scene/volumetricMesh/cubicMesh.h`
- `src/core/scene/volumetricMesh/cubicMesh.cpp`
- `src/tools/cubicMesher/cubicMesher.cpp`
- `src/api/runSimCore.h`
- `src/api/runSimCore.cpp`
- `src/core/energy/solidDeformationModel/simulationMesh.h`
- `src/core/energy/solidDeformationModel/simulationMesh.cpp`
- `src/core/energy/solidDeformationModel/cubicMeshDeformationModel.h`
- `src/core/energy/solidDeformationModel/cubicMeshDeformationModel.cpp`
- `src/core/energy/solidDeformationModel/deformationModelManager.cpp`
- `src/core/energy/solidDeformationModel/deformationModelAssembler.cpp`
- `src/core/energy/solidDeformationModel/deformationModelEnergy.cpp`

对应的验证入口则主要是：

- `tests/tools/cubicMesher_test.cpp`
- `tests/api/runSimConfig_parse_test.cpp`
- `tests/core/energy/simulationMesh_material_binding_test.cpp`
- `tests/core/energy/cubicMeshDeformationModel_test.cpp`

## 建议阅读顺序

如果你是第一次系统看这条线，建议严格按下面顺序读：

1. [Phase 01: 网格表示与资产生成](phase01_mesh_assets.md)
2. [Phase 02: 配置入口与 shared runtime 前置整理](phase02_config_runtime_entry.md)
3. [Phase 03: `CubicMesh -> SimulationMesh` handoff](phase03_simulation_mesh.md)
4. [Phase 04: `CubicMeshDeformationModel` 单元模型](phase04_element_model.md)
5. [Phase 05: 能量装配链](phase05_energy_assembly.md)
6. [Phase 06: full runtime 与 examples](phase06_full_runtime_examples.md)

这个顺序和代码里的实际依赖关系一致：

- 没有 asset，就没有 config。
- 没有 config 和 solver-side mesh 语义整理，runtime 很难共享。
- 没有 `SimulationMesh` handoff，就没有 DMM。
- 没有 element model，就没有 assembler / energy。
- 没有 runtime setup，就没法解释 example 目录为什么有意义。

## 这组文档的写作约束

- 只使用 docs 站点内可用的相对链接。
- 只使用 repo-relative 代码路径，不暴露本地绝对路径。
- 文中提到的“已完成”都以当前 committed 源码、测试文件和 example 资产为准。
- 文中可以引用 committed coverage，但不会声称“本次会话已重新跑通全部测试或 example”。

## 一个最简 mental model

如果只想先记一句话，可以先记下面这条：

```text
asset generation
  -> config/runtime entry
  -> SimulationMesh
  -> CubicMeshDeformationModel
  -> DMM / Assembler / Energy
  -> runSimCore shared volumetric runtime
  -> cubic examples
```

接下来的 6 个 phase，就是把这条链逐段拆开。
