# Phase 02：配置入口与 Runtime 分发

本阶段把 cubic 从"资产存在"推进到"能从 JSON 配置进入 shared runtime"。

## `RunSimConfig` 结构总览

`runSimCore.h` 定义了一组 config struct，`parseRunSimConfig(...)` 负责从 JSON 填充它们：

```cpp
struct RunSimConfig {
    RunSimMeshConfig       mesh;        // 网格类型、文件路径、scale
    RunSimSceneConfig      scene;       // 重力、初速度、固定顶点、外部物体
    RunSimSimulationConfig simulation;  // timestep、帧数、dump 间隔、sim-type
    RunSimContactConfig    contact;     // contact model、IPC 参数
    RunSimSolverConfig     solver;      // 收敛阈值、最大迭代、弹性材料类型
    RunSimRuntimeConfig    runtime;     // 输出目录、deterministic 模式
};
```

cubic 直接相关的是 `RunSimMeshConfig`：

```cpp
struct RunSimMeshConfig {
    std::string             tetMeshFilename;
    std::string             cubicMeshFilename;
    VolumetricMeshInputType volumetricMeshType = VolumetricMeshInputType::TET;
    std::string             surfaceMeshFilename;
    double                  scale = 1.0;
};

enum class VolumetricMeshInputType { TET, CUBIC };
```

体网格类型通过显式 enum 表达，不靠文件扩展名推断。

## `parseRunSimConfig(...)` 逻辑

### Mesh 类型判定

```text
hasTetMesh   = jconfig.exist("tet-mesh")
hasCubicMesh = jconfig.exist("cubic-mesh")

if both present  → throw（互斥）
if neither       → throw（必须恰好一个）
if hasTetMesh    → volumetricMeshType = TET
if hasCubicMesh  → volumetricMeshType = CUBIC
```

可选字段 `volumetric-mesh-type` 是一致性提示：如果存在且与实际 key 矛盾（比如写了 `"cubic-mesh"` 但 hint 是 `"tet"`），直接报错。冗余一致时静默通过。

### 其他字段解析

所有文件路径通过 `ConfigPathResolver(configFilePath)` 做相对路径解析。其余字段按类型直接读取：

| Config 分组 | 关键字段 | 说明 |
| --- | --- | --- |
| `mesh` | `surface-mesh`, `scale` | surface mesh 独立于 volumetric mesh，scale 统一应用到两者 |
| `scene` | `g`, `init-vel`, `init-disp` | 三维向量 |
| `scene` | `fixed-vertices`, `external-objects` | 数组，每项含 filename + movement |
| `simulation` | `timestep`, `num-timestep`, `sim-type` | `sim-type` 决定 dynamic / static 分支 |
| `contact` | `contact-model`, `contact-stiffness`, ... | `"penalty"` 或 `"ipc-barrier"` |
| `solver` | `solver-eps`, `solver-max-iter`, `elastic-material` | `"stable-neo"` 或 `"stvk-vol"` |
| `runtime` | `output`, `deterministic` | deterministic 模式强制单线程 |

IPC 参数（`external-ipc-dhat` 等）只在 `contact-model = "ipc-barrier"` 时解析，且有正值校验和 friction=0 校验。

## `runSimFromConfig(...)` 分发

这是 shared runtime 的主函数。cubic 在这里的分叉点只有一处 — 读取 mesh 并构造 `SimulationMesh`：

```cpp
if (mesh.volumetricMeshType == VolumetricMeshInputType::TET) {
    auto tetMesh = make_unique<TetMesh>(mesh.tetMeshFilename.c_str());
    for (vi : vertices) tetMesh->setVertex(vi, tetMesh->getVertex(vi) * mesh.scale);
    simMesh        = SimulationMesh::createFromTetMesh(*tetMesh);
    volumetricMesh = std::move(tetMesh);
} else {
    auto cubicMesh = make_unique<CubicMesh>(mesh.cubicMeshFilename.c_str());
    for (vi : vertices) cubicMesh->setVertex(vi, cubicMesh->getVertex(vi) * mesh.scale);
    simMesh        = SimulationMesh::createFromCubicMesh(*cubicMesh);
    volumetricMesh = std::move(cubicMesh);
}
```

分叉之后，两条路径汇合到同一条主线：

```text
volumetricMesh (持有 concrete mesh，供 mass matrix、barycentric 查询)
simMesh        (solver-side 统一容器，供 DeformationModelManager)
    ↓
surface mesh 加载 + scale
    ↓
BarycentricCoordinates(surfaceMesh, volumetricMesh)  → 插值矩阵 W
    ↓
DeformationModelManager → Assembler → DeformationModelEnergy
    ↓
mass matrix M = computeMassMatrix(volumetricMesh)
    ↓
fixed vertices / external objects / contact handler
    ↓
时间积分器 (dynamic) 或静态求解 (static)
```

关键点：

- **scale 在 concrete mesh 上做**：先缩放顶点，再构造 `SimulationMesh`，保证 rest position 和 surface mesh 在同一坐标系
- **`volumetricMesh` 保留到 runtime 结束**：`BarycentricCoordinates` 和 `computeMassMatrix` 都需要原始 `VolumetricMesh` 接口
- **分叉之后无 tet/cubic 判断**：`simMesh->getElementType()` 返回 `CUBIC`，由 `DeformationModelManager` 在 init 时自动选择 `CubicMeshDeformationModel`，见 [Phase 05](phase05_energy_assembly.md)

## 验证

`runSimConfig_parse_test.cpp` 覆盖：cubic 入口选择正确类型、tet+cubic 互斥报错、缺少 mesh key 报错、type hint 不一致报错、hint 一致时通过。

---

上一阶段：[Phase 01](phase01_mesh_assets.md)
下一阶段：[Phase 03](phase03_simulation_mesh.md)
