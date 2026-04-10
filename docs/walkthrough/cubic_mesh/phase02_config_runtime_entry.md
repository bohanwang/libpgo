# Phase 02：配置入口与 Shared Runtime 前置整理

## 1. 这阶段要回答什么

这一阶段专门回答下面三个问题：

1. cubic mesh 是怎样从 JSON 配置进入 `runSimCore` 的？
2. 为什么这一阶段不能只把它看成“加了一个 parser 字段”？
3. 当前 repo 为什么能让 tet 和 cubic 共用同一条 volumetric runtime 主线，而不是把 cubic 做成一条独立旁路？

如果只看 “`parseRunSimConfig(...)` 支持 `cubic-mesh` 了” 这一点，会低估这一阶段的真正意义。它不只是配置层的增量，而是：

> config 语义、solver-side mesh 语义和 runtime handoff 语义同时被整理到足够一致，cubic 才能干净地接进 shared volumetric path。

## 2. 这阶段在整条主线里解决什么

Phase 01 解决的是“cubic 资产如何存在”。  
Phase 02 要解决的则是“这些资产怎样进入仿真入口，而且进入以后语义不再混乱”。

从代码职责上看，这一步包含三块：

- 扩展 `RunSimConfig` 和 JSON 解析
- 整理 `SimulationMesh` 材料绑定接口
- 整理 `SimulationMesh` 成员布局和 typed construction

如果站在“代码提交顺序”上看，它们像是三步。  
但如果站在“cubic runtime 为什么终于能接得顺”这个问题上看，它们其实是在共同完成一件事：

> 把原来 tet-only 的入口语义和 mesh handoff 语义整理成适合同时承载 tet 与 cubic 的形状。

## 3. 当前 repo 已经完成了什么

当前仓库里，这一阶段已经完成的内容包括：

- `RunSimMeshConfig` 同时持有：
  - `tetMeshFilename`
  - `cubicMeshFilename`
  - `volumetricMeshType`
  - `surfaceMeshFilename`
  - `scale`
- `VolumetricMeshInputType` 已经显式区分：
  - `TET`
  - `CUBIC`
- `parseRunSimConfig(...)` 已经支持：
  - `tet-mesh`
  - `cubic-mesh`
  - 可选 `volumetric-mesh-type`
- `SimulationMesh` 当前已经具备：
  - `ElementMaterialBinding`
  - typed factory
  - `getPrimaryMaterial(...)`
  - `getSecondaryMaterial(...)`
  - `SimulationMeshType::CUBIC`
- `runSimFromConfig(...)` 当前已经能在 runtime 里：
  - 分叉读取 tet 或 cubic concrete mesh
  - 在 concrete mesh 上应用统一缩放
  - 立刻创建 solver-side `SimulationMesh`
  - 之后把两者收口到 shared volumetric handle

因此当前 repo 的真实状态应表述为：

- 不是只有 parser 接受了 `cubic-mesh`
- 而是 cubic 配置语义已经真正进入 runtime 主链

## 4. 当前设计里，配置层到底新增了什么

### 4.1 `RunSimMeshConfig` 不再假设只有 tet

`src/api/runSimCore.h` 里的 `RunSimMeshConfig` 当前字段是：

```cpp
std::string             tetMeshFilename;
std::string             cubicMeshFilename;
VolumetricMeshInputType volumetricMeshType = VolumetricMeshInputType::TET;
std::string             surfaceMeshFilename;
double                  scale = 1.0;
```

这个结构本身已经体现出两个重要判断：

1. 体网格输入类型不是通过文件扩展名猜，而是通过显式语义字段表达。
2. surface mesh 不是附属物，而是和 volumetric mesh 并列出现在 mesh config 里。

### 4.2 `VolumetricMeshInputType` 的意义

当前 repo 并没有让后续代码自己猜“到底是 tet 还是 cubic”。相反，parser 一开始就把这一层语义固定成：

```cpp
enum class VolumetricMeshInputType {
    TET,
    CUBIC,
};
```

这件事看起来简单，但它对后面 runtime 结构的影响很大：

- 读取 concrete mesh 的分支可以保持清晰
- shared runtime 重新汇合的位置也能保持清晰
- 测试可以直接断言输入类型，而不是绕过一堆字符串推断

## 5. `parseRunSimConfig(...)` 现在怎样工作

### 5.1 它先解决的是输入类型，不是材料

配置文件进入 `parseRunSimConfig(...)` 之后，最先要确定的不是材料模型，也不是 timestep，而是：

> 这次运行的 volumetric mesh 究竟是哪一种。

当前实现采用的是“显式二选一”策略：

1. 检查 `tet-mesh` 是否存在
2. 检查 `cubic-mesh` 是否存在
3. 如果两者都存在，直接报错
4. 如果两者都不存在，也直接报错
5. 只有一个存在时，解析路径并设置 `VolumetricMeshInputType`

这一步把很多后续歧义提前消掉了。

### 5.2 它为什么还支持 `volumetric-mesh-type`

`volumetric-mesh-type` 在当前设计里不是主入口，而是一个**一致性提示字段**。  
它允许配置显式写出：

- `tet`
- `cubic`

但前提是这个提示值必须和真正使用的 mesh key 一致。也就是说：

- `"cubic-mesh": "..."`
- `"volumetric-mesh-type": "tet"`

这种写法会直接报错。

这说明 repo 当前的 parser 设计并不是“多给一个字段让用户自己填着玩”，而是：

> 允许冗余表达，但绝不允许互相矛盾。

### 5.3 path resolver 在这里的重要性

parser 里所有和文件相关的输入都会通过 `FileService::ConfigPathResolver` 解析。对 cubic 路径来说，这意味着：

- `cubic-mesh`
- `surface-mesh`
- `fixed-vertices` 里的文件
- `external-objects` 里的文件
- `output`

都会相对于配置文件路径做稳定解析。

这个细节对例子目录尤其重要，因为很多 example 都依赖相对路径组织。

## 6. 为什么材料绑定和成员布局不能只当“内部清理”

`SimulationMesh` 里的材料绑定和成员布局看起来像内部重构，但从当前 repo 角度看，它们直接影响 cubic runtime 是否能干净接入。

### 6.1 `ElementMaterialBinding`

如果没有 `ElementMaterialBinding`，后面至少会出现两个问题：

- cubic 路径很难和 tet 路径统一 primary / secondary material 语义
- `DeformationModelManager` 读取材料时会继续依赖更脆弱的隐式约定

当前结构把 per-element material 语义整理成：

```cpp
struct ElementMaterialBinding {
    int primary = -1;
    int secondary = -1;
};
```

这让 volumetric 与 shell 路径都能在统一接口下表达材料绑定。

### 6.2 typed factory

当前 `SimulationMesh` 并不是随手 new 一个对象，而是通过：

- `createTet(...)`
- `createCubic(...)`
- `createTriangle(...)`
- `createEdgeQuad(...)`
- `createShell(...)`

统一落到 `createTyped<N>(...)`。

这一步对 cubic 很重要，因为它意味着：

- cubic 并不是临时 patch 进去的 mesh type
- 而是正式进入了 `SimulationMesh` 的 typed construction 体系

### 6.3 更直接的成员布局

当前 `SimulationMesh` 直接持有：

- `vertices_`
- `elements_`
- `elementMaterialBindings_`
- `materials_`
- `meshType_`

这让后面 runtime 和文档解释都更清晰。  
对当前 walkthrough 来说，这一点也很重要，因为 cubic 路径需要被解释成一条稳定主线，而不是一堆绕来绕去的间接层。

## 7. `runSimFromConfig(...)` 是怎样把两条 concrete path 收口成 shared path 的

当前 runtime setup 的关键结构是：

1. 按 `volumetricMeshType` 分叉读取 `TetMesh` 或 `CubicMesh`
2. 在 concrete mesh 上做统一缩放
3. 分别调用：
   - `SimulationMesh::createFromTetMesh(...)`
   - `SimulationMesh::createFromCubicMesh(...)`
4. 再把 concrete mesh move 到：

```cpp
std::unique_ptr<VolumetricMesh> volumetricMesh;
```

后面的 shared runtime 统一只看：

- `simMesh`
- `volumetricMesh.get()`

这一步正是 shared volumetric runtime 的关键结构，只不过当前 repo 已经把它真正落成代码了。

## 8. 为什么说 cubic 不是“复制一份 tet 主链”

如果 cubic 只是靠复制 tet 主链实现，代码结构通常会变成这样：

- `runTetSim(...)`
- `runCubicSim(...)`
- 大量 setup 逻辑重复

当前 repo 没有走这条路。它采取的是：

- concrete mesh load 分叉
- `SimulationMesh` handoff 分叉
- 之后尽量收敛到统一 `VolumetricMesh*` 和统一 solver-side pipeline

因此更准确的描述应该是：

> cubic 通过最少必要分叉接进 shared runtime，而不是自立门户。

## 9. 当前 repo 的验证证据

`tests/api/runSimConfig_parse_test.cpp` 已经覆盖了这阶段最关键的 parser 语义，包括：

- `CubicMeshConfigSelectsCubicInputType`
- `RejectsConfigWithBothTetAndCubicMeshKeys`
- `RejectsConfigWithMissingVolumetricMeshKey`
- `RejectsVolumetricMeshTypeHintMismatch`
- `AcceptsVolumetricMeshTypeHintWhenConsistent`

这些测试共同固定了下面几个事实：

- cubic 入口不是文档约定，而是 parser 已经显式支持的运行语义
- `volumetric-mesh-type` 的一致性检查不是口头约束，而是 committed 行为
- tet / cubic 的“互斥而且至少有一个”的规则已经被测试固定下来

## 10. 这一阶段的边界

这一阶段仍然不展开：

- `createFromCubicMesh(...)` 的具体 handoff 细节
- cubic 单元如何计算 `energy / gradient / Hessian`
- 全局装配链的内部角色分工
- dynamic/contact/IPC 细节

它解决的是“入口语义”和“shared runtime 的前置整理”，不是后面所有层的数学和 runtime 行为本身。

上一阶段： [Phase 01](phase01_mesh_assets.md)  
下一阶段： [Phase 03](phase03_simulation_mesh.md)
