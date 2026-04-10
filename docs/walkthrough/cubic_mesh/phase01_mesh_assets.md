# Phase 01：网格表示与资产生成

## 1. 这阶段要回答什么

这一阶段专门回答三个基础问题：

1. libpgo 里的 cubic/hexa 体网格究竟是什么对象，它和 tet 的共性、差异分别在哪里？
2. 一个 cubic element 的局部顶点顺序怎样固定下来，为什么这个顺序后面所有层都必须沿用？
3. 当前 repo 怎样生成 cubic 资产，而不是只停留在“如果外部已经有 `.veg` 就继续往下跑”的状态？

如果这三个问题没有先讲清楚，后面所有关于 `SimulationMesh`、`CubicMeshDeformationModel`、`runSimCore` 的讨论都会失焦，因为你会不断把下面几层语义混在一起：

- 资产层的体单元表示
- 求解层的自由度容器
- 表面几何层的 `.obj`
- contact / embedding 需要的边界表示

## 2. 这一阶段在整个主线里解决什么

这一阶段解决的是：

> 把 cubic/hexa 从“概念上存在的体网格类型”，推进成 repo 里可生成、可落盘、可被例子消费的资产层事实。

它还没有进入 solver，但它已经完成了后面所有阶段的前提：

- 没有 `CubicMesh` 的稳定局部顶点顺序，Phase 03 的 handoff 无从谈起。
- 没有 `cubicMesher`，Phase 06 的 `dragon-cubic` 这类例子就只能变成仓库外的手工资产。
- 没有 `.veg + .obj` 的成对资产语义，运行时的内部力学和表面接触也无法分工。

## 3. 当前 repo 已经完成了什么

当前仓库里，Phase 01 已经完成的不只是“有一个 `CubicMesh` 类”，而是下面几层能力都已经存在：

- `CubicMesh` 作为 `VolumetricMesh` 的具体子类，表示规则、等边、轴对齐的 hexa 体单元集合
- `CubicMesh::createFromUniformGrid(...)`，把激活 voxel 索引集转换成共享顶点表和 8 点单元连通性
- `cubicMesher` CLI 工具
- `uniform` 资产生成路径
- `triangle-mesh` 资产生成路径
- 可选导出 surface `.obj`
- 多组 committed cubic 资产目录和例子配置

把当前状态概括成一句话：

> repo 已经不只是“计划做一个 cubic mesher”，而是已经有一个能生成规则体块、也能体素化输入三角网格的资产层工具链。

## 4. `CubicMesh` 本体到底是什么

### 4.1 类的职责

`src/core/scene/volumetricMesh/cubicMesh.h` 里的 `CubicMesh` 继承自 `VolumetricMesh`。它负责的不是任意 hexa mesh，而是非常明确的一类：

- 所有 cube 都边长相等
- 默认与坐标轴对齐
- 拓扑上来自规则网格的 voxel 子集

这意味着它不是一个“任意 8 点 hexa 容器”。它带有很强的网格生成假设，而这恰好也是后面很多简化成立的原因。

### 4.2 局部顶点顺序

当前 repo 里 cubic element 的 8 个局部顶点顺序已经固定为：

```text
0: 000
1: 100
2: 110
3: 010
4: 001
5: 101
6: 111
7: 011
```

源码里对应的构造辅助数组是：

```cpp
int vtxI[8] = {0, 1, 1, 0, 0, 1, 1, 0};
int vtxJ[8] = {0, 0, 1, 1, 0, 0, 1, 1};
int vtxK[8] = {0, 0, 0, 0, 1, 1, 1, 1};
```

这个顺序不是可随便替换的实现细节。后面至少有 4 层逻辑都依赖它：

- `CubicMesh::computeBarycentricWeights(...)`
- `SimulationMesh::createFromCubicMesh(...)`
- `CubicMeshDeformationModel`
- 表面提取和 example 资产的一致性

所以在当前 repo 里，它应该被当成一个**跨层 contract**，而不是某个 `.cpp` 文件里的局部约定。

### 4.3 `createFromUniformGrid(...)` 在做什么

`createFromUniformGrid(...)` 的输入不是 surface mesh，而是激活 voxel 的 `(i, j, k)` 网格索引。它完成的工作可以拆成四步：

1. 遍历每个 voxel，把 8 个角点加入一个共享顶点集合
2. 为这些角点建立 `(i, j, k) -> vertex id` 的映射
3. 按固定局部顶点顺序组装每个 cube 的 8 点连通性
4. 用顶点数组和 element 数组构造最终 `CubicMesh`

这里一个很重要但很容易被忽略的实现细节是：

- 它默认把所有坐标放进规范域 `[-0.5, 0.5]^3`
- 也就是说，它先生成一个“规范 cubic mesh”
- 具体世界尺度、平移和额外用户变换，交给工具层后处理

这个设计让 asset generator 和后续变换逻辑分工更清晰。

## 5. `cubicMesher` 现在到底有哪些入口

### 5.1 `uniform` 子命令

`uniform` 是最直接的 cubic 资产生成入口，但当前实现已经不只是一个最小示例。它支持：

- `--resolution`
- `--output-mesh`
- 可选 `--output-surface`
- `--size`
- `--offset`
- `--E`
- `--nu`
- `--density`

也就是说，当前 `uniform` 路径不只是“生成一个规范域里的 N x N x N 立方体体网格”，它还已经支持：

- 直接控制几何尺寸
- 直接控制中心偏移
- 直接写入材料参数

最小示例：

```bash
cubicMesher uniform \
  --resolution 4 \
  --output-mesh cubic-box.veg \
  --output-surface cubic-box.obj
```

如果再加上：

```bash
--size 2.0 --offset 1.5 -2.0 0.25
```

生成出来的 `.veg` bounding box 就会反映这些几何变换，而不是停留在规范域。

### 5.2 `triangle-mesh` 子命令

这是当前 repo 更通用的一条资产生成入口。它把输入三角网格体素化成 cubic mesh，并支持：

- `--input-mesh`
- `--resolution`
- `--output-mesh`
- 可选 `--output-surface`
- `--padding-voxels`
- `--scale`
- `--offset`
- `--classify-mode`
- `--E`
- `--nu`
- `--density`

它内部做的事情比 `uniform` 明显更复杂：

1. 读取并清洗 `.obj`
2. 检查 manifold、closed、self-intersection 等前提
3. 根据 `resolution` 和 `padding-voxels` 构造 cubic voxelization domain
4. 按 `center` 或 `conservative` 模式做 occupancy 判定
5. 调 `CubicMesh::createFromUniformGrid(...)`
6. 把规范域网格变换回世界坐标域
7. 再应用用户提供的 `scale` / `offset`

这一步的意义不只是“支持另一种输入格式”，而是：

> cubic asset 已经不必局限于规则盒子，而可以从任意满足约束的闭合 triangle mesh 派生。

## 6. `triangle-mesh` 路径解决了什么真实问题

如果只看 `uniform` 这条最基础路径，cubic 最容易给人一种印象：

- 似乎只适合做规则盒子
- 似乎只是为了验证 hexa FEM 能不能跑起来

当前 repo 已经用 `triangle-mesh` 体素化打破了这个限制。它说明 cubic 路径不是只为“教学方块”存在，而是已经能够服务于更通用的资产生产。

`examples/dragon-cubic/dragon-cubic.json` 的注释就明确说明：

- 这组资产来自 `triangle-mesh --resolution 64`

因此从资产链角度看，当前 repo 已经完成了下面这个闭环：

```text
triangle mesh (.obj)
  -> voxelization
  -> cubic volumetric mesh (.veg)
  -> runtime example
```

## 7. 为什么 `.veg` 和 `.obj` 要同时存在

即使 cubic 是体网格，运行时依然会单独消费一个 surface `.obj`。这不是历史包袱，而是职责分工：

- `.veg` 负责内部力学离散、体材料、质量和自由度
- `.obj` 负责显示、接触边界、表面嵌入和边界 primitive 组织

这和当前 runtime 里的职责划分是一致的：体网格负责内部力学，表面网格负责边界几何。

对当前 repo 来说，这个分工已经体现在工具层：

- `cubicMesher` 可选导出 surface `.obj`
- `runSimCore` 总是显式读取 surface mesh

因此不要把 surface `.obj` 误解成“只是为了可视化方便留的一份副本”。在运行时，它仍然承担实际的边界语义。

## 8. 当前 repo 里有哪些 committed 资产

这一阶段不只是工具实现，还已经沉淀出多组 committed example 资产：

- `examples/cubic-box/`
- `examples/cubic-box-ipc/`
- `examples/pulled-cubic-box-self-ipc/`
- `examples/dragon-cubic/`

这些目录的存在本身就说明：

- cubic 资产不是临时生成后手工丢弃
- 它们已经被当作 runtime example 的正式输入

其中：

- `cubic-box*` 系列说明规则体块资产生成已经稳定可用
- `dragon-cubic` 说明 triangle-mesh voxelization 产物已经被运行层消费

## 9. 当前 repo 的验证证据

`tests/tools/cubicMesher_test.cpp` 已经把这阶段最核心的工具行为写成 committed coverage，包括：

- `UniformResolution2ProducesExpectedCubicMesh`
- `UniformResolution4CanExportSurfaceMesh`
- `UniformSizeScalesOutputBoundingBox`
- `UniformOffsetMovesOutputBoundingBoxCenter`
- `TriangleMeshCenterVoxelizationProducesExpectedCubicMesh`

这些测试分别固定了：

- `uniform` 路径的 vertex / element 数量
- surface `.obj` 导出是否成功且可读
- `size` 是否正确影响输出 bounding box
- `offset` 是否正确移动输出中心
- `triangle-mesh` 路径至少有 committed 行为基线

因此这阶段最稳妥的表述不是“仓库里有一个 cubic mesher”，而是：

> 仓库里已经有一个带测试约束的 cubic 资产生成层。

## 10. 这一阶段不覆盖什么

为了保持边界清晰，这一阶段不讨论：

- `SimulationMesh`
- `CubicMeshDeformationModel`
- `DeformationModelManager`
- `DeformationModelAssembler`
- `DeformationModelEnergy`
- `runSimCore` 里的 shared runtime setup

这些都建立在当前阶段已经完成的资产层之上，但不属于资产层本身。

下一阶段进入配置入口与 runtime 前置整理： [Phase 02](phase02_config_runtime_entry.md)
