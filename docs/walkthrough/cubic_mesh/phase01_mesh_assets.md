# Phase 01：网格表示与资产生成

## `CubicMesh` 本体

`CubicMesh`（继承自 `VolumetricMesh`）表示一类受限的 hexa mesh：

- 所有 cube 边长相等、轴对齐
- 拓扑来自规则网格的 voxel 子集
- 不是任意 8 点 hexa 容器

### 局部顶点顺序

8 个局部顶点的固定顺序（跨层 contract，后续所有层都依赖）：

```text
0: (0,0,0)   1: (1,0,0)   2: (1,1,0)   3: (0,1,0)
4: (0,0,1)   5: (1,0,1)   6: (1,1,1)   7: (0,1,1)
```

对应构造辅助数组：

```cpp
int vtxI[8] = {0, 1, 1, 0, 0, 1, 1, 0};
int vtxJ[8] = {0, 0, 1, 1, 0, 0, 1, 1};
int vtxK[8] = {0, 0, 0, 0, 1, 1, 1, 1};
```

`computeBarycentricWeights`、`createFromCubicMesh`、`CubicMeshDeformationModel` 和表面提取都依赖这个顺序。

### `createFromUniformGrid(...)`

输入激活 voxel 的 `(i,j,k)` 索引集，输出共享顶点的 `CubicMesh`：

1. 遍历 voxel，8 个角点加入共享顶点集
2. 建立 `(i,j,k) → vertex id` 映射
3. 按固定顶点顺序组装 8 点单元连通性
4. 坐标默认放在规范域 `[-0.5, 0.5]^3`，世界尺度由工具层后处理

## `cubicMesher` 工具

### `uniform` 子命令

生成规则体块 cubic mesh：

```bash
cubicMesher uniform \
  --resolution 4 \
  --output-mesh cubic-box.veg \
  --output-surface cubic-box.obj \
  --size 2.0 --offset 1.5 -2.0 0.25 \
  --E 1e6 --nu 0.45 --density 1000
```

支持 `--resolution`、`--size`、`--offset`、`--E`/`--nu`/`--density`、可选 `--output-surface`。

### `triangle-mesh` 子命令

将输入三角网格体素化为 cubic mesh：

```bash
cubicMesher triangle-mesh \
  --input-mesh dragon.obj \
  --resolution 64 \
  --output-mesh dragon-cubic.veg \
  --output-surface dragon-cubic-surface.obj
```

执行链：

```text
CLI 参数解析 → TriangleMeshVoxelizerOptions
  → createTriangleMeshCubicMesh(...)
    → loadAndSanitizeTriangleMesh(...)     # 读取 + 清洗
    → validateTriangleMeshForVoxelization(...)  # 几何合法性
    → buildVoxelDomain(...)                # 构造 voxel domain
    → buildOccupiedVoxelIndices(...)       # 占据判定
    → CubicMesh::createFromUniformGrid(...)  # 规范域网格
    → applyDomainTransform(...)            # 还原世界坐标
    → applyUserTransform(...)              # 用户 scale/offset
```

CLI 层只做参数收集和基本检查（`resolution^3` 溢出、`scale > 0`），核心逻辑全在 `triangleMeshVoxelizer.cpp`。

#### 输入 mesh 清洗与校验

`loadAndSanitizeTriangleMesh(...)` 先做轻量清理：`removeInvalidTriangles` + `removeIsolatedVertices`。

`validateTriangleMeshForVoxelization(...)` 做三项严格检查：

| 检查 | 原因 |
| --- | --- |
| `isManifold(mesh)` | 流形前提 |
| `getExteriorEdges(...).empty()` | 必须封闭（无边界边） |
| `!isSelfIntersected(mesh)` | 不能自交 |

三项都是 `CGAL::Side_of_triangle_mesh` 的隐含前提 — 它需要一个闭合、非自交、方向一致的实体边界才能给出稳定的 inside/outside 判定。

#### Voxel Domain 构造

`buildVoxelDomain(...)` 的公式：

$$\text{interiorRes} = N - 2K$$

$$\text{voxelSide} = \frac{\text{baseSide}}{\text{interiorRes}}$$

$$\text{domainSide} = N \times \text{voxelSide}$$

- $N$ = `resolution`（整个 domain 每轴 voxel 数）
- $K$ = `padding-voxels`
- `baseSide` = 输入 mesh 包围盒 regularize 后的边长

`padding` 不是"bbox 外扩一点"，而是从固定的 $N$ 中预留 $2K$ 层空 voxel 给模型周围。$K$ 越大 → `interiorRes` 越小 → `voxelSide` 越大 → 网格越粗。

`VoxelDomain` 提供 `voxelCenter(i,j,k)` 和 `voxelBoundingBox(i,j,k)` 供后续占据判定使用。

#### 占据判定

`buildOccupiedVoxelIndices(...)` 三层循环遍历 $N^3$ 个 voxel，对每个 `(i,j,k)`：

```text
center_inside = insideQuery(voxelCenter(i,j,k))

if classify_mode == center:
    occupied = center_inside
elif classify_mode == conservative:
    occupied = center_inside || voxelIntersectsTriangleSurface(voxelBox)
```

- `insideQuery`：`CGAL::Side_of_triangle_mesh` 点内外查询
- `voxelIntersectsTriangleSurface`：两级过滤 — 先用 `BoundingBoxBVTree` AABB 粗筛候选三角形，再做精确 `intersectTriAABB`

`center` 模式快但会漏薄特征；`conservative` 模式对薄壳、尖角、擦边表面更稳，但多一次 voxel-triangle 相交查询。

#### 坐标还原

`createFromUniformGrid(...)` 输出的坐标固定在规范域 $[-0.5, 0.5]^3$。还原分两步，都通过 `applyLinearTransformation(pos, R)` 即 $X_{\text{new}} = \text{pos} + R \cdot X$ 实现：

1. **`applyDomainTransform`**：$R = \text{domainSide} \cdot I$，$\text{pos} = \text{worldBox.center}$  
   → 从规范域映射回 voxelization domain 的世界坐标

2. **`applyUserTransform`**：$R = \text{scale} \cdot I$，$\text{pos} = \text{offset}$  
   → 应用用户指定的额外缩放和平移

#### 具体例子

`cubicMesher_test.cpp` 中的 unit cube 测试：`resolution=4, padding=1, mode=center`

$$\text{interiorRes} = 4 - 2 = 2, \quad \text{voxelSide} = \frac{1}{2} = 0.5, \quad \text{domainSide} = 4 \times 0.5 = 2$$

模型占据内部 $2 \times 2 \times 2 = 8$ 个 voxel，共享顶点 $3 \times 3 \times 3 = 27$ 个。还原后包围盒 $1 \times 1 \times 1$，中心 $(0.5, 0.5, 0.5)$ — 与测试断言完全吻合。

#### 实现边界

- 遍历是 per-voxel $O(N^3)$，不是 per-triangle rasterization
- 只接受 `.obj` 格式
- 依赖 `PGO_HAS_CGAL`（geometry stack），没有时直接失败

## `.veg` + `.obj` 分工

- `.veg`：内部力学离散、体材料、质量、自由度
- `.obj`：显示、接触边界、surface embedding

`cubicMesher` 可选导出 surface `.obj`，`runSimCore` 总是显式读取 surface mesh。

## Committed 资产

| 目录 | 资产来源 |
| --- | --- |
| `examples/cubic-box/` | `uniform` |
| `examples/cubic-box-ipc/` | `uniform` |
| `examples/pulled-cubic-box-self-ipc/` | `uniform` |
| `examples/dragon-cubic/` | `triangle-mesh --resolution 64` |

## 验证

`cubicMesher_test.cpp` 覆盖：uniform vertex/element 数量、surface 导出、size/offset 影响 bounding box、triangle-mesh voxelization 基线。

---

下一阶段：[Phase 02](phase02_config_runtime_entry.md)
