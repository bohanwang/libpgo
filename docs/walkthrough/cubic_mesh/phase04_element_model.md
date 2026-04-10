# Phase 04：`CubicMeshDeformationModel` 单元模型

## 1. 这阶段要回答什么

这一阶段专门回答下面四个问题：

1. cubic/hexa 单元在 solver 里到底使用哪个局部 deformation model？
2. 这个模型怎样把 hexa 理论真正翻译成 repo 当前的代码？
3. 为什么实现里使用 `[0, 1]^3` 参数域，而教材里常写 `[-1, 1]^3`？
4. 当前 repo 到底只实现了 energy，还是 gradient、Hessian、参数导数也已经在位？

这一步如果不展开，后面很容易出现两种误判：

- 以为 cubic 只是把 tet 的单元模型硬扩成 8 个点
- 以为 repo 里只有 element energy 能算，global assembly 还没有真正可用

## 2. 这阶段在整条主线里解决什么

Phase 03 解决的是 `CubicMesh -> SimulationMesh` 的 solver-side handoff。  
Phase 04 要解决的是：

> 在 solver 看见一个 `SimulationMeshType::CUBIC` element 之后，究竟由哪个局部 FEM 对象负责计算它的能量、力和切线刚度。

当前 repo 给出的答案已经很明确：

- cubic/hexa 单元的局部 deformation model 是 `CubicMeshDeformationModel`

## 3. 当前 repo 已经完成了什么

当前实现已经完成的不只是“新增了一个类文件”，而是：

- `CubicMeshDeformationModel` 类已存在
- 它遵守 `DeformationModel` 抽象接口
- 它针对 8 节点 cubic 单元工作
- 它使用 `2 x 2 x 2 = 8` 个 Gauss 点
- 它已经实现：
  - `prepareData(...)`
  - `computeEnergy(...)`
  - `compute_dE_dx(...)`
  - `compute_d2E_dx2(...)`
  - 参数导数相关接口
- 它复用了现有的：
  - `ElasticModel3DDeformationGradient`
  - `PlasticModel3DDeformationGradient`

因此当前 repo 的真实状态不是“有一个待接线的 hexa 单元草稿”，而是：

> 已经有一个可供 manager / assembler / energy 直接消费的 cubic element model。

## 4. 这个类的固定 contract 是什么

当前 `CubicMeshDeformationModel` 的 element contract 可以概括成：

- 顶点数：`8`
- 位移自由度：`24`
- 积分点数：`8`
- 单元类型：三线性 8 节点 cubic / hexa
- 材料接口：deformation gradient 形式
- plastic 接口：deformation gradient 的 plastic split

这组 contract 非常重要，因为后面至少三层逻辑会依赖它：

- `DeformationModelManager` 创建对象时要喂 `restPositions[24]`
- `DeformationModelAssembler` 需要知道局部自由度规模是 24
- 测试会直接断言 `getNumVertices()`、`getNumDOFs()`、`getNumMaterialLocations()`

## 5. 构造函数在预计算什么

构造函数接收的是单元 rest positions `restPositions[24]`。之后它会为 8 个积分点预计算：

- `dN_dabc`
- `DmInv`
- `weightDetJ`
- `restBm`
- `rest_dFdx`

### 5.1 `dN_dabc`

这是 8 个形函数对 `(alpha, beta, gamma)` 的梯度组织结果。它对应三线性 hexa 单元的形函数梯度矩阵，只不过参数域已经换成了 `[0, 1]^3`。

### 5.2 `DmInv`

`Dm` 是材料空间 Jacobian，`DmInv` 是它的逆。  
它把 rest 构型几何编码成一个局部线性映射，是后续恢复 `F` 的关键。

### 5.3 `weightDetJ`

这是：

```text
Gauss 权重 * |det(Dm)|
```

也就是把参考域积分转到材料空间后，当前积分点对单元总能量的几何权重。

### 5.4 `restBm`

`restBm` 是后续局部力积分反复使用的缓存项。把它预计算出来的意义在于：

- 单元 rest 几何对所有后续状态都相同
- 没必要在每次 `prepareData(...)` 时重复构造

### 5.5 `rest_dFdx`

这表示形变梯度对节点自由度的 Jacobian 在 rest 几何下的基础形式。  
它是后面 Hessian 计算的关键缓存。

因此构造函数真正做的事情可以概括成：

> 把“与当前状态无关、只依赖单元 rest 构型”的几何量全部提前编译进对象。

## 6. `prepareData(...)` 到底在恢复什么

runtime 中，`prepareData(...)` 会读取：

- 当前 8 个节点位置 `x`
- plastic 参数 `param`
- material 参数 `materialParam`

然后对每个积分点依次计算：

1. `Fp`
2. `FpInv`
3. `detFp`
4. `Fref`
5. `Fe = Fref * FpInv`
6. `SVD(Fe)`
7. `dFdx`
8. `Bm`

这一步最重要的不是变量名，而是它把 hexa 几何主链和 repo 现有的 plastic split 统一了起来。

### 6.1 如果只看 hexa 几何部分

三线性 hexa 单元最核心的几何公式可以写成：

$$
F_q = D_s H_q (D_m H_q)^{-1}
$$

当前实现里，对应的是：

- `xMat` 代表当前顶点矩阵 `D_s`
- `quad.dN_dabc` 对应形函数梯度
- `quad.DmInv` 对应 `(D_m H_q)^{-1}`

### 6.2 为什么还要分 `Fref` 和 `Fe`

repo 当前的 volumetric 材料接口并不直接吃几何意义上的总形变梯度，而是延续了 tet 主线已有的 deformation-gradient + plastic split 抽象。

因此实现里先计算：

- `Fref`

然后再得到：

- `Fe = Fref * FpInv`

真正送给材料模型的是 `Fe`。  
所以更准确的表述应当是：

- hexa 几何主链首先落成 `Fref`
- 再通过 plastic split 落成材料层消费的 `Fe`

## 7. 形函数、顶点顺序和参数域怎样对齐

### 7.1 顶点顺序仍然严格沿用 `CubicMesh`

当前实现没有再发明另一套 hexa 局部编号，而是继续使用：

```text
000, 100, 110, 010, 001, 101, 111, 011
```

这使得：

- `CubicMesh` 的局部顶点顺序
- `SimulationMesh::createFromCubicMesh(...)` 的 element 顺序
- `CubicMeshDeformationModel` 的形函数顺序

三者完全一致。

### 7.2 为什么用 `[0, 1]^3`

标准 hexa 推导常使用参考域 `[-1, 1]^3`。当前实现改用 `[0, 1]^3`，不是因为理论变了，而是为了和现有 `CubicMesh::computeBarycentricWeights(...)` 一致。

二者通过仿射换元联系：

$$
\xi_1 = 2\alpha - 1,\quad
\xi_2 = 2\beta - 1,\quad
\xi_3 = 2\gamma - 1
$$

因此：

```text
标准参考域表达
  <=> repo 中的 alpha / beta / gamma 参数化
```

是完全等价的。

## 8. `computeEnergy(...)`、`compute_dE_dx(...)`、`compute_d2E_dx2(...)` 分别在做什么

### 8.1 `computeEnergy(...)`

这一层做的是 8 点积分求和：

$$
E_e = \sum_{q=1}^{8} \Psi(F_q^e)\, w_q \lvert \det(D_m H_q) \rvert \det(F_p)
$$

因此 element total energy 不是某个单点近似，而是完整地在 8 个积分点上累计。

### 8.2 `compute_dE_dx(...)`

这一层逐积分点计算一阶 Piola 应力 `P`，再与 `Bm` 相乘累加局部梯度：

$$
f^e = \sum_q P_q B_{m,q}
$$

这意味着当前实现并不是“先有 energy，gradient 以后再补”，而是：

- energy
- force / gradient

都已经在同一个 element contract 内闭合。

### 8.3 `compute_d2E_dx2(...)`

这一层逐积分点计算 `dPdF`，再通过 `dFdx` 累积局部切线刚度：

$$
K_e = \sum_q \left(\frac{\partial F_q}{\partial x}\right)^T
\frac{\partial P_q}{\partial F_q}
\left(\frac{\partial F_q}{\partial x}\right)
$$

这一步的意义非常大，因为后面 assembler 和 implicit runtime 都依赖它。

## 9. 当前 repo 到底实现到了什么深度

这里有必要明确一下，不然后面很容易低估现状。

当前 `CubicMeshDeformationModel` 已经不只实现了：

- `computeEnergy`

它实际上还已经实现了：

- `compute_dE_dx`
- `compute_d2E_dx2`
- `compute_dE_da`
- `compute_d2E_da2`
- `compute_d2E_dxda`
- `compute_dE_db`
- `compute_d2E_db2`
- `compute_d2E_dxdb`
- `compute_d2E_dadb`

也就是说，当前 repo 中 cubic element model 的 contract 深度已经和主线 deformation model 抽象对齐，而不是一个只够最小静态示例的半成品。

## 10. 当前 repo 的验证证据

`tests/core/energy/cubicMeshDeformationModel_test.cpp` 已经覆盖了两类最核心的 element-level 行为：

- `AffineDeformationMatchesMaterialEnergyDensityTimesVolume`
  - 检查仿射形变下的能量是否与材料能量密度乘体积一致
- `FiniteDifferenceMatchesGradientAndHessian`
  - 用有限差分同时对照一阶和二阶导数

这些测试的重要性在于：

- 它们验证的不是 manager 或 runtime glue
- 而是 cubic element model 自身的局部数学 contract

因此当前 repo 最稳妥的表述是：

> `CubicMeshDeformationModel` 已经具备 committed 的 element-level 数值一致性证据。

## 11. 这一阶段的边界

这一阶段仍然不展开：

- `DeformationModelManager` 怎样为 cubic 创建对象
- assembler 怎样把局部量装配回全局
- `runSimCore` 怎样构造 `W`、`M` 和 runtime integrator

它解决的是单元级 deformation model 本身，而不是之后的全局装配和 runtime。

上一阶段： [Phase 03](phase03_simulation_mesh.md)  
下一阶段： [Phase 05](phase05_energy_assembly.md)
