# Phase 04：`CubicMeshDeformationModel` 单元模型

## 在系统中的位置

`CubicMeshDeformationModel` 是**单个六面体单元**的局部有限元模型，组合几何（8 顶点参考/当前位置）、弹性本构（`ElasticModel3DDeformationGradient`）和塑性本构（`PlasticModel3DDeformationGradient`），提供局部能量、梯度、Hessian 和参数导数：

```text
SimulationMesh (8 顶点 + 材料)
  → ElasticModel (能密度 ψ, 应力 P, 切线 dP/dF)
  → PlasticModel (Fp, Fp⁻¹, detFp 及对参数 a 的导数)
  → CubicMeshDeformationModel (积分 + 求导)
  → DeformationModelAssembler (全局装配)
```

## Element Contract

| 属性 | 值 |
| --- | --- |
| 顶点数 | 8 |
| 位移自由度 | 24 |
| 积分点数 | 8（$2 \times 2 \times 2$ Gauss） |
| 单元类型 | 三线性 8 节点 cubic/hexa |
| 材料接口 | deformation gradient (`ElasticModel3DDeformationGradient`) |
| Plastic 接口 | deformation gradient 的 plastic split (`PlasticModel3DDeformationGradient`) |

## 1. 参考单元与三线性形函数

### 参考参数域

使用单位立方体 $(\alpha,\beta,\gamma) \in [0,1]^3$（非教材常见的 $[-1,1]^3$），与 `CubicMesh::computeBarycentricWeights(...)` 一致。8 个角点按 `CubicMesh` 的局部顶点顺序：

$$
(0,0,0),\ (1,0,0),\ (1,1,0),\ (0,1,0),\ (0,0,1),\ (1,0,1),\ (1,1,1),\ (0,1,1)
$$

### 形函数

8 个 shape function 统一写成：

$$
N_{abc}(\alpha,\beta,\gamma) = \alpha^a(1-\alpha)^{1-a} \cdot \beta^b(1-\beta)^{1-b} \cdot \gamma^c(1-\gamma)^{1-c}, \quad a,b,c \in \{0,1\}
$$

例如：$N_{000} = (1-\alpha)(1-\beta)(1-\gamma)$，$N_{110} = \alpha\beta(1-\gamma)$，$N_{111} = \alpha\beta\gamma$。

### 形函数梯度 `dN_dabc`

实现不存 $N_i$ 本身，而是直接存梯度。定义 $a_0 = 1-\alpha$，$b_0 = 1-\beta$，$g_0 = 1-\gamma$，把 8 个 shape function 组成列向量 $N(\boldsymbol{\xi}) \in \mathbb{R}^8$，其对局部坐标的 Jacobian 为：

$$
J_N(\boldsymbol{\xi}) = \frac{\partial N}{\partial \boldsymbol{\xi}} \in \mathbb{R}^{8 \times 3}
$$

代码中存的是它的转置：

$$
\texttt{dN\_dabc} = J_N^T \in \mathbb{R}^{3 \times 8}
$$

其中：

- `dN_dabc.row(0)` = $(\partial N / \partial \alpha)^T$
- `dN_dabc.row(1)` = $(\partial N / \partial \beta)^T$
- `dN_dabc.row(2)` = $(\partial N / \partial \gamma)^T$

对应代码 `fillShapeGradients(...)`：

```cpp
dN_dabc.row(0) << -b0*g0, b0*g0, beta*g0, -beta*g0, -b0*gamma, b0*gamma, beta*gamma, -beta*gamma;
dN_dabc.row(1) << -a0*g0, -alpha*g0, alpha*g0, a0*g0, -a0*gamma, -alpha*gamma, alpha*gamma, a0*gamma;
dN_dabc.row(2) << -a0*b0, -alpha*b0, -alpha*beta, -a0*beta, a0*b0, alpha*b0, alpha*beta, a0*beta;
```

存 $J_N^T$ 而非 $J_N$ 的原因：后续几何映射 $D_m = X \cdot J_N$ 写成代码时正好是 `X * dN_dabc.transpose()`。

## 2. 参考几何与形变梯度

### 参考 Jacobian $D_m$

把 8 个顶点参考坐标按列排成 $X = [X_0, \ldots, X_7] \in \mathbb{R}^{3 \times 8}$，参考映射 Jacobian：

$$
D_m = \frac{\partial X}{\partial \boldsymbol{\xi}} = X \cdot J_N(\boldsymbol{\xi}) = X \cdot \texttt{dN\_dabc}^T
$$

对应代码：

```cpp
void computeDm(const M3x8d& X, const M3x8d& dN_dabc, ES::M3d& Dm) const {
    Dm.noalias() = X * dN_dabc.transpose();
}
```

### 总形变梯度 $F_{\text{ref}}$

当前构型下 $x = [x_0, \ldots, x_7] \in \mathbb{R}^{3 \times 8}$，通过链式法则：

$$
F_{\text{ref}} = \frac{\partial x}{\partial X} = \frac{\partial x}{\partial \boldsymbol{\xi}} \cdot D_m^{-1} = x \cdot \texttt{dN\_dabc}^T \cdot D_m^{-1}
$$

推导：$x$ 和 $X$ 都通过 $\boldsymbol{\xi}$ 参数化，$\partial x / \partial \boldsymbol{\xi} = F_{\text{ref}} \cdot D_m$，因此 $F_{\text{ref}} = (\partial x / \partial \boldsymbol{\xi}) D_m^{-1}$。

对应代码：

```cpp
void computeF(const M3x8d& x, const M3x8d& dN_dabc, const ES::M3d& DmInv, ES::M3d& F) const {
    F.noalias() = x * dN_dabc.transpose() * DmInv;
}
```

## 3. $2 \times 2 \times 2$ 高斯积分

### 积分点坐标

标准 $[-1,1]$ 上的二点 Gauss-Legendre 节点是 $\pm 1/\sqrt{3}$。映射到 $[0,1]$：$\xi = (\hat{\xi}+1)/2$，得到：

$$
\xi_{1,2} = 0.5 \mp \frac{0.5}{\sqrt{3}}
$$

三维张量积 $2 \times 2 \times 2 = 8$ 个积分点。

### 积分权重

$[-1,1]$ 上每点权重 $\hat{w} = 1$，映射到 $[0,1]$ 后乘 Jacobian $1/2$，三维张量积：

$$
w_q = \frac{1}{2} \cdot \frac{1}{2} \cdot \frac{1}{2} = \frac{1}{8}
$$

对应代码：

```cpp
constexpr double kQuadratureWeight = 0.125;
```

二点 Gauss-Legendre 可精确积分到三次多项式，对三线性六面体足够。

## 4. 每个积分点的预计算：`QuadratureData`

构造函数中，对 8 个积分点各预计算一份：

```cpp
struct QuadratureData {
    M3x8d   dN_dabc;     // 该点的形函数梯度
    ES::M3d DmInv;       // 参考 Jacobian 的逆
    M3x8d   restBm;      // 应力→节点力的装配模板
    M9x24d  rest_dFdx;   // vec(Fref) 对 24 DOF 的导数模板
    double  weightDetJ;  // 积分权重 × |det Dm|
};
```

构造代码：

```cpp
fillShapeGradients(quadratureCoord[ia], quadratureCoord[ib], quadratureCoord[ig], quad.dN_dabc);
computeDm(restX, quad.dN_dabc, Dm);

quad.DmInv      = Dm.fullPivLu().inverse();
quad.weightDetJ = kQuadratureWeight * std::abs(Dm.determinant());
quad.restBm     = quad.weightDetJ * quad.DmInv.transpose() * quad.dN_dabc;
compute_dF_dx(quad.dN_dabc, quad.DmInv, quad.rest_dFdx);
```

### `restBm` 的含义

$$
\texttt{restBm} = w_q |\det D_m| \cdot D_m^{-T} \cdot \texttt{dN\_dabc}
$$

其中 $D_m^{-T} \texttt{dN\_dabc}$ 就是形函数在参考物理坐标下的梯度（链式法则 $\nabla_X N_i = D_m^{-T} \nabla_{\boldsymbol{\xi}} N_i$），再乘积分体积权重。它是纯参考几何量，运行时乘上塑性因子得到当前 `Bm`。

### `rest_dFdx` 的含义

$$
\texttt{rest\_dFdx} = \frac{\partial \operatorname{vec}(F_{\text{ref}})}{\partial x} \in \mathbb{R}^{9 \times 24}
$$

只依赖 `dN_dabc` 和 $D_m^{-1}$。运行时乘上 $F_p^{-1}$ 得到弹性形变梯度的导数。

## 5. `dF/dx` 的推导

定义 $G = \texttt{dN\_dabc}^T \cdot A \in \mathbb{R}^{8 \times 3}$（$A$ 是 $D_m^{-1}$ 或 $D_m^{-1} F_p^{-1}$），则 $F = x \cdot G$。

$F$ 对 $x$ 是线性的。扰动单个自由度 $x_{vi,dim}$：

$$
\frac{\partial F}{\partial x_{vi,dim}} = E_{dim,vi} \cdot G
$$

其中 $E_{dim,vi} \in \mathbb{R}^{3 \times 8}$ 是只在 $(dim, vi)$ 位置为 1 的基矩阵。效果：结果矩阵只有第 $dim$ 行非零，等于 $G$ 的第 $vi$ 行。

对应代码：

```cpp
const Eigen::Matrix<double, 8, 3> G = dN_dabc.transpose() * A;
for (int vi = 0; vi < 8; vi++) {
    for (int dim = 0; dim < 3; dim++) {
        ES::M3d dF = ES::M3d::Zero();
        dF.row(dim) = G.row(vi);
        dFdx.col(vi * 3 + dim) = Eigen::Map<const ES::V9d>(dF.data());
    }
}
```

把每个 $3 \times 3$ 的 $\partial F / \partial x_{vi,dim}$ 向量化为 $9 \times 1$，存入 `dFdx` 的一列，构成完整的 $9 \times 24$ 矩阵。

## 6. 弹塑性分解

标准乘法分解：

$$
F = F_e F_p, \quad F_e = F_{\text{ref}} F_p^{-1}
$$

`prepareData(...)` 中从塑性参数 $a$ 恢复：

```cpp
plasticModel->computeA(param, Fp.data());        // Fp
plasticModel->computeAInv(param, FpInv.data());   // Fp⁻¹
cacheData->detFp = plasticModel->compute_detA(param);  // det(Fp)
```

## 7. `prepareData(...)`：运行时缓存

输入当前 8 节点位置 $x[24]$、塑性参数 $a$、材料参数 $b$。

### 缓存字段

```cpp
ES::M3d Fp, FpInv;    double detFp;
ES::M3d Fref[8];      // 每个积分点的总形变梯度
ES::M3d Fe[8];        // 弹性形变梯度
ES::M3d U[8], V[8];   ES::V3d S[8];  // Fe 的 SVD
M9x24d  dFdx[8];      // Fe 对 24 DOF 的导数
M3x8d   Bm[8];        // 应力→节点力装配矩阵
```

### 逐积分点计算

```cpp
computeF(xMat, quad.dN_dabc, quad.DmInv, cacheData->Fref[qi]);
cacheData->Fe[qi] = cacheData->Fref[qi] * cacheData->FpInv;
computeSVD(cacheData->Fe[qi], cacheData->U[qi], cacheData->V[qi], cacheData->S[qi]);

compute_dF_dx(quad.dN_dabc, quad.DmInv * cacheData->FpInv, cacheData->dFdx[qi]);
cacheData->Bm[qi] = cacheData->detFp * cacheData->FpInv.transpose() * quad.restBm;
```

关键点：

- `dFdx[qi]`：取 $A = D_m^{-1} F_p^{-1}$，得到的是 $\partial \operatorname{vec}(F_e) / \partial x$
- `Bm[qi]`：$\det(F_p) \cdot F_p^{-T} \cdot \texttt{restBm}$，把参考几何模板更新到当前弹塑性构型

### 为什么做 SVD

```cpp
Eigen::JacobiSVD<ES::M3d, Eigen::NoQRPreconditioner> svd(Fe, ComputeFullU | ComputeFullV);
```

Stable Neo-Hookean 等本构模型在奇异值空间中计算更稳健。SVD 预缓存避免能量/梯度/Hessian 重复分解。代码还做符号修正（`det(U) < 0` 或 `det(V) < 0` 时翻转），保证朝向一致性。

## 8. 能量

$$
E_e = \sum_{q=1}^{8} \psi(F_{e,q};\, b) \cdot w_q |\det D_{m,q}| \cdot \det F_p
$$

对应代码：

```cpp
for (int qi = 0; qi < 8; qi++) {
    energy += elasticModel->compute_psi(materialParam, Fe[qi], U[qi], V[qi], S[qi])
              * quad[qi].weightDetJ * cacheData->detFp;
}
```

- `compute_psi`：弹性能密度 $\psi(F_e)$，由 `ElasticModel` 提供
- `weightDetJ`：$w_q |\det D_m|$
- `detFp`：塑性映射的体积因子

## 9. 对几何自由度 $x$ 的一阶导（节点力）

$$
\frac{\partial E_e}{\partial x} = \sum_q P_q \cdot B_{m,q}
$$

其中 $P_q = \partial \psi / \partial F_e$ 是第一类 Piola 应力。

对应代码：

```cpp
ES::M3d P;
elasticModel->compute_P(materialParam, Fe[qi], U[qi], V[qi], S[qi], P.data());
gradMap.noalias() += P * cacheData->Bm[qi];
```

这就是标准非线性有限元中"应力积分得到内力"的形式。`Bm` 负责把应力从积分点映射回 24 个节点自由度。

## 10. 对几何自由度 $x$ 的二阶导（切线刚度矩阵）

$$
\frac{\partial^2 E_e}{\partial x^2} = \sum_q \left(\frac{\partial \operatorname{vec}(F_{e,q})}{\partial x}\right)^T \frac{\partial P_q}{\partial F_{e,q}} \left(\frac{\partial \operatorname{vec}(F_{e,q})}{\partial x}\right)
$$

对应代码：

```cpp
ES::M9d dPdF;
elasticModel->compute_dPdF(materialParam, Fe[qi], U[qi], V[qi], S[qi], dPdF.data());
dPdF *= quad[qi].weightDetJ * cacheData->detFp;
hessMap.noalias() += cacheData->dFdx[qi].transpose() * dPdF * cacheData->dFdx[qi];
```

- `dPdF`：$9 \times 9$ 材料切线矩阵 $\partial \operatorname{vec}(P) / \partial \operatorname{vec}(F_e)$
- `dFdx`：$9 \times 24$ 几何 Jacobian
- 结果：$24 \times 24$ 局部切线刚度矩阵

## 11. 对塑性参数 $a$ 的导数

这是 `CubicMeshDeformationModel` 与纯弹性六面体最大的区别。

### 链式法则

$$
\frac{\partial F_e}{\partial a_i} = F_{\text{ref}} \cdot \frac{\partial F_p^{-1}}{\partial a_i}
$$

### 一阶导 `compute_dE_da`

$$
\frac{\partial E_e}{\partial a_i} = \sum_q \left( \frac{\partial V_q}{\partial a_i} \psi_q + V_q \cdot P_q : \frac{\partial F_e}{\partial a_i} \right)
$$

其中 $V_q = w_q |\det D_m| \det F_p$，$\partial V_q / \partial a_i$ 来自 $\partial \det F_p / \partial a_i$。

### 已实现的全部参数导数

| 方法 | 含义 |
| --- | --- |
| `compute_dE_da` | $\partial E / \partial a$ |
| `compute_d2E_da2` | $\partial^2 E / \partial a^2$ |
| `compute_d2E_dxda` | $\partial^2 E / \partial x \partial a$ |
| `compute_dE_db` | $\partial E / \partial b$（材料参数） |
| `compute_d2E_db2` | $\partial^2 E / \partial b^2$ |
| `compute_d2E_dxdb` | $\partial^2 E / \partial x \partial b$ |
| `compute_d2E_dadb` | $\partial^2 E / \partial a \partial b$ |

Contract 深度与 `DeformationModel` 抽象完全对齐。

## 验证

`cubicMeshDeformationModel_test.cpp` 覆盖：

- `AffineDeformationMatchesMaterialEnergyDensityTimesVolume` — 仿射形变下 energy = 能密度 × 体积（验证积分正确性）
- `FiniteDifferenceMatchesGradientAndHessian` — 有限差分对照一阶/二阶导数（验证解析求导正确性）

---

上一阶段：[Phase 03](phase03_simulation_mesh.md)
下一阶段：[Phase 05](phase05_energy_assembly.md)
