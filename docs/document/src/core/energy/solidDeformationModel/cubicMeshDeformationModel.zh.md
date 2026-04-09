# CubicMeshDeformationModel

源码文件：

- `src/core/energy/solidDeformationModel/cubicMeshDeformationModel.h`
- `src/core/energy/solidDeformationModel/cubicMeshDeformationModel.cpp`

## 1. 这个类在系统里的位置

`CubicMeshDeformationModel` 是 **单个六面体体单元** 的局部有限元模型。它把下面三类信息组合起来：

- 几何：
  - 8 个顶点的参考位置 `restPositions[24]`
  - 当前 8 个顶点的当前位置 `x[24]`
- 弹性本构：
  - `ElasticModel3DDeformationGradient`
- 塑性本构：
  - `PlasticModel3DDeformationGradient`

然后对这个单元提供：

- 局部能量 `E`
- 对几何自由度的梯度 `dE/dx`
- 对几何自由度的 Hessian `d²E/dx²`
- 对塑性参数 `a` 的导数
- 对材料参数 `b` 的导数
- von Mises 应力和最大应变等后处理量

从对象关系看，它位于下面这条链路中：

```text
SimulationMesh
  -> ElasticModel
  -> PlasticModel
  -> CubicMeshDeformationModel
  -> DeformationModelAssembler
```

其中：

- `SimulationMesh` 提供单元的 8 个顶点及其参考几何
- `ElasticModel` 提供弹性能密度 $\psi(F_e)$、一阶应力 $P$ 和材料切线 $dP/dF$
- `PlasticModel` 提供塑性映射 $F_p$ 及其参数化导数
- `CubicMeshDeformationModel` 把几何和本构组合成一个真正可积分、可求导的六面体单元
- `DeformationModelAssembler` 再把所有单元的局部量装配成全局系统

---

## 2. 对应的离散对象

头文件接口见：

```cpp
class CubicMeshDeformationModel : public DeformationModel {
public:
    CubicMeshDeformationModel(const double restPositions[24], ElasticModel* elasticModel, PlasticModel* plasticModel);

    virtual double computeEnergy(const CacheData* cacheData) const override;
    virtual void   compute_dE_dx(const CacheData* cacheData, double* grad) const override;
    virtual void   compute_d2E_dx2(const CacheData* cacheData, double* hess) const override;

    virtual void compute_dE_da(const CacheData* cacheData, double* grad) const override;
    virtual void compute_d2E_da2(const CacheData* cacheData, double* hess) const override;
    virtual void compute_d2E_dxda(const CacheData* cacheData, double* hess) const override;

    virtual void compute_dE_db(const CacheData* cacheData, double* grad) const override;
    virtual void compute_d2E_db2(const CacheData* cacheData, double* hess) const override;
    virtual void compute_d2E_dxdb(const CacheData* cacheData, double* hess) const override;

    virtual void compute_d2E_dadb(const CacheData* cacheData, double* hess) const override;

    virtual int getNumVertices() const override { return 8; }
    virtual int getNumDOFs() const override { return 24; }
    virtual int getNumMaterialLocations() const override { return 8; }
};
```

这个接口说明了三个基本事实：

1. 一个 cubic 单元固定有 8 个顶点。
2. 每个顶点有 3 个几何自由度，因此总自由度数是 `24`。
3. 单元内部有 `8` 个材料采样点，也就是 `2 x 2 x 2` 的 Gauss 积分点。

---

## 3. 理论模型：8 节点三线性六面体

### 3.1 参考单元

这个类使用的参考参数域是单位立方体：

$$
(\alpha,\beta,\gamma) \in [0,1]^3.
$$

8 个角点对应：

$$
(0,0,0), (1,0,0), (1,1,0), (0,1,0), (0,0,1), (1,0,1), (1,1,1), (0,1,1).
$$

这与 `CubicMesh` 的局部顶点顺序一致。

### 3.2 三线性形函数

8 个 shape functions 可以统一写成：

$$
N_{abc}(\alpha,\beta,\gamma)
=
\alpha^a (1-\alpha)^{1-a}
\beta^b (1-\beta)^{1-b}
\gamma^c (1-\gamma)^{1-c},
$$

其中 $a,b,c \in \{0,1\}$。

例如：

$$
\begin{aligned}
N_{000} &= (1-\alpha)(1-\beta)(1-\gamma), \\
N_{100} &= \alpha(1-\beta)(1-\gamma), \\
N_{110} &= \alpha\beta(1-\gamma), \\
N_{111} &= \alpha\beta\gamma.
\end{aligned}
$$

实现里并不显式保存 $N_i$，而是直接保存它们的梯度 $\frac{\partial N}{\partial \alpha}, \frac{\partial N}{\partial \beta}, \frac{\partial N}{\partial \gamma}$，对应函数：

```cpp
void CubicMeshDeformationModelInternal::fillShapeGradients(double alpha, double beta, double gamma,
                                                           M3x8d& dN_dabc) const
```

关键实现：

```cpp
dN_dabc.row(0) << -b0 * g0, b0 * g0, beta * g0, -beta * g0, -b0 * gamma, b0 * gamma, beta * gamma,
    -beta * gamma;
dN_dabc.row(1) << -a0 * g0, -alpha * g0, alpha * g0, a0 * g0, -a0 * gamma, -alpha * gamma, alpha * gamma,
    a0 * gamma;
dN_dabc.row(2) << -a0 * b0, -alpha * b0, -alpha * beta, -a0 * beta, a0 * b0, alpha * b0, alpha * beta,
    a0 * beta;
```

这里：

- `row(0)` 是对 `alpha` 的偏导
- `row(1)` 是对 `beta` 的偏导
- `row(2)` 是对 `gamma` 的偏导

### 3.3 用 Jacobian 的方式理解 `dN_dabc`

如果把 8 个 shape function 视为一个 8 维向量，而把局部参考坐标视为一个 3 维向量，那么这部分可以非常自然地写成一个 Jacobian。

先定义局部坐标向量：

$$
\boldsymbol{\xi}
=
\begin{bmatrix}
\alpha \\
\beta \\
\gamma
\end{bmatrix}
\in \mathbb{R}^3.
$$

再把 8 个 shape function 组织成一个列向量：

$$
N(\boldsymbol{\xi})
=
\begin{bmatrix}
N_0 \\
N_1 \\
N_2 \\
N_3 \\
N_4 \\
N_5 \\
N_6 \\
N_7
\end{bmatrix}
\in \mathbb{R}^8.
$$

按照本文件使用的局部顶点顺序

$$
000,\ 100,\ 110,\ 010,\ 001,\ 101,\ 111,\ 011,
$$

这个向量可以写成：

$$
N(\boldsymbol{\xi})
=
\begin{bmatrix}
(1-\alpha)(1-\beta)(1-\gamma) \\
\alpha(1-\beta)(1-\gamma) \\
\alpha\beta(1-\gamma) \\
(1-\alpha)\beta(1-\gamma) \\
(1-\alpha)(1-\beta)\gamma \\
\alpha(1-\beta)\gamma \\
\alpha\beta\gamma \\
(1-\alpha)\beta\gamma
\end{bmatrix}.
$$

于是它对局部坐标的 Jacobian 就是：

$$
J_N(\boldsymbol{\xi})
=
\frac{\partial N}{\partial \boldsymbol{\xi}}
=
\begin{bmatrix}
\frac{\partial N_0}{\partial \alpha} & \frac{\partial N_0}{\partial \beta} & \frac{\partial N_0}{\partial \gamma} \\
\frac{\partial N_1}{\partial \alpha} & \frac{\partial N_1}{\partial \beta} & \frac{\partial N_1}{\partial \gamma} \\
\vdots & \vdots & \vdots \\
\frac{\partial N_7}{\partial \alpha} & \frac{\partial N_7}{\partial \beta} & \frac{\partial N_7}{\partial \gamma}
\end{bmatrix}
\in \mathbb{R}^{8 \times 3}.
$$

也可以理解为：

$$
J_N(\boldsymbol{\xi})
=
\begin{bmatrix}
\frac{\partial N}{\partial \alpha} &
\frac{\partial N}{\partial \beta} &
\frac{\partial N}{\partial \gamma}
\end{bmatrix},
$$

其中每一列都是一个 8 维向量。

为了简化记号，定义：

$$
a_0 = 1-\alpha,\qquad b_0 = 1-\beta,\qquad g_0 = 1-\gamma.
$$

那么三列 Jacobian 分别是：

$$
\frac{\partial N}{\partial \alpha}
=
\begin{bmatrix}
-b_0 g_0 \\
\;\;b_0 g_0 \\
\;\;\beta g_0 \\
-\beta g_0 \\
-b_0 \gamma \\
\;\;b_0 \gamma \\
\;\;\beta \gamma \\
-\beta \gamma
\end{bmatrix},
$$

$$
\frac{\partial N}{\partial \beta}
=
\begin{bmatrix}
-a_0 g_0 \\
-\alpha g_0 \\
\;\;\alpha g_0 \\
\;\;a_0 g_0 \\
-a_0 \gamma \\
-\alpha \gamma \\
\;\;\alpha \gamma \\
\;\;a_0 \gamma
\end{bmatrix},
$$

$$
\frac{\partial N}{\partial \gamma}
=
\begin{bmatrix}
-a_0 b_0 \\
-\alpha b_0 \\
-\alpha \beta \\
-a_0 \beta \\
\;\;a_0 b_0 \\
\;\;\alpha b_0 \\
\;\;\alpha \beta \\
\;\;a_0 \beta
\end{bmatrix}.
$$

因此完整 Jacobian 可以写成：

$$
J_N(\boldsymbol{\xi})
=
\begin{bmatrix}
-b_0 g_0 & -a_0 g_0 & -a_0 b_0 \\
\;\;b_0 g_0 & -\alpha g_0 & -\alpha b_0 \\
\;\;\beta g_0 & \;\;\alpha g_0 & -\alpha \beta \\
-\beta g_0 & \;\;a_0 g_0 & -a_0 \beta \\
-b_0 \gamma & -a_0 \gamma & \;\;a_0 b_0 \\
\;\;b_0 \gamma & -\alpha \gamma & \;\;\alpha b_0 \\
\;\;\beta \gamma & \;\;\alpha \gamma & \;\;\alpha \beta \\
-\beta \gamma & \;\;a_0 \gamma & \;\;a_0 \beta
\end{bmatrix}.
$$

源码中的 `dN_dabc` 不是这个 Jacobian 本身，而是它的转置：

$$
dN\_dabc = J_N(\boldsymbol{\xi})^T \in \mathbb{R}^{3 \times 8}.
$$

所以：

- `dN_dabc.row(0)` 对应 $\left(\frac{\partial N}{\partial \alpha}\right)^T$
- `dN_dabc.row(1)` 对应 $\left(\frac{\partial N}{\partial \beta}\right)^T$
- `dN_dabc.row(2)` 对应 $\left(\frac{\partial N}{\partial \gamma}\right)^T$

这就和实现完全对应起来了：

```cpp
dN_dabc.row(0) << -b0 * g0, b0 * g0, beta * g0, -beta * g0, -b0 * gamma, b0 * gamma, beta * gamma,
    -beta * gamma;
dN_dabc.row(1) << -a0 * g0, -alpha * g0, alpha * g0, a0 * g0, -a0 * gamma, -alpha * gamma, alpha * gamma,
    a0 * gamma;
dN_dabc.row(2) << -a0 * b0, -alpha * b0, -alpha * beta, -a0 * beta, a0 * b0, alpha * b0, alpha * beta,
    a0 * beta;
```

为什么代码更适合存 $J_N^T$ 而不是 $J_N$？因为后面几何映射写成矩阵形式时更自然。

若把 8 个顶点坐标写成：

$$
X = [X_0, X_1, \dots, X_7] \in \mathbb{R}^{3 \times 8},
$$

则几何映射可以写成：

$$
x(\boldsymbol{\xi}) = X N(\boldsymbol{\xi}).
$$

对局部坐标求导：

$$
\frac{\partial x}{\partial \boldsymbol{\xi}}
=
X \frac{\partial N}{\partial \boldsymbol{\xi}}
=
X J_N(\boldsymbol{\xi}).
$$

而在代码中使用的是：

$$
X\, dN\_dabc^T,
$$

因为：

$$
dN\_dabc^T = J_N(\boldsymbol{\xi}).
$$

这也解释了为什么源码里会同时出现 `computeDm(...)` 和 `computeF(...)`：

```cpp
Dm.noalias() = X * dN_dabc.transpose();
F.noalias()  = x * dN_dabc.transpose() * DmInv;
```

其中 `dN_dabc.transpose()` 正是上面推导出来的 Jacobian
$$
\frac{\partial N}{\partial \boldsymbol{\xi}}.
$$

---

## 4. 参考几何、形变梯度与高斯积分

### 4.1 参考映射

设参考构型下 8 个顶点坐标为：

$$
X = [X_0, X_1, \dots, X_7] \in \mathbb{R}^{3 \times 8},
$$

当前构型下 8 个顶点坐标为：

$$
x = [x_0, x_1, \dots, x_7] \in \mathbb{R}^{3 \times 8}.
$$

参考单元到参考构型的 Jacobian 为：

$$
D_m = \frac{\partial X}{\partial(\alpha,\beta,\gamma)}
= X \left(\frac{\partial N}{\partial(\alpha,\beta,\gamma)}\right)^T.
$$

当前构型相对参考构型的总 deformation gradient 为：

$$
F_{\mathrm{ref}}
=
\frac{\partial x}{\partial X}
=
x \left(\frac{\partial N}{\partial(\alpha,\beta,\gamma)}\right)^T D_m^{-1}.
$$

下面把这两条公式按矩阵微分的方式完整展开。

#### 第一步：写出参考单元到物理单元的插值映射

参考构型中的位置场由 shape function 插值得到：

$$
X(\boldsymbol{\xi})
=
\sum_{i=0}^{7} X_i N_i(\boldsymbol{\xi}).
$$

如果把 8 个顶点坐标按列排成矩阵

$$
X =
\begin{bmatrix}
\vert & \vert &        & \vert \\
X_0   & X_1   & \cdots & X_7 \\
\vert & \vert &        & \vert
\end{bmatrix}
\in \mathbb{R}^{3\times 8},
$$

而把 shape function 组织成列向量

$$
N(\boldsymbol{\xi})
=
\begin{bmatrix}
N_0(\boldsymbol{\xi}) \\
N_1(\boldsymbol{\xi}) \\
\vdots \\
N_7(\boldsymbol{\xi})
\end{bmatrix}
\in \mathbb{R}^{8},
$$

那么插值映射就可以写成

$$
X(\boldsymbol{\xi}) = X\,N(\boldsymbol{\xi}).
$$

同理，当前构型中的位置场也满足

$$
x(\boldsymbol{\xi}) = x\,N(\boldsymbol{\xi}),
$$

其中

$$
x = [x_0, x_1, \dots, x_7] \in \mathbb{R}^{3\times 8}.
$$

#### 第二步：对局部坐标求导

由于 `X` 和 `x` 都只是顶点坐标矩阵，对局部坐标 $\boldsymbol{\xi}$ 来说是常量，所以：

$$
\frac{\partial X}{\partial \boldsymbol{\xi}}
=
X \frac{\partial N}{\partial \boldsymbol{\xi}},
$$

$$
\frac{\partial x}{\partial \boldsymbol{\xi}}
=
x \frac{\partial N}{\partial \boldsymbol{\xi}}.
$$

而前面已经定义过

$$
J_N(\boldsymbol{\xi}) = \frac{\partial N}{\partial \boldsymbol{\xi}} \in \mathbb{R}^{8\times 3}.
$$

因此：

$$
\frac{\partial X}{\partial \boldsymbol{\xi}}
=
X J_N(\boldsymbol{\xi}),
$$

$$
\frac{\partial x}{\partial \boldsymbol{\xi}}
=
x J_N(\boldsymbol{\xi}).
$$

这两个量都是 `3 x 3` 矩阵。

#### 第三步：得到参考 Jacobian $D_m$

把局部坐标写成

$$
\boldsymbol{\xi} =
\begin{bmatrix}
\alpha \\
\beta \\
\gamma
\end{bmatrix},
$$

那么参考映射 Jacobian 就是

$$
D_m
=
\frac{\partial X}{\partial \boldsymbol{\xi}}
=
X J_N(\boldsymbol{\xi}).
$$

由于代码里保存的是

$$
dN\_dabc = J_N(\boldsymbol{\xi})^T \in \mathbb{R}^{3\times 8},
$$

所以也可以写成

$$
D_m
=
X\, dN\_dabc^T.
$$

因此更统一的记号应该写成

$$
D_m = \frac{\partial X}{\partial(\alpha,\beta,\gamma)}
= X \frac{\partial N}{\partial(\alpha,\beta,\gamma)}.
$$



只有在和源码中的变量 `dN_dabc` 对应时，才需要写成

$$
dN\_dabc^T = \frac{\partial N}{\partial(\alpha,\beta,\gamma)}.
$$

所以：

- 数学记号层面：
  $$
  D_m = X \frac{\partial N}{\partial(\alpha,\beta,\gamma)}
  $$
- 代码变量层面：
  $$
  D_m = X\, dN\_dabc^T
  $$

这两种写法表达的是同一件事。

#### 第四步：由链式法则得到 `F_ref`

总 deformation gradient 是“当前坐标对参考坐标”的导数：

$$
F_{\mathrm{ref}} = \frac{\partial x}{\partial X}.
$$

但当前实现里，`x` 和 `X` 都不是直接用对方做自变量，而都是通过局部坐标 $\boldsymbol{\xi}$ 参数化的，所以要用链式法则：

$$
\frac{\partial x}{\partial \boldsymbol{\xi}}
=
\frac{\partial x}{\partial X}
\frac{\partial X}{\partial \boldsymbol{\xi}}.
$$

也就是：

$$
\frac{\partial x}{\partial \boldsymbol{\xi}}
=
F_{\mathrm{ref}} D_m.
$$

因此：

$$
F_{\mathrm{ref}}
=
\frac{\partial x}{\partial \boldsymbol{\xi}} D_m^{-1}.
$$

再把

$$
\frac{\partial x}{\partial \boldsymbol{\xi}}
=
x J_N(\boldsymbol{\xi})
=
x\, dN\_dabc^T
$$

代进去，就得到：

$$
F_{\mathrm{ref}}
=
x J_N(\boldsymbol{\xi}) D_m^{-1}
=
x\, dN\_dabc^T D_m^{-1}.
$$

这正是代码里 `computeF(...)` 时的计算式。

#### 第五步：和实现逐字对应

源码中对应的是两个更明确的 helper：

```cpp
void CubicMeshDeformationModelInternal::computeDm(
    const M3x8d& X, const M3x8d& dN_dabc, ES::M3d& Dm) const {
    Dm.noalias() = X * dN_dabc.transpose();
}
void CubicMeshDeformationModelInternal::computeF(
    const M3x8d& x, const M3x8d& dN_dabc, const ES::M3d& DmInv, ES::M3d& F) const {
    F.noalias() = x * dN_dabc.transpose() * DmInv;
}
```

也就是说：

- `computeDm(...)` 只负责计算
  $$
  D_m = X\, dN\_dabc^T = \frac{\partial X}{\partial \boldsymbol{\xi}}
  $$
- `computeF(...)` 只负责计算
  $$
  F_{\mathrm{ref}} = x\, dN\_dabc^T D_m^{-1} = \frac{\partial x}{\partial X}
  $$

所以这一块的数学结构可以概括成：

$$
\text{参考坐标参数化}
\;\Longrightarrow\;
\frac{\partial(\cdot)}{\partial \boldsymbol{\xi}}
\;\Longrightarrow\;
\text{链式法则乘 } D_m^{-1}
\;\Longrightarrow\;
F_{\mathrm{ref}}.
$$

这也和现在的命名更一致：`computeDm` 只算参考 Jacobian，`computeF` 才真正算 deformation gradient。

### 4.2 2x2x2 高斯积分

构造函数中使用 8 个积分点：

```cpp
const double quadratureCoord[2] = {
    0.5 - 0.5 / std::sqrt(3.0),
    0.5 + 0.5 / std::sqrt(3.0),
};
```

这是把标准区间 `[-1,1]` 上的两点 Gauss-Legendre 积分映射到 `[0,1]` 得到的结果。

#### 第一步：从一维两点 Gauss-Legendre 公式出发

在标准区间 `[-1,1]` 上，二点 Gauss-Legendre 积分的节点是：

$$
\hat{\xi}_1 = -\frac{1}{\sqrt{3}},\qquad
\hat{\xi}_2 = \frac{1}{\sqrt{3}},
$$

对应权重都是：

$$
\hat{w}_1 = \hat{w}_2 = 1.
$$

因此一维积分近似可以写成：

$$
\int_{-1}^{1} f(\hat{\xi})\, d\hat{\xi}
\approx
f\!\left(-\frac{1}{\sqrt{3}}\right)
+
f\!\left(\frac{1}{\sqrt{3}}\right).
$$

这个二点公式可以精确积分到三次多项式，因此对于 trilinear 六面体这类低阶单元非常常见。

#### 第二步：把积分区间从 `[-1,1]` 映射到 `[0,1]`

本文件的局部参考坐标不是 `[-1,1]`，而是单位立方体 `[0,1]^3`。  
所以先对一维区间做线性映射：

$$
\xi = \frac{\hat{\xi} + 1}{2}.
$$

把两个标准高斯点代进去，就得到：

$$
\xi_1
=
\frac{1 - \frac{1}{\sqrt{3}}}{2}
=
0.5 - \frac{0.5}{\sqrt{3}},
$$

$$
\xi_2
=
\frac{1 + \frac{1}{\sqrt{3}}}{2}
=
0.5 + \frac{0.5}{\sqrt{3}}.
$$

这正是代码中的：

```cpp
const double quadratureCoord[2] = {
    0.5 - 0.5 / std::sqrt(3.0),
    0.5 + 0.5 / std::sqrt(3.0),
};
```

也就是说，这两个数并不是经验参数，而是标准高斯节点经过区间变换后的结果。

#### 第三步：一维权重为什么会变成 `1/2`

因为做变量替换后，

$$
\xi = \frac{\hat{\xi} + 1}{2}
\quad\Longrightarrow\quad
d\xi = \frac{1}{2} d\hat{\xi}.
$$

所以原来 `[-1,1]` 上的积分：

$$
\int_{-1}^{1} f(\hat{\xi})\, d\hat{\xi}
$$

在 `[0,1]` 上会变成：

$$
\int_0^1 f(\xi)\, d\xi
=
\frac{1}{2}
\int_{-1}^{1} f\!\left(\frac{\hat{\xi}+1}{2}\right) d\hat{\xi}.
$$

因此映射到 `[0,1]` 后，每个一维高斯点的权重都变成：

$$
w^{(1D)} = \frac{1}{2}.
$$

#### 第四步：从一维 2 点扩展到三维 8 点

三维积分使用张量积构造：

- `alpha` 方向取 2 个点
- `beta` 方向取 2 个点
- `gamma` 方向取 2 个点

因此总共有：

$$
2 \times 2 \times 2 = 8
$$

个积分点。

每个三维积分点都可以写成：

$$
(\alpha_i,\beta_j,\gamma_k),
\qquad i,j,k\in\{1,2\}.
$$

源码里正是通过三重循环枚举：

```cpp
for (int ia = 0; ia < 2; ia++) {
    for (int ib = 0; ib < 2; ib++) {
        for (int ig = 0; ig < 2; ig++) {
            ...
        }
    }
}
```

所以这 8 个点实际上是以下三维节点的全排列组合：

$$
\left(0.5 \pm \frac{0.5}{\sqrt{3}},
      0.5 \pm \frac{0.5}{\sqrt{3}},
      0.5 \pm \frac{0.5}{\sqrt{3}}\right).
$$

#### 第五步：为什么每个三维高斯点的权重都是 `1/8`

三维张量积积分的权重，就是三个一维权重的乘积：

$$
w_{ijk}
=
w_i^{(1D)} w_j^{(1D)} w_k^{(1D)}.
$$

而这里每个一维权重都是：

$$
w_i^{(1D)} = \frac{1}{2},
$$

因此三维每个点的权重就是：

$$
w_{ijk}
=
\frac{1}{2}\cdot\frac{1}{2}\cdot\frac{1}{2}
=
\frac{1}{8}.
$$

这就是源码里的：

```cpp
constexpr double kQuadratureWeight = 0.125;
```

#### 第六步：这 8 个高斯点在这个类里是怎么用的

得到每个积分点坐标后，构造函数会在每个点上预计算：

- shape function 梯度 `dN_dabc`
- 参考 Jacobian 的逆 `DmInv`
- 参考体积权重 `weightDetJ`
- `restBm`
- `rest_dFdx`

也就是说，`2 x 2 x 2` 高斯积分不只是“取了 8 个点做数值积分”，而是整个单元局部能量、梯度和 Hessian 计算的离散基础。

每个一维积分权重都会乘上一个区间缩放因子 `1/2`，因此三维中每个积分点的总权重是：

$$
w_q = \frac{1}{2}\cdot\frac{1}{2}\cdot\frac{1}{2} = \frac{1}{8}.
$$

源码里就是：

```cpp
constexpr int kNumQuadraturePoints = 8;
constexpr double kQuadratureWeight = 0.125;
```

### 4.3 `QuadratureData`

每个积分点预存的数据是：

```cpp
struct QuadratureData {
    M3x8d  dN_dabc;
    ES::M3d DmInv;
    M3x8d  restBm;
    M9x24d rest_dFdx;
    double weightDetJ;
};
```

各字段含义：

- `dN_dabc`
  - 该积分点的三线性 shape function 梯度
- `DmInv`
  - 参考几何 Jacobian 的逆
- `restBm`
  - 以参考几何为基础的应力装配矩阵
- `rest_dFdx`
  - `Fref` 对 24 个几何自由度的导数模板
- `weightDetJ`
  - 高斯权重乘 `|det Dm|`

构造函数中的预计算代码：

```cpp
ind->fillShapeGradients(quadratureCoord[ia], quadratureCoord[ib], quadratureCoord[ig], quad.dN_dabc);

ES::M3d Dm;
ind->computeDm(restX, quad.dN_dabc, Dm);

quad.DmInv      = Dm.fullPivLu().inverse();
quad.weightDetJ = kQuadratureWeight * std::abs(Dm.determinant());
quad.restBm     = quad.weightDetJ * quad.DmInv.transpose() * quad.dN_dabc;
ind->compute_dF_dx(quad.dN_dabc, quad.DmInv, quad.rest_dFdx);
```

这段做的事情是：

1. 用参考几何计算每个积分点的 $Dm$
2. 预存 $Dm^{-1}$
3. 预存积分权重和参考体积缩放
4. 预存后续梯度和 Hessian 会反复用到的矩阵模板

这些都是只依赖 rest geometry 的量，因此非常适合在构造函数里一次性算好。

### 4.4 `restBm` 的数学含义

`restBm` 对应的代码是：

```cpp
quad.restBm = quad.weightDetJ * quad.DmInv.transpose() * quad.dN_dabc;
```

先看不带积分权重的部分：

$$
D_m^{-T} dN\_{dabc}.
$$

这里 `dN_dabc` 存的是

$$
\left(\frac{\partial N}{\partial \boldsymbol{\xi}}\right)^T,
$$

因此每一列对应的是某个 shape function 在局部参考坐标下的梯度。  
根据坐标变换下的链式法则：

$$
\nabla_X N_i = D_m^{-T}\nabla_{\boldsymbol{\xi}} N_i.
$$

把 8 个 shape function 一起写成矩阵形式，就是：

$$
\left(\frac{\partial N}{\partial X}\right)^T
=
D_m^{-T} dN\_{dabc}.
$$

因此：

$$
\texttt{restBm}
=
\texttt{weightDetJ}\, D_m^{-T} dN\_{dabc}
=
\texttt{weightDetJ}\left(\frac{\partial N}{\partial X}\right)^T.
$$

也就是说，`restBm` 是：

- 参考构型物理坐标下的 shape function 梯度
- 再乘上该积分点对应的参考体积权重

它最终会进入节点力装配。  
后面在 `prepareData(...)` 中，真正用于当前弹塑性构型的一阶导模板是：

```cpp
cacheData->Bm[qi] = cacheData->detFp * cacheData->FpInv.transpose() * quad.restBm;
```

所以可以把关系理解成：

$$
\texttt{Bm}
=
\det(F_p)\,F_p^{-T}\,\texttt{restBm}.
$$

`restBm` 是纯参考几何模板，`Bm` 则是把塑性映射也纳入后的当前模板。

### 4.5 `rest_dFdx` 的数学含义

`rest_dFdx` 对应的代码是：

```cpp
ind->compute_dF_dx(quad.dN_dabc, quad.DmInv, quad.rest_dFdx);
```

它表示参考几何下 deformation gradient 对 24 个几何自由度的导数模板：

$$
\texttt{rest\_dFdx}
=
\frac{\partial \operatorname{vec}(F_{\mathrm{ref}})}{\partial x}
\in \mathbb{R}^{9\times 24}.
$$

这里之所以叫 “rest”，是因为它只使用：

- `dN_dabc`
- `DmInv`

还没有乘上当前塑性映射 `FpInv`。  
因此它完全由参考构型几何决定。

---

## 5. 弹塑性分解：为什么会出现 `Fref`、`Fe` 和 `Fp`

这个类实现的是基于 deformation gradient 的弹塑性模型。它使用的基本分解是：

$$
F = F_e F_p.
$$

因此：

$$
F_e = F F_p^{-1}.
$$

在代码中：

- `Fref` 表示当前相对 rest 配置的总形变梯度
- `Fp` 表示塑性映射
- `FpInv` 表示 $F_p^{-1}$
- `Fe = Fref * FpInv`

这部分在 `prepareData(...)` 中完成：

```cpp
ind->plasticModel->computeA(param, cacheData->Fp.data());
ind->plasticModel->computeAInv(param, cacheData->FpInv.data());
cacheData->detFp = ind->plasticModel->compute_detA(param);
```

以及：

```cpp
ind->computeF(xMat, quad.dN_dabc, quad.DmInv, cacheData->Fref[qi]);
cacheData->Fe[qi] = cacheData->Fref[qi] * cacheData->FpInv;
```

这里的 `param` 就是单元的塑性参数向量 `a`。  
`PlasticModel3DDeformationGradient` 负责把参数 `a` 解码成：

- $F_p$
- $F_p^{-1}$
- $\det F_p$
- 以及它们对参数 `a` 的导数

---

## 6. `prepareData(...)`：一次求值前的缓存准备

`prepareData(...)` 是这个类最关键的步骤。所有后续 `computeEnergy`、`compute_dE_dx`、`compute_d2E_dx2` 等函数，都是基于这里缓存的结果工作。

### 6.1 输入含义

```cpp
void CubicMeshDeformationModel::prepareData(const double* x,
                                            const double* param,
                                            const double* materialParam,
                                            CacheData* cacheDataBase) const
```

这里：

- `x`
  - 当前单元 8 个顶点的位置，共 24 个数
- `param`
  - 当前单元的塑性参数 `a`
- `materialParam`
  - 当前单元的材料参数 `b`

### 6.2 缓存的主要内容

`CacheData` 中的核心字段包括：

```cpp
Vp      plasticParam;
ES::M3d Fp, FpInv;
double  detFp;

ES::V3d x[kNumVertices];
ES::M3d Fref[kNumQuadraturePoints];
ES::M3d Fe[kNumQuadraturePoints];
ES::M3d U[kNumQuadraturePoints], V[kNumQuadraturePoints];
ES::V3d S[kNumQuadraturePoints];
M9x24d  dFdx[kNumQuadraturePoints];
M3x8d   Bm[kNumQuadraturePoints];
```

### 6.3 主要步骤

#### 第一步：读取当前顶点位置

```cpp
for (int vi = 0; vi < kNumVertices; vi++) {
    cacheData->x[vi] = ES::V3d(x[vi * 3 + 0], x[vi * 3 + 1], x[vi * 3 + 2]);
}
```

#### 第二步：从塑性参数恢复 `Fp`

```cpp
ind->plasticModel->computeA(param, cacheData->Fp.data());
ind->plasticModel->computeAInv(param, cacheData->FpInv.data());
cacheData->detFp = ind->plasticModel->compute_detA(param);
```

#### 第三步：预计算塑性参数的一阶/二阶导数

```cpp
ind->plasticModel->compute_ddetA_da(...);
ind->plasticModel->compute_d2detA_da2(...);
ind->plasticModel->compute_dAInv_da(...);
ind->plasticModel->compute_d2AInv_da2(...);
```

这些量后面会用于 `dE/da`、`d²E/da²` 和 `d²E/dxda`。

#### 第四步：逐积分点计算 `Fref`、`Fe`、SVD、`dFdx` 和 `Bm`

```cpp
ind->computeF(xMat, quad.dN_dabc, quad.DmInv, cacheData->Fref[qi]);
cacheData->Fe[qi] = cacheData->Fref[qi] * cacheData->FpInv;
ind->computeSVD(cacheData->Fe[qi], cacheData->U[qi], cacheData->V[qi], cacheData->S[qi]);

ind->compute_dF_dx(quad.dN_dabc, quad.DmInv * cacheData->FpInv, cacheData->dFdx[qi]);
cacheData->Bm[qi] = cacheData->detFp * cacheData->FpInv.transpose() * quad.restBm;
```

这里的意义是：

- `Fe[qi]`
  - 第 `qi` 个积分点的弹性 deformation gradient
- `U,S,V`
  - `Fe` 的 SVD
  - 某些弹性模型在奇异值空间中计算更稳健
- `dFdx[qi]`
  - `Fe` 对 24 个 nodal DOFs 的导数
- `Bm[qi]`
  - 用于装配 `dE/dx`

### 6.4 为什么要做 SVD

```cpp
Eigen::JacobiSVD<ES::M3d, Eigen::NoQRPreconditioner> svd(Fe, Eigen::ComputeFullU | Eigen::ComputeFullV);
```

很多本构模型，例如 Stable Neo-Hookean 或基于不变量/奇异值的模型，都会直接用：

- `F`
- 或 `F` 的奇异值 `S`
- 以及左右奇异向量 `U,V`

因此这里把 SVD 预先缓存下来，避免能量、梯度、Hessian 重复分解。

代码里还做了符号修正：

```cpp
if (U.determinant() < 0.0) { ... }
if (V.determinant() < 0.0) { ... }
```

这样能保证 SVD 的朝向一致性，减少反射带来的数值问题。

---

## 7. 几何导数：`dF/dx`

`compute_dF_dx(...)` 负责构造：

$$
\frac{\partial \mathrm{vec}(F_e)}{\partial x}
\in \mathbb{R}^{9 \times 24}.
$$

实现：

```cpp
void CubicMeshDeformationModelInternal::compute_dF_dx(const M3x8d& dN_dabc, const ES::M3d& A, M9x24d& dFdx) const {
    dFdx.setZero();

    const Eigen::Matrix<double, 8, 3> G = dN_dabc.transpose() * A;
    for (int vi = 0; vi < kNumVertices; vi++) {
        for (int dim = 0; dim < 3; dim++) {
            ES::M3d dF = ES::M3d::Zero();
            dF.row(dim) = G.row(vi);
            dFdx.col(vi * 3 + dim) = Eigen::Map<const ES::V9d>(dF.data());
        }
    }
}
```

这个构造的直觉是：

- 每个局部自由度 `x_{vi,dim}` 改变时
- 只会影响 `F` 的一行
- 因为 `F = x * (...)`

因此对每个自由度都可以写出一个 `3x3` 的 `dF`，再向量化填进 `dFdx` 的一列。

### 7.1 用全微分从 $F = xG$ 推导 `dFdx`

为了更清楚地理解 `compute_dF_dx(...)`，把

$$
G = dN\_{dabc}^T A \in \mathbb{R}^{8\times 3}
$$

记成一个固定矩阵。  
那么 deformation gradient 可以写成：

$$
F = xG,
$$

其中：

- $x \in \mathbb{R}^{3\times 8}$
- $G \in \mathbb{R}^{8\times 3}$
- $F \in \mathbb{R}^{3\times 3}$

在 `compute_dF_dx(...)` 里，`A` 被视为常量，因此 `G` 也是常量。  
于是全微分直接给出：

$$
dF = d(xG) = (dx)G + x(dG) = (dx)G.
$$

这说明 `F` 对 `x` 是线性的。

如果只扰动一个自由度 $x_{vi,dim}$，记对应的基矩阵为 $E_{dim,vi}\in\mathbb{R}^{3\times 8}$，则

$$
dx = E_{dim,vi}\, d x_{vi,dim}.
$$

代入全微分：

$$
dF = E_{dim,vi} G \, d x_{vi,dim}.
$$

因此：

$$
\frac{\partial F}{\partial x_{vi,dim}} = E_{dim,vi} G.
$$

而 $E_{dim,vi}G$ 的效果非常简单：

- 只保留 `G` 的第 `vi` 行
- 把它放到结果矩阵的第 `dim` 行
- 其余位置全为 0

这正是代码里的构造方式：

```cpp
ES::M3d dF = ES::M3d::Zero();
dF.row(dim) = G.row(vi);
```

最后再把这个 `3x3` 的 `dF` 按列或按内存布局展平成 `9x1` 向量，存入 `dFdx` 的一列：

```cpp
dFdx.col(vi * 3 + dim) = Eigen::Map<const ES::V9d>(dF.data());
```

于是整张矩阵

$$
dFdx = \frac{\partial \operatorname{vec}(F)}{\partial x}
\in \mathbb{R}^{9\times 24}
$$

就被显式构造出来了。

对 `rest_dFdx` 来说，取的是

$$
A = D_m^{-1},
$$

因此

$$
F_{\mathrm{ref}} = x\, dN\_{dabc}^T D_m^{-1},
$$

而

$$
\texttt{rest\_dFdx}
=
\frac{\partial \operatorname{vec}(F_{\mathrm{ref}})}{\partial x}.
$$

后面在 `prepareData(...)` 里真正当前态使用的是：

```cpp
ind->compute_dF_dx(quad.dN_dabc, quad.DmInv * cacheData->FpInv, cacheData->dFdx[qi]);
```

也就是把塑性映射 `FpInv` 进一步乘进去，得到当前弹性 deformation gradient 的导数模板。

---

## 8. 能量公式与代码实现

### 8.1 单元能量

代码：

```cpp
double energy = 0.0;
for (int qi = 0; qi < kNumQuadraturePoints; qi++) {
    energy += ind->elasticModel->compute_psi(cacheData->materialParam.data(), cacheData->Fe[qi].data(),
                                             cacheData->U[qi].data(), cacheData->V[qi].data(),
                                             cacheData->S[qi].data()) *
              ind->quad[qi].weightDetJ * cacheData->detFp;
}
```

可写成：

$$
E_e
=
\sum_{q=1}^{8}
\psi(F_{e,q}; b)\;
w_q |\det D_{m,q}| \det F_p.
$$

这里：

- `psi(...)`
  - 弹性能密度，由 `ElasticModel3DDeformationGradient` 提供
- `weightDetJ`
  - `w_q |\det D_m|`
- `detFp`
  - 塑性映射的体积因子

也就是说，这个实现的局部能量来自：

- 8 个高斯点的密度积分
- 每个高斯点都基于弹性部分 `Fe`

---

## 9. 对几何自由度的导数

### 9.1 一阶导 `compute_dE_dx`

代码：

```cpp
ES::M3d P;
ind->elasticModel->compute_P(..., cacheData->Fe[qi].data(), ..., P.data());
gradMap.noalias() += P * cacheData->Bm[qi];
```

对应理论上：

$$
\frac{\partial E_e}{\partial x}
=
\sum_q
\left(\frac{\partial F_e}{\partial x}\right)^T \mathrm{vec}(P_q).
$$

这里：

- $P_q = \partial \psi / \partial F_e$
  是第一类 Piola 应力
- `Bm` 是把应力映射回节点自由度的装配矩阵

这和标准非线性有限元里“应力积分得到内力”的形式完全一致。

### 9.2 二阶导 `compute_d2E_dx2`

代码：

```cpp
ES::M9d dPdF;
ind->elasticModel->compute_dPdF(..., dPdF.data());
dPdF *= ind->quad[qi].weightDetJ * cacheData->detFp;
hessMap.noalias() += cacheData->dFdx[qi].transpose() * dPdF * cacheData->dFdx[qi];
```

对应理论上：

$$
\frac{\partial^2 E_e}{\partial x^2}
=
\sum_q
\left(\frac{\partial F_e}{\partial x}\right)^T
\frac{\partial P_q}{\partial F_e}
\left(\frac{\partial F_e}{\partial x}\right).
$$

其中：

- `dPdF`
  是材料切线
- `dFdx`
  是几何部分的雅可比

这个表达式就是局部切线刚度矩阵。

---

## 10. 对塑性参数 `a` 的导数

这是 `CubicMeshDeformationModel` 和纯弹性六面体最大的区别。

### 10.1 基本链式法则

因为：

$$
F_e = F_{\mathrm{ref}} F_p^{-1},
$$

所以对某个塑性参数 $a_i$ 有：

$$
\frac{\partial F_e}{\partial a_i}
=
F_{\mathrm{ref}} \frac{\partial F_p^{-1}}{\partial a_i}.
$$

实现里对应：

```cpp
void CubicMeshDeformationModelInternal::compute_dFe_dai(const ES::M3d& Fref, const ES::M3d& dAInvdai,
                                                        ES::M3d& dFdai) const {
    dFdai = Fref * dAInvdai;
}
```

### 10.2 能量对 $a_i$ 的一阶导

代码中的核心式子：

```cpp
const double dVda = ...
const double dpsiDa = ...
gradMap[i] += dVda * psi + quad.weightDetJ * cacheData->detFp * dpsiDa;
```

对应理论：

$$
\frac{\partial E_e}{\partial a_i}
=
\sum_q
\left(
\frac{\partial V_q}{\partial a_i}\,\psi_q

+ V_q \frac{\partial \psi_q}{\partial a_i}
\right),
$$

其中：

$$
V_q = w_q |\det D_m| \det F_p.
$$

而：

$$
\frac{\partial \psi_q}{\partial a_i}
=
P_q : \frac{\partial F_e}{\partial a_i}.
$$

实现里这个双点积由：

```cpp
return P.cwiseProduct(dFe_dai).sum();
```

给出。

### 10.3 对 `a` 的二阶导和混合导

实现还支持：

- `compute_d2E_da2`
- `compute_d2E_dxda`

这两部分把：

- $\det(F_p)$ 对参数的导数
- $F_p^{-1}$ 对参数的导数
- $P$ 对 $F_e$ 的导数
- $F_e$ 对 $x$ 和 $a$ 的导数

全部通过链式法则展开。

这是它能支持塑性参数优化、耦合求解和灵敏度分析的关键。

---

## 11. 对材料参数 `b` 的导数

如果弹性模型有可调参数，例如某些参数化材料模型，那么这个类还支持对材料参数 `b` 求导。

### 一阶导

```cpp
gradMap[i] +=
    vol * ind->elasticModel->compute_dpsi_dparam(...);
```

对应：

$$
\frac{\partial E_e}{\partial b_i}
=
\sum_q V_q \frac{\partial \psi_q}{\partial b_i}.
$$

### 二阶导和混合导

类似地还有：

- `compute_d2E_db2`
- `compute_d2E_dxdb`
- `compute_d2E_dadb`

所以 `CubicMeshDeformationModel` 的接口不是只服务于常规位移求解，也兼容材料参数反演与联合优化。

---

## 12. 后处理：von Mises 应力与最大应变

### 12.1 von Mises 应力

代码：

```cpp
const double detF = cacheData->Fe[qi].determinant();
ES::M3d cauchyStress = P * cacheData->Fe[qi].transpose() / detF;
```

即先把一阶 Piola 应力转换成 Cauchy 应力：

$$
\sigma = \frac{1}{\det F_e} P F_e^T.
$$

然后再按经典公式计算 von Mises 标量。

### 12.2 最大应变

代码：

```cpp
ES::M3d E = 0.5 * (cacheData->Fe[qi].transpose() * cacheData->Fe[qi] - ES::M3d::Identity());
Eigen::SelfAdjointEigenSolver<ES::M3d> eigSolver(E);
stresses[qi] = eigSolver.eigenvalues().maxCoeff();
```

这里使用的是 Green-Lagrange strain：

$$
\mathcal{E} = \frac{1}{2}(F_e^T F_e - I),
$$

并取其最大特征值作为最大主应变指标。

---

## 13. 代码实现流程总览

可以把 `CubicMeshDeformationModel` 的执行流程压缩成下面这张表：

| 阶段 | 输入 | 主要输出 | 对应函数 |
| --- | --- | --- | --- |
| 构造期预计算 | `restPositions` | `QuadratureData[8]` | 构造函数 |
| 当前状态准备 | `x, a, b` | `Fref, Fe, dFdx, Bm, Fp, FpInv` | `prepareData` |
| 能量计算 | `CacheData` | `E_e` | `computeEnergy` |
| 力计算 | `CacheData` | `dE/dx` | `compute_dE_dx` |
| 刚度计算 | `CacheData` | `d²E/dx²` | `compute_d2E_dx2` |
| 塑性导数 | `CacheData` | `dE/da`, `d²E/da²`, `d²E/dxda` | `compute_dE_da` 等 |
| 材料导数 | `CacheData` | `dE/db`, `d²E/db²`, `d²E/dxdb`, `d²E/dadb` | `compute_dE_db` 等 |
| 后处理 | `CacheData` | 应力/应变标量 | `vonMisesStress`, `maxStrain` |

---

## 14. 和 `TetMeshDeformationModel` 的异同

和 `TetMeshDeformationModel` 相比，它们的核心框架很像：

- 都是 `DeformationModel`
- 都使用 deformation-gradient-based elastic / plastic model
- 都提供能量、一阶导、二阶导以及参数导数

但 `CubicMeshDeformationModel` 的特点是：

- 元素是 8 节点六面体，而不是 4 节点四面体
- shape function 是 trilinear，而不是线性 tetra basis
- 积分使用 8 个 Gauss 点，而不是 tetra 常见的低阶积分模板
- 需要在每个积分点单独计算 `Fref`、`Fe`、`dFdx`

因此，它更接近一个标准的非线性 hexahedral FEM 单元实现。

---

## 15. 阅读这份源码时最值得抓住的三条主线

如果你第一次读 `cubicMeshDeformationModel.cpp`，建议只盯住下面三条主线：

### 主线 1：几何主线

$$
X \rightarrow D_m \rightarrow D_m^{-1} \rightarrow F_{\mathrm{ref}}
$$

看懂参考几何、形函数梯度和 deformation gradient 的关系。

### 主线 2：弹塑性主线

$$
F_{\mathrm{ref}} \rightarrow F_p^{-1} \rightarrow F_e = F_{\mathrm{ref}} F_p^{-1}
$$

看懂为什么真正送进 elastic model 的不是 `Fref`，而是 `Fe`。

### 主线 3：求导主线

$$
E \rightarrow \frac{\partial E}{\partial x}, \frac{\partial^2 E}{\partial x^2},
\frac{\partial E}{\partial a}, \frac{\partial E}{\partial b}
$$

看懂它是如何把本构导数 `P`、`dP/dF` 与几何导数 `dF/dx` 组合起来的。

---

## 16. 相关文件

建议结合下面几个文件一起看：

- [deformationModel.md](deformationModel.md)
  - `DeformationModel` 抽象接口
- [deformationModelManager.md](deformationModelManager.md)
  - `CubicMeshDeformationModel` 是如何被创建的
- [deformationModelAssembler.md](deformationModelAssembler.md)
  - 单元级量如何装配成全局量
- [plasticModel3DDeformationGradient.md](plasticModel3DDeformationGradient.md)
  - `F_p` 参数化与导数接口
- [elasticModel3DDeformationGradient.md](elasticModel3DDeformationGradient.md)
  - `psi`, `P`, `dP/dF` 的上层接口
- [simulationMesh.md](simulationMesh.md)
  - cubic 单元的 8 个顶点从哪里来

---

## 17. 一句话总结

`CubicMeshDeformationModel` 可以概括为：

**一个基于 8 节点 trilinear 六面体、8 点高斯积分、并支持 deformation-gradient 弹塑性本构及其高阶导数的局部单元模型。**

它是整个 solid deformation pipeline 中“局部能量和局部导数”最核心的计算单元之一。
