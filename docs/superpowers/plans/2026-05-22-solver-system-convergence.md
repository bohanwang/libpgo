# Solver System Convergence Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [x]`) syntax for tracking.

**Goal:** 将 solver 体系从混合 `int` 返回码收口到统一的 `SolverResult` / `SolveStatus` API，并清理外部 solver 宏、状态映射和 timestep 接受策略。

**Architecture:** `NewtonSolver` 已经是 typed-result API；本计划继续把 `EnergyOptimizer::minimize*`、`TimeIntegratorSolver::solveDirect`、legacy callers 和日志全部迁到同一结果对象。外部 solver 的原始返回码保留在 `rawStatusCode`，语义状态通过集中映射函数转换，dynamic/static 的接受策略通过公共 helper 表达。

**Tech Stack:** C++17, CMake, Eigen, optional IPOPT/Knitro, GoogleTest.

---

## 文件结构与职责

- `src/core/nonlinearOptimization/solverResult.h/.cpp`
  - 增加 solver 结果构造、IPOPT/Knitro raw status 映射、dynamic/static 接受策略 helper。
- `src/core/nonlinearOptimization/minimizeEnergy.h/.cpp`
  - 将 `EnergyOptimizer::minimize*` 系列返回值从 `int` 改为 `SolverResult`。
  - `minimizeUsingNewton` 保留 Newton 的完整 diagnostics，不再只返回 raw code。
  - IPOPT/Knitro wrapper 在拿到 raw code 后立即调用集中映射函数。
- `src/core/simulation/timeIntegratorSolver.h/.cpp`
  - 移除本地 `makeExternalSolverResult`，复用 `EnergyOptimizer` 和 `solverResult` 的统一 helper。
  - `SO_IPOPT` 使用 `PGO_HAS_IPOPT`，不再使用 `USE_IPOPT`。
  - `solveDirect` 变成 thin wrapper，不再手动桥接 `int`。
- `src/core/simulation/implicitBackwardEulerTimeIntegrator.cpp`
  - 将 dynamic timestep 接受策略替换为公共 helper。
- `src/core/genericPotentialEnergies/laplacianProblem.cpp`
  - 调用 typed `minimizeUsingKnitro`，日志输出 `formatSolverResultSummary`。
- `src/core/nonlinearOptimization/CMakeLists.txt`, `CMakeLists.txt`
  - 统一 optional solver 宏为 `PGO_HAS_IPOPT` / `PGO_HAS_KNITRO`。
  - 修正 IPOPT imported target 使用 `Ipopt::Core`。
- `tests/src/core/NewtonSolver_gtest.cpp`
  - 增加 `solverResult` helper 的单元测试，覆盖 raw mapping 和 timestep 接受策略。

## API 目标

目标签名：

```cpp
SolverResult minimize(...);
SolverResult minimizeUsingIpopt(...);
SolverResult minimizeUsingNewton(...);
SolverResult minimizeUsingKnitro(...);
SolverResult minimizeUsingKnitroDense(...);
SolverResult minimizeUsingApproximateActiveSet(...);
```

保留语义：

- `SolverResult::status` 是统一语义。
- `SolverResult::rawStatusCode` 保存原始 backend code。
- `SolverResult::diagnostics` 对 Newton 有完整信息；外部 solver 默认 reset。
- static solve 只接受 `Converged`。
- dynamic timestep 接受 `Converged`、`MaxIterations`、`StepTooSmall`。

## Task 1: 扩展 `SolverResult` 公共 helper

**Files:**
- Modify: `src/core/nonlinearOptimization/solverResult.h`
- Modify: `src/core/nonlinearOptimization/solverResult.cpp`
- Modify: `tests/src/core/NewtonSolver_gtest.cpp`

- [x] **Step 1: 增加声明**

在 `solverResult.h` 增加：

```cpp
SolverResult makeSolverResult(SolveStatus status, int rawStatusCode, int iterations = 0);
SolverResult makeIpoptSolverResult(int rawStatusCode);
SolverResult makeKnitroSolverResult(int rawStatusCode);
bool acceptsDynamicSolveStatus(SolveStatus status);
bool acceptsStrictSolveStatus(SolveStatus status);
```

- [x] **Step 2: 增加实现**

在 `solverResult.cpp` 中实现：

```cpp
SolverResult makeSolverResult(SolveStatus status, int rawStatusCode, int iterations)
{
  SolverResult result;
  result.status = status;
  result.rawStatusCode = rawStatusCode;
  result.iterations = iterations;
  result.diagnostics.reset();
  return result;
}

SolverResult makeIpoptSolverResult(int rawStatusCode)
{
  switch (rawStatusCode) {
    case 0:
    case 1:
    case 6:
      return makeSolverResult(SolveStatus::Converged, rawStatusCode);
    case -1:
      return makeSolverResult(SolveStatus::MaxIterations, rawStatusCode);
    case 3:
      return makeSolverResult(SolveStatus::StepTooSmall, rawStatusCode);
    case -13:
      return makeSolverResult(SolveStatus::NonFinite, rawStatusCode);
    case -3:
      return makeSolverResult(SolveStatus::LinearSolveFailed, rawStatusCode);
    default:
      return makeSolverResult(SolveStatus::ExternalSolverFailure, rawStatusCode);
  }
}

SolverResult makeKnitroSolverResult(int rawStatusCode)
{
  if (rawStatusCode == 0 || rawStatusCode == -100 || rawStatusCode == -101 || rawStatusCode == -102) {
    return makeSolverResult(SolveStatus::Converged, rawStatusCode);
  }
  if (rawStatusCode == -400 || rawStatusCode == -401 || rawStatusCode == -402) {
    return makeSolverResult(SolveStatus::MaxIterations, rawStatusCode);
  }
  if (rawStatusCode <= -500 && rawStatusCode > -600) {
    return makeSolverResult(SolveStatus::LinearSolveFailed, rawStatusCode);
  }
  return makeSolverResult(SolveStatus::ExternalSolverFailure, rawStatusCode);
}

bool acceptsDynamicSolveStatus(SolveStatus status)
{
  return status == SolveStatus::Converged ||
    status == SolveStatus::MaxIterations ||
    status == SolveStatus::StepTooSmall;
}

bool acceptsStrictSolveStatus(SolveStatus status)
{
  return status == SolveStatus::Converged;
}
```

- [x] **Step 3: 增加测试**

在 `NewtonSolver_gtest.cpp` 的 status string 测试附近增加：

```cpp
TEST(SolverResultGTest, MapsExternalRawStatusCodes)
{
  EXPECT_EQ(makeIpoptSolverResult(0).status, SolveStatus::Converged);
  EXPECT_EQ(makeIpoptSolverResult(1).status, SolveStatus::Converged);
  EXPECT_EQ(makeIpoptSolverResult(-1).status, SolveStatus::MaxIterations);
  EXPECT_EQ(makeIpoptSolverResult(3).status, SolveStatus::StepTooSmall);
  EXPECT_EQ(makeIpoptSolverResult(-13).status, SolveStatus::NonFinite);
  EXPECT_EQ(makeIpoptSolverResult(-3).status, SolveStatus::LinearSolveFailed);
  EXPECT_EQ(makeIpoptSolverResult(-199).status, SolveStatus::ExternalSolverFailure);
  EXPECT_EQ(makeIpoptSolverResult(-199).rawStatusCode, -199);

  EXPECT_EQ(makeKnitroSolverResult(0).status, SolveStatus::Converged);
  EXPECT_EQ(makeKnitroSolverResult(-100).status, SolveStatus::Converged);
  EXPECT_EQ(makeKnitroSolverResult(-400).status, SolveStatus::MaxIterations);
  EXPECT_EQ(makeKnitroSolverResult(-500).status, SolveStatus::LinearSolveFailed);
  EXPECT_EQ(makeKnitroSolverResult(-200).status, SolveStatus::ExternalSolverFailure);
  EXPECT_EQ(makeKnitroSolverResult(-200).rawStatusCode, -200);
}

TEST(SolverResultGTest, EncodesDynamicAndStrictAcceptancePolicies)
{
  EXPECT_TRUE(acceptsStrictSolveStatus(SolveStatus::Converged));
  EXPECT_FALSE(acceptsStrictSolveStatus(SolveStatus::MaxIterations));

  EXPECT_TRUE(acceptsDynamicSolveStatus(SolveStatus::Converged));
  EXPECT_TRUE(acceptsDynamicSolveStatus(SolveStatus::MaxIterations));
  EXPECT_TRUE(acceptsDynamicSolveStatus(SolveStatus::StepTooSmall));
  EXPECT_FALSE(acceptsDynamicSolveStatus(SolveStatus::LineSearchFailed));
  EXPECT_FALSE(acceptsDynamicSolveStatus(SolveStatus::ExternalSolverFailure));
}
```

- [x] **Step 4: 验证**

Run:

```bash
cmake --build --preset base_no_mkl_release --target NewtonSolver_gtest
./build/base_no_mkl/tests/src/core/NewtonSolver_gtest
```

Expected: all `NewtonSolver_gtest` tests pass.

## Task 2: 将 `EnergyOptimizer` API 迁到 `SolverResult`

**Files:**
- Modify: `src/core/nonlinearOptimization/minimizeEnergy.h`
- Modify: `src/core/nonlinearOptimization/minimizeEnergy.cpp`
- Modify: `src/core/genericPotentialEnergies/laplacianProblem.cpp`

- [x] **Step 1: 头文件引入 typed result**

在 `minimizeEnergy.h` 增加：

```cpp
#include "solverResult.h"
```

并将 `EnergyOptimizer::minimize*` 系列返回值从 `int` 改为 `SolverResult`。

- [x] **Step 2: 实现文件统一 optional solver 宏**

将 `minimizeEnergy.cpp` 中所有 `USE_IPOPT` 改为 `PGO_HAS_IPOPT`，所有 `USE_KNITRO` 改为 `PGO_HAS_KNITRO`。

- [x] **Step 3: 修改 `minimizeUsingIpopt`**

每次 `solver.solve()` 后返回：

```cpp
return makeIpoptSolverResult(solverRet);
```

缺失 IPOPT 分支继续 `throw std::runtime_error("No Ipopt Module");`，不再返回 magic `1`。

- [x] **Step 4: 修改 `minimizeUsingKnitro` 和 dense 版本**

每次 `solver.solve()` 后返回：

```cpp
return makeKnitroSolverResult(solverRet);
```

Knitro 不存在时，volume 版本沿用旧 fallback 到 IPOPT，但返回 `SolverResult`：

```cpp
return minimizeUsingIpopt(...);
```

dense 版本缺失 Knitro 时改为：

```cpp
throw std::runtime_error("No Knitro Module");
```

- [x] **Step 5: 修改 `minimizeUsingApproximateActiveSet`**

内部 lambda `minf` 返回 `SolverResult`；原先 `solverRet < 0` 判断改为 `!acceptsStrictSolveStatus(result.status)`。
最终返回最后一次 `SolverResult`。

- [x] **Step 6: 修改 `minimizeUsingNewton`**

保留 `SolverResult ret = solver.solve(...)`，函数直接返回 `ret`，不再返回 `ret.rawStatusCode`。

- [x] **Step 7: 修改调用方**

`laplacianProblem.cpp` 改为：

```cpp
NonlinearOptimization::SolverResult result = NonlinearOptimization::EnergyOptimizer::minimizeUsingKnitro(...);
std::cout << "Solver result: " << NonlinearOptimization::formatSolverResultSummary(result) << std::endl;
```

- [x] **Step 8: 验证**

Run:

```bash
cmake --build --preset base_no_mkl_release --target nonlinearOptimization
```

Expected: `nonlinearOptimization` target builds.

## Task 3: 收口 `TimeIntegratorSolver`

**Files:**
- Modify: `src/core/simulation/timeIntegratorSolver.h`
- Modify: `src/core/simulation/timeIntegratorSolver.cpp`
- Modify: `src/core/simulation/implicitBackwardEulerTimeIntegrator.cpp`
- Modify: `tests/src/core/implicitBackwardEulerTimeIntegrator_gtest.cpp`

- [x] **Step 1: 删除本地 `makeExternalSolverResult`**

从 `timeIntegratorSolver.cpp` 删除匿名 namespace 中的 local helper。

- [x] **Step 2: SO_KNITRO 直接使用集中映射**

将：

```cpp
da->lastSolveResult = makeExternalSolverResult(solverRet);
```

改为：

```cpp
da->lastSolveResult = makeKnitroSolverResult(solverRet);
```

- [x] **Step 3: SO_IPOPT 宏改为 `PGO_HAS_IPOPT`**

将 `#if defined(USE_IPOPT)` 改为 `#if defined(PGO_HAS_IPOPT)`。

- [x] **Step 4: `solveDirect` 返回 typed result**

将 `int solverRet = 0;` 改为：

```cpp
SolverResult result = makeSolverResult(SolveStatus::UnsupportedBackend,
  static_cast<int>(SolveStatus::UnsupportedBackend));
```

每个分支直接赋值 `EnergyOptimizer::minimize(...)` 或 `EnergyOptimizer::minimizeUsingKnitro(...)` 的返回值，最后 `return result;`。

- [x] **Step 5: dynamic 接受策略使用公共 helper**

删除 `implicitBackwardEulerTimeIntegrator.cpp` 中本地 `acceptsDynamicSolveStatus`，直接使用 `NonlinearOptimization::acceptsDynamicSolveStatus`。

- [x] **Step 6: 验证**

Run:

```bash
cmake --build --preset base_no_mkl_release --target implicitBackwardEulerTimeIntegrator_gtest
./build/base_no_mkl/tests/src/core/implicitBackwardEulerTimeIntegrator_gtest
```

Expected: all tests pass.

## Task 4: 统一 optional solver CMake 宏和 target

**Files:**
- Modify: `CMakeLists.txt`
- Modify: `src/core/nonlinearOptimization/CMakeLists.txt`

- [x] **Step 1: 启用 IPOPT option 的发现逻辑**

在 root `CMakeLists.txt` 的 Knitro block 前增加：

```cmake
if(PGO_OPT_USE_IPOPT)
  find_package(Ipopt)

  if(TARGET Ipopt::Core)
    add_compile_definitions(PGO_HAS_IPOPT)
  endif()
endif()
```

- [x] **Step 2: 修正 nonlinearOptimization 的 IPOPT target**

将 `TARGET Ipopt::Ipopt` 全部改为 `TARGET Ipopt::Core`，link 也使用 `Ipopt::Core`。

- [x] **Step 3: 宏残留扫描**

Run:

```bash
rg -n "USE_IPOPT|USE_KNITRO|Ipopt::Ipopt|PGO_HAS_IPOPT|PGO_HAS_KNITRO" CMakeLists.txt CMakeModules src tests
```

Expected:
- No `USE_IPOPT`
- No `USE_KNITRO`
- No `Ipopt::Ipopt`
- `PGO_HAS_IPOPT` and `PGO_HAS_KNITRO` only appear in intentional feature gates.

## Task 5: 回归验证与文档

**Files:**
- Modify: `README.md`
- Modify: `examples/ipc/README.md`

- [x] **Step 1: 文档补充**

在 solver status 说明附近补充：`EnergyOptimizer`、`TimeIntegratorSolver`、`NewtonSolver` 均返回/保存 `SolverResult`；外部 solver raw code 在 `rawStatusCode` 中保留。

- [x] **Step 2: API 残留扫描**

Run:

```bash
rg -n "int minimize|int minimizeUsing|USE_IPOPT|USE_KNITRO|makeExternalSolverResult|getSolverReturn\\(|solverRet=|Solver Ret" src tests README.md examples/ipc/README.md
```

Expected:
- No old `int minimize*` declarations.
- No `USE_*` optional solver macros.
- No old timestep return logs.

- [x] **Step 3: 全量聚焦验证**

Run:

```bash
cmake --build --preset base_no_mkl_release --target NewtonSolver_gtest implicitBackwardEulerTimeIntegrator_gtest runIPCSim_gtest runSimShared_gtest
./build/base_no_mkl/tests/src/core/NewtonSolver_gtest
./build/base_no_mkl/tests/src/core/implicitBackwardEulerTimeIntegrator_gtest
./build/base_no_mkl/tests/src/tools/runSimShared_gtest
./build/base_no_mkl/tests/src/tools/runIPCSim_gtest
git diff --check
```

Expected:
- All listed tests pass.
- `git diff --check` emits no output.

## Self Review

### 1. 需求覆盖

- “solver 体系全部收口”：Task 2 改 `EnergyOptimizer` API，Task 3 改 `TimeIntegratorSolver`，Task 1 集中 result helper，覆盖核心链路。
- “外部 solver 状态映射”：Task 1 增加 IPOPT/Knitro mapper，Task 2/3 使用 mapper。
- “宏体系统一”：Task 4 明确 `PGO_HAS_*` 和 `Ipopt::Core`。
- “dynamic/static 策略”：Task 1 定义公共 helper；Task 3 迁 dynamic；static 已经使用 strict `Converged` 语义，helper 可供后续调用。
- “验证正确性”：Task 1/3/5 有构建、单测、入口回归和 diff check。

### 2. 占位符扫描

计划中没有 `TBD`、`TODO`、`implement later`，每个代码改动任务都列了具体文件、目标签名、命令和期望结果。

### 3. 类型一致性

- 所有 `minimize*` 返回类型统一为 `SolverResult`。
- `TimeIntegratorSolver::solveDirect` 已经是 `SolverResult`，不会再做 `int` bridge。
- `LaplacianProblem` 从 `int ret` 迁到 `SolverResult result`。
- `acceptsDynamicSolveStatus` 命名与现有 IBE 本地函数一致，迁移成本低。

### 4. 风险与调整

- IPOPT target 目前 CMake 中存在 `FindIpopt.cmake` 创建 `Ipopt::Core`，但 nonlinearOptimization 检查的是 `Ipopt::Ipopt`；Task 4 明确修复这个不一致。
- Knitro raw code 映射基于常见 Knitro return code 区间；保留 `rawStatusCode`，避免损失诊断信息。
- `minimizeUsingApproximateActiveSet` 原先用负数判断错误，现在改用 typed status 后语义更清楚；这块风险中等，需要编译覆盖。
- 本计划不引入兼容 `int` wrapper，符合“不向前兼容，直接新 API”的决策。
