# Phase 5：Feasible Line Search

仅有 barrier energy 还不够 — solver 仍可能沿让接触距离恶化的方向走 full step。本篇描述 contact-aware feasible alpha upper bound 如何在不重写 Newton 的前提下接入 solver。

## Callback 安装条件

两层检查：

1. **配置层**：`contact-model == "ipc-barrier"` 且 `ipc-enable-feasible-line-search == true`
2. **帧级别**：当前帧至少有一条路径存在 active set（`hasExternalAlphaUpperBound || hasSelfAlphaUpperBound`）

没有 active set 的帧不安装 callback。

## Merged Callback 组装

`runSimCore.cpp` 中安装的 callback 流程：

```cpp
intg->setAlphaTestFunc([...](const ES::VXd& z, const ES::VXd& dz) -> double {
    ES::VXd currentU = u + z;      // 当前绝对位移
    // currentU 和 dz 用于查询 handler
    double alphaUpper = 1.0;
    if (hasExternalAlphaUpperBound)
        alphaUpper = min(alphaUpper, externalHandler->computeEmbeddedAlphaUpperBound(
            currentU, dz, externalDSafe, alphaSafety));
    if (hasSelfAlphaUpperBound)
        alphaUpper = min(alphaUpper, selfCD->computeEmbeddedAlphaUpperBound(
            currentU, dz, selfDSafe, alphaSafety));
    return alphaUpper;
});
```

关键点：

- Callback 输入是 Newton 内部变量 `z`（增量）和 `dz`（搜索方向），需转换成绝对位移
- `dSafe` = `normalizePositiveZero(dhat, -1.0)`，与 barrier 数值安全域对齐
- External / self 上界取 min — solver 不需要理解两条路径的语义

## `ipc-alpha-safety` 的作用

只影响 alpha upper bound 的最终结果（乘以 `clamp(alphaSafety, 0, 1)`）。不参与 active set、不参与 barrier energy、不参与 `kappa`。

- `alphaSafety = 1.0`：使用 conservative 几何上界
- `alphaSafety < 1.0`：在上界基础上进一步收紧

## Solver 侧消费

透传链路：`TimeIntegrator::setAlphaTestFunc(...)` → `TimeIntegratorSolver::solve(...)` → `NewtonRaphsonSolver`。

`NewtonRaphsonSolver::solve(...)` 中的关键逻辑：

```cpp
alpha = alphaTestFunc ? alphaTestFunc(x, deltax) : 1.0;

if (alpha <= 0.0) {
    status = SR_FEASIBLE_STEP_ZERO;
    break;
}
// 否则从 alpha 开始 LSM_SIMPLE step acceptance（alpha *= 0.75 回缩）
```

`SR_FEASIBLE_STEP_ZERO` 是可观测 solver 状态（`statusToString` → `"feasible_step_zero"`）。

默认主路径是 `SO_NEWTON` + `LSM_SIMPLE`：先用 contact callback 给出几何上界，再在上界内做 step acceptance。

## 验证

- External alpha：`contact_embedding_test.cpp` — 解析 plane bound、safety scaling、outward recovery、safe-band 行为
- Self alpha：`self_contact_handler_test.cpp` — inward shrink、safety scaling、outward recovery、safe-band block
- Runtime：API smoke 验证日志中存在 `IPC feasible alpha callback active.` 和 `< 1` 的 alpha upper bound

---

上一阶段：[Phase 4](phase-4-dynamic-incremental-potential.md)
下一阶段：[Phase 6](phase-6-tests-and-validation.md)
