# Phase 1：配置层与 Runtime 分发

## 1. 这阶段要回答什么

这一阶段专门回答下面四个问题：

1. `ipc-barrier` 是怎样作为一种 contact model 进入 `RunSimConfig` 的。
2. external / self 为什么在当前实现里使用拆开的 `dhat` 和 `kappa`。
3. parser 现在对 IPC 模式做了哪些显式限制和默认行为。
4. dynamic runtime 怎样根据配置决定某一帧是否创建 contact handler、barrier energy 和 feasible alpha callback。

如果只看“多了几个 JSON key”，会低估这一阶段的作用。  
Phase 1 真正完成的是：

> 把 IPC 从一组零散参数，提升成 dynamic runtime 里有明确语义的一级分发开关。

## 2. 这阶段在整条 IPC 主线里解决什么

后面几个 phase 会分别讲：

- external active sample 怎样构造、怎样变成 barrier energy
- self active pair 怎样从 BVH query band 细化出来
- external/self barrier 怎样进入 timestep 总势能
- feasible alpha 怎样进入 Newton

但这些阶段共享同一套入口前提。  
如果 Phase 1 没把下面三件事先固定下来，后面的 runtime 语义都会发散：

- 当前是 penalty 还是 barrier
- external/self 是否允许独立启停
- solver-side feasibility filter 是否启用

因此 Phase 1 在当前 repo 里不只是 parser 层增量，而是整条 IPC 主线的 dispatch phase。

## 3. 当前 repo 已经完成了什么

当前 `RunSimContactConfig` 已经不再只有 penalty 接触参数；它包含：

- `contactModel`
- `contactStiffness`
- `contactSamples`
- `enableSelfContact`
- `contactFrictionCoeff`
- `contactVelEps`
- `externalIpcDhat`
- `selfIpcDhat`
- `externalIpcKappa`
- `selfIpcKappa`
- `ipcAlphaSafety`
- `ipcEnableFeasibleLineSearch`

同时，当前仓库已经把下面几件事固定成 runtime 事实：

- `parseRunSimConfig(...)` 支持 `"ipc-barrier"`
- external / self IPC 参数是解耦表达，而不是共享一组 `ipc-dhat/ipc-kappa`
- friction 在 IPC 模式下被显式禁止
- dynamic runtime 不再只看 `contactStiffness > 0`
- 每帧 contact runtime 都会先 `clearGeneralImplicitForceModel()` 与 `clearAlphaTestFunc()`

这说明当前 repo 里的 IPC 配置已经不只是“可写进 JSON”，而是已经真正控制 runtime 行为。

## 4. 当前配置结构里新增了什么

### 4.1 `RunSimContactConfig` 已经显式区分 penalty 和 IPC

当前 `src/api/runSimCore.h` 中的 contact 配置结构是：

```cpp
struct RunSimContactConfig {
    std::string contactModel         = "penalty";
    double      contactStiffness     = 0.0;
    int         contactSamples       = 1;
    bool        enableSelfContact    = true;
    double      contactFrictionCoeff = 0.0;
    double      contactVelEps        = 0.0;
    double      externalIpcDhat      = 1e-3;
    double      selfIpcDhat          = 1e-3;
    double      externalIpcKappa     = 1e4;
    double      selfIpcKappa         = 1e4;
    double      ipcAlphaSafety       = 0.9;
    bool        ipcEnableFeasibleLineSearch = true;
};
```

这里最重要的不是字段数量变多，而是 contact 语义已经拆成两层：

- 通用 contact 配置
  例如 `contactSamples`、`enableSelfContact`
- barrier 专属配置
  例如 external/self `dhat`、external/self `kappa`、alpha filter 开关

这样后面的 runtime 可以共享 sampling 与 self-contact enable 语义，但不必再假设 penalty 和 IPC 共用同一套数值参数。

### 4.2 external / self `dhat` 与 `kappa` 被拆开表达

当前实现没有使用单一 `ipc-dhat` 或单一 `ipc-kappa`，而是拆成：

- `external-ipc-dhat`
- `self-ipc-dhat`
- `external-ipc-kappa`
- `self-ipc-kappa`

这件事在当前 repo 里很关键，因为 external 和 self 的 runtime 结构并不对称：

- external 走 sampled point vs fixed external target snapshot
- self 走 sampled point-triangle pair

把两条路径的 barrier 参数拆开，意味着 runtime 不需要再假设 external 与 self 必须同尺度调参。

## 5. `parseRunSimConfig(...)` 当前怎样固定 IPC 语义

### 5.1 parser 先读通用 contact 字段，再读 `contact-model`

当前 parser 会先读：

- `contact-stiffness`
- `contact-sample`
- `enable-self-contact`
- `contact-friction-coeff`
- `contact-vel-eps`

然后再用 `jconfig.handle().value(...)` 读取：

```cpp
config.contact.contactModel = jconfig.handle().value("contact-model", std::string("penalty"));
```

这意味着当前 contact model 的默认行为是：

- 配置里没写 `contact-model`
  仍然走 `penalty`
- 配置里写了 `contact-model = "ipc-barrier"`
  才进入 IPC 专属解析分支

### 5.2 parser 对 contact model 做显式校验

当前实现只接受两种 contact model：

- `"penalty"`
- `"ipc-barrier"`

任何其它字符串都会直接抛错。  
因此在当前 repo 里，IPC 不是隐式模式切换，而是明确的配置选择。

### 5.3 IPC 专属字段只在 barrier 模式下解析

只有当 `contact-model == "ipc-barrier"` 时，parser 才会继续读取：

- `external-ipc-dhat`
- `self-ipc-dhat`
- `external-ipc-kappa`
- `self-ipc-kappa`
- `ipc-alpha-safety`
- `ipc-enable-feasible-line-search`

这一步把当前 repo 的 contact 语义固定成：

- penalty 模式
  继续按旧 penalty contact 配置运行
- IPC 模式
  进入 frictionless、decoupled external/self barrier runtime

### 5.4 旧共享键已经被显式拒绝

当前 parser 会主动拒绝：

- `ipc-dhat`
- `ipc-kappa`

并要求使用 external/self 分开的新键。  
这意味着仓库当前已经不允许把新 runtime 继续描述成“一组共享的 IPC 参数”。

### 5.5 IPC 模式下 friction 被显式禁止

当前 parser 还固定了一条重要限制：

- `contact-model == "ipc-barrier"`
- 且 `contact-friction-coeff > 0`

会直接抛错。

所以仓库当前对 IPC 的真实表述应当是：

> frictionless barrier contact only

而不是“已经有 frictional IPC，只是文档没写”。

### 5.6 默认值已经被写入 parser 语义

IPC 模式下当前默认值包括：

- `external-ipc-dhat = 1e-3`
- `self-ipc-dhat = 1e-3`
- `external-ipc-kappa = 1e4`
- `self-ipc-kappa = 1e4`
- `ipc-alpha-safety = 0.9`
- `ipc-enable-feasible-line-search = true`

这意味着 feasible alpha filter 在当前设计中不是可有可无的附属项，而是 IPC runtime 的默认组成部分。

## 6. dynamic runtime 当前怎样根据配置做分发

### 6.1 `isDynamicContactEnabled(...)` 不再只看 penalty stiffness

当前 `runSimCore.cpp` 使用下面这个判据决定 dynamic path 是否需要 contact runtime：

```cpp
bool isDynamicContactEnabled(const RunSimContactConfig& contact) {
    return (contact.contactModel == "penalty" && contact.contactStiffness > 0.0) ||
           (contact.contactModel == "ipc-barrier" &&
            (contact.externalIpcKappa > 0.0 || contact.selfIpcKappa > 0.0));
}
```

这一步非常关键，因为它让两条 contact 模式在 runtime 入口上真正分家：

- penalty 由 `contactStiffness` 驱动
- IPC 由 external/self barrier scale 驱动

也就是说，当前 repo 已经不把 `contactStiffness` 当作 barrier 的别名。

### 6.2 handler 的构建条件也跟着 contact model 一起切换

在 dynamic path 中，external handler 与 self handler 只会在 `contactEnabled` 成立时创建。  
然后再叠加：

- external 还需要存在 `external-objects`
- self 还需要 `enable-self-contact = true`

因此 Phase 1 不只是 parser 层逻辑，还决定了 runtime 是否会真正持有：

- `TriangleMeshExternalContactHandler`
- `TriangleMeshSelfContactHandler`

### 6.3 每帧 contact runtime 是按 frame 重建的

dynamic 帧循环开始时，当前实现固定先做：

```cpp
intg->clearGeneralImplicitForceModel();
intg->clearAlphaTestFunc();
```

这个顺序很重要。它说明当前 IPC runtime 的生命周期是：

- 配置是跨帧持有的
- active set、barrier object、alpha callback 是每帧重建的

因此 Phase 1 里的 dispatch 不是一次性初始化，而是：

> 用全局配置控制每一帧 contact runtime 应当怎样重建

### 6.4 Phase 1 实际同时控制三条后续链路

当 `contact-model == "ipc-barrier"` 时，Phase 1 实际上已经决定了下面三类后续行为：

- external path 用 `externalIpcDhat` / `externalIpcKappa`
- self path 用 `selfIpcDhat` / `selfIpcKappa`
- solver callback 用 `ipcAlphaSafety` / `ipcEnableFeasibleLineSearch`

也就是说，后面的 active-set、energy、dynamic-assembly 和 line-search 几篇虽然会分开讲，但它们的 runtime 开关都已经在这里固定。

## 7. example 侧当前怎样使用这些配置

当前仓库里最直接的两个 IPC 配置入口是：

- `examples/cubic-box-ipc/cubic-box-ipc-ls.json`
- `examples/pulled-cubic-box-self-ipc/pulled-cubic-box-self-ipc.json`

前者主要体现：

- external IPC
- `enable-self-contact = false`
- `contact-sample = 6`
- 更小的 `external-ipc-dhat`
- 更高的 `external-ipc-kappa`

后者主要体现：

- self IPC 与 merged path
- `enable-self-contact = true`
- `external-ipc-dhat` / `self-ipc-dhat` 同时启用
- `external-ipc-kappa` / `self-ipc-kappa` 同时启用
- `ipc-enable-feasible-line-search = true`

也就是说，example 目录已经直接把 Phase 1 的配置语义落成可运行 scene。

## 8. 当前 repo 的验证证据

`tests/api/runSimConfig_parse_test.cpp` 已经把这阶段最关键的行为写成 committed coverage，包括：

- `DefaultsToPenaltyContactModelAndDefaultIpcParameters`
- `ParsesIpcBarrierContactParameters`
- `ParsesDecoupledExternalAndSelfIpcBarrierParameters`
- `RejectsLegacySharedIpcBarrierParameters`
- `IpcBarrierFeasibleLineSearchDefaultsToEnabled`
- `RejectsUnsupportedContactModel`
- `RejectsNonPositiveExternalIpcDhat`
- `RejectsNonPositiveSelfIpcKappa`
- `RejectsOutOfRangeIpcAlphaSafety`
- `RejectsFrictionWithIpcBarrierContactModel`

这些测试共同固定了下面这组 runtime 事实：

- IPC 参数有默认值
- external/self 参数分离
- 旧共享配置被拒绝
- feasible line search 默认开启
- friction 在 IPC 模式下被禁止

## 9. 这一阶段不展开什么

这一阶段还不展开下面这些内容：

- external active sample 怎样生成
- self active pair 怎样生成
- barrier energy 的 `func / gradient / hessian`
- feasible alpha upper bound 的几何计算
- `ImplicitBackwardEuler` 内部怎样消费 general `PotentialEnergy`

这些分别留给：

- [Phase 2](phase-2-external-barrier-energy.md)
- [Phase 3](phase-3-self-near-contact-active-set.md)
- [Phase 4](phase-4-dynamic-incremental-potential.md)
- [Phase 5](phase-5-feasible-line-search.md)

上一阶段： [IPC Overview](index.md)  
下一阶段： [Phase 2](phase-2-external-barrier-energy.md)
