# Phase 1：配置层与 Runtime 分发

Phase 1 把 IPC 从一组零散参数提升成 dynamic runtime 的一级分发开关。后续所有 phase 的 runtime 行为（active-set、energy、line-search）都由这里的配置决定。

## `RunSimContactConfig` 结构

```cpp
struct RunSimContactConfig {
    std::string contactModel         = "penalty";
    double      contactStiffness     = 0.0;
    int         contactSamples       = 1;
    bool        enableSelfContact    = true;
    double      contactFrictionCoeff = 0.0;
    double      contactVelEps        = 0.0;
    // --- IPC barrier 专属 ---
    double      externalIpcDhat      = 1e-3;
    double      selfIpcDhat          = 1e-3;
    double      externalIpcKappa     = 1e4;
    double      selfIpcKappa         = 1e4;
    double      ipcAlphaSafety       = 0.9;
    bool        ipcEnableFeasibleLineSearch = true;
};
```

关键设计：external / self 的 `dhat` 和 `kappa` 分开表达，因为两条路径的几何结构不对称（external 是 point vs fixed plane，self 是 point-triangle pair），不应共享调参尺度。

## 配置解析逻辑

`parseRunSimConfig(...)` 的 IPC 分支：

1. **读取 `contact-model`**：只接受 `"penalty"` 或 `"ipc-barrier"`，其它值直接抛错
2. **拒绝旧共享键**：`ipc-dhat` / `ipc-kappa` 会触发错误，要求使用 external/self 分离的新键
3. **IPC 专属字段仅在 barrier 模式下解析**：`external-ipc-dhat`、`self-ipc-dhat`、`external-ipc-kappa`、`self-ipc-kappa`、`ipc-alpha-safety`、`ipc-enable-feasible-line-search`
4. **Friction 禁止**：`ipc-barrier` + `contact-friction-coeff > 0` → 抛错
5. **参数校验**：`dhat` 必须正、`kappa` 必须正、`alpha-safety` 必须在 `(0, 1]`

## Runtime 分发

### `isDynamicContactEnabled(...)`

```cpp
bool isDynamicContactEnabled(const RunSimContactConfig& contact) {
    return (contact.contactModel == "penalty" && contact.contactStiffness > 0.0) ||
           (contact.contactModel == "ipc-barrier" &&
            (contact.externalIpcKappa > 0.0 || contact.selfIpcKappa > 0.0));
}
```

两种 contact model 在 runtime 入口真正分家：penalty 由 `contactStiffness` 驱动，IPC 由 barrier scale 驱动。

### Handler 创建条件

`contactEnabled` 成立后，还需叠加：

- external handler：需要存在 `external-objects`
- self handler：需要 `enable-self-contact = true`

### 每帧重建协议

帧循环开始时固定先执行：

```cpp
intg->clearGeneralImplicitForceModel();
intg->clearAlphaTestFunc();
```

配置是跨帧持有的，但 active set、barrier object、alpha callback 都是**每帧重建**。Phase 1 的 dispatch 不是一次性初始化，而是用全局配置控制每帧 contact runtime 的重建方式。

## Example 配置

- `examples/cubic-box-ipc/cubic-box-ipc-ls.json`：external-only IPC，`enable-self-contact = false`
- `examples/pulled-cubic-box-self-ipc/pulled-cubic-box-self-ipc.json`：ext+self merged IPC

## 验证

`tests/api/runSimConfig_parse_test.cpp` 覆盖了全部关键 parser 语义：默认值、参数分离、旧键拒绝、feasible line search 默认开启、friction 禁止、无效参数拒绝。

---

上一阶段：[IPC Overview](index.md)
下一阶段：[Phase 2](phase-2-external-barrier-energy.md)
