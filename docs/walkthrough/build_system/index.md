# 构建系统 Walkthrough

## 核心结论

`libpgo` 只有一套 C++ 后端（Conan + CMake），加两个入口壳：

```text
原生 C++ 入口: cmake --preset ...  ─┐
                                      ├→ BootstrapConan.cmake
Python 入口: uv sync / uv build    ─┘    → conanfile.py + recipes
                                          → CMakeToolchain + CMakeDeps
                                          → 根 CMakeLists.txt
                                              ├→ 原生库 / 工具 / tests
                                              └→ pypgo module / wheel
```

## Phase Map

| Phase | 主题 | 适合什么时候读 |
| --- | --- | --- |
| [Phase 01](phase01_stack_and_entrypoints.md) | 技术栈：Conan 与 scikit-build-core | 刚进入代码库，需要理解构建工具链 |
| [Phase 02](phase02_cmake_files.md) | 核心 CMake 文件：BootstrapConan、LibpgoConanDeps、根 CMakeLists | 想搞清楚 configure/build 阶段的完整控制流 |
| [Phase 03](phase03_recipes_and_ci.md) | Repo-local recipes 与 CI/CD | 想改依赖、添加 recipe、修改 CI 覆盖范围 |

## 核心文件地图

| 文件 | 角色 |
| --- | --- |
| `CMakePresets.json` | 原生 C++ 入口，统一 generator/binaryDir/build type/feature 组合 |
| `cmake/BootstrapConan.cmake` | Configure 入口：feature 规范化 → Conan install → 注入 toolchain |
| `cmake/LibpgoConanDeps.cmake` | `find_package()` 暴露 Conan 依赖 + 旧 target 名兼容别名 |
| `CMakeLists.txt` | Feature 开关、编译参数、target 工厂、安装逻辑 |
| `conanfile.py` | 依赖图声明、feature→option 映射、`CMakeToolchain`/`CMakeDeps` 生成 |
| `conan/recipes/` | 仓库自带 recipe，补齐 ConanCenter 上缺失或不匹配的包 |
| `pyproject.toml` | Python 构建入口，配置 `scikit-build-core` 参数 |
| `.github/workflows/ci-cd.yml` | CI/CD，复用本地同一套构建模型 |
