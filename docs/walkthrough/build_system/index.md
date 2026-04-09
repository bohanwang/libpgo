# 构建系统 Walkthrough

这组文档面向已经开始阅读 `libpgo` 源码、需要理解“构建系统到底是怎么串起来的”贡献者和维护者，而不是只想快速跑通一次安装的用户。

如果你只需要快速上手，请先看：

- [Getting Started](../../getting-started.md)
- [Build System Documentation](../../build-system.md)

如果你想知道下面这些问题，这组 walkthrough 才是主入口：

- 为什么这个仓库同时用了 `Conan`、`CMake`、`CMakePresets`、`scikit-build-core` 和 `uv`
- `cmake --preset ...`、`uv sync`、`uv build` 最终分别会落到哪些文件
- `BootstrapConan.cmake` 是怎么决定要不要重新跑 `conan install`
- Conan 生成的 `toolchain`、`generators`、`find_package()` 配置是怎么被根 `CMakeLists.txt` 接进来的
- Python 打包为什么没有绕开 C++ 构建链，而是回到了同一套 Conan + CMake 后端

## 先读什么

| 文档 | 主题 | 适合什么时候读 |
| --- | --- | --- |
| [Phase 01: 技术栈与构建入口](phase01_stack_and_entrypoints.md) | `Conan`、`scikit-build-core`、`CMakePresets`、`uv` 的分工 | 刚进入代码库，先建立整体心智模型 |
| [Phase 02: Conan Bootstrap 控制流](phase02_conan_bootstrap.md) | `BootstrapConan.cmake` 如何规范化 feature、导出 recipe、运行 `conan install`、注入 toolchain | 想搞清楚 configure 阶段到底做了什么 |
| [Phase 03: 根 CMakeLists.txt 构建图](phase03_root_cmake_graph.md) | 根 `CMakeLists.txt` 如何定义 feature、编译参数、依赖暴露、target 工厂和安装分支 | 想读懂主构建图，或准备改顶层 CMake |
| [Phase 04: Recipes、Python 打包与 CI/CD](phase04_recipes_python_ci.md) | 本地 Conan recipe、`conanfile.py`、`pyproject.toml`、GitHub Actions | 想改依赖、Python wheel、CI 覆盖范围 |

## 三个常见入口的落点

| 用户入口 | 第一跳 | 中间层 | 最终作用点 |
| --- | --- | --- | --- |
| `cmake --preset core-release` | `CMakePresets.json` | `cmake/BootstrapConan.cmake` | Conan 解析依赖，根 `CMakeLists.txt` 生成原生构建图 |
| `uv sync` | `pyproject.toml` | `scikit-build-core` + `cmake/BootstrapConan.cmake` | 构建 `pypgo` 可编辑开发环境 |
| `uv build` | `pyproject.toml` | `scikit-build-core` + `cmake/BootstrapConan.cmake` | 构建 Python 分发产物 |

一个重要结论是：`libpgo` 其实只有一套 C++ 后端，只是有多个前门。

- 原生 C++ 开发从 `CMakePresets.json` 进入。
- Python 开发和发包从 `pyproject.toml` 进入。
- 但两条路都会回到同一个 `cmake/BootstrapConan.cmake`，再进入同一个根 `CMakeLists.txt`。

## 核心文件地图

| 文件 | 角色 |
| --- | --- |
| `CMakePresets.json` | 为原生 C++ 开发提供稳定入口，统一 `generator`、`binaryDir`、`build type` 和常见 feature 组合 |
| `cmake/BootstrapConan.cmake` | configure 入口；计算有效 feature、规范化 profile、决定是否重跑 Conan，并注入生成结果 |
| `cmake/LibpgoConanDeps.cmake` | 调用 `find_package()` 暴露 Conan 依赖，并给旧 target 名做兼容性 alias |
| `CMakeLists.txt` | 定义 feature 开关、编译参数、target 工厂、输出目录、测试挂载和安装逻辑 |
| `conanfile.py` | 声明依赖图、feature 到 Conan option 的映射，并生成 `CMakeToolchain` / `CMakeDeps` |
| `conan/recipes/` | 仓库自带 recipe，补齐 ConanCenter 上没有或不完全匹配 repo 需要的包 |
| `pyproject.toml` | Python 构建入口，配置 `scikit-build-core`、`uv` 和 `pypgo` 构建参数 |
| `.github/workflows/ci-cd.yml` | 把本地同一套构建模型搬到 GitHub Actions 上执行 |

## 这组文档采用什么边界

- 以当前 repo 已提交实现为准，不照搬旧设计草图中的过时结论。
- 重点解释控制流、文件职责、feature 传播和输出布局。
- 不重复 `Getting Started` 里的安装步骤，只在需要时引用入口命令。
- 不引入额外的文档渲染插件，流程图全部用纯 Markdown 表格、步骤和代码块表达。

## 建议阅读顺序

1. [Phase 01: 技术栈与构建入口](phase01_stack_and_entrypoints.md)
2. [Phase 02: Conan Bootstrap 控制流](phase02_conan_bootstrap.md)
3. [Phase 03: 根 CMakeLists.txt 构建图](phase03_root_cmake_graph.md)
4. [Phase 04: Recipes、Python 打包与 CI/CD](phase04_recipes_python_ci.md)
