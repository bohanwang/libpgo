# Phase 01: 技术栈与构建入口

这一篇先回答两个最基础的问题：

1. 为什么 `libpgo` 的构建系统不是“只有一个 CMakeLists.txt”那么简单。
2. 为什么这个仓库把 `Conan`、`CMakePresets`、`scikit-build-core`、`uv` 叠在了一起。

## 整体分层

先看最粗粒度的构建分层：

| 层 | 入口文件 / 工具 | 责任 |
| --- | --- | --- |
| 用户入口 | `cmake --preset ...`、`uv sync`、`uv build` | 选择“原生 C++”还是“Python 打包”入口 |
| 入口配置 | `CMakePresets.json`、`pyproject.toml` | 提供常见 feature 组合和 Python 构建参数 |
| 依赖 bootstrap | `cmake/BootstrapConan.cmake` | 规范化 feature，运行 `conan install`，接入 Conan 生成结果 |
| 依赖图定义 | `conanfile.py`、`conan/recipes/` | 声明包依赖、本地 recipe、`CMakeToolchain`、`CMakeDeps` |
| 主构建图 | `CMakeLists.txt`、`src/`、`tests/` | 编译静态库、工具、C API、Python module |
| Python 构建后端 | `scikit-build-core` | 把 Python 打包流程接回 CMake |

所以真正的心智模型不是“有一个大 CMake 项目”，而是：

- CMake 负责主构建图。
- Conan 负责 C++ 依赖解析与 ABI 对齐。
- `CMakePresets.json` 负责原生构建入口。
- `scikit-build-core` 负责 Python 打包入口。
- `uv` 负责 Python 环境与打包命令。

## Conan 在这个仓库里负责什么

### 它是什么

`Conan` 是 C++ 包管理器。放在这个仓库的语境里，可以把它近似理解为：

- `pip` / `uv` 在 Python 世界里负责“解析依赖、下载包、处理环境”
- `Conan` 在 C++ 世界里负责“解析依赖、选择 ABI、生成给构建系统消费的依赖元数据”

但这个类比只能帮助入门，不能把它当成完全等价关系。`Conan` 额外要管一件 Python 包管理器通常不用显式处理的事：ABI。

### 为什么这里不用“纯系统包管理器 + find_package”

如果只依赖系统包管理器，常见问题是：

- 版本在不同平台上不一致
- `find_package()` 找到的是系统目录下的任意一个安装版本，结果不稳定
- 开发者本机、CI、发布环境之间容易出现“头文件找到了，但链接的是另一套二进制”的问题
- 对 Windows/MSVC、Linux/GCC、macOS/Apple Clang 这几种 ABI 组合，不容易做统一管理

`libpgo` 现在的实现选择把这些问题前置到 Conan：

- 依赖版本在 `conanfile.py` 里统一声明
- 平台/编译器/`cppstd`/runtime 通过 Conan profile 约束
- `find_package()` 最终尽量从 Conan 生成目录里找包，而不是优先去系统目录碰运气

### 为什么这里不用 submodule / FetchContent 当主方案

它们当然能工作，但这里不是主线，原因主要有三类：

1. 依赖复用差

- `FetchContent` / `add_subdirectory()` 更适合把第三方源码直接编进当前工程。
- 同一个依赖很难在多个工程之间共享一份缓存好的二进制结果。

2. 主工程和子工程的 CMake 容易互相污染

- `add_subdirectory()` 把第三方项目直接并进当前 CMake 作用域。
- 第三方 CMake 的全局变量、编译选项、policy、option 名称可能和主工程互相干扰。

3. 构建目录和配置复杂度容易失控

- 源码级依赖意味着每个新 build tree 都可能重新配置和编译一遍第三方。
- 对大型依赖组合，构建目录和 configure/build 时间会明显放大。

这也是为什么 `libpgo` 现在把第三方依赖的“构建”和“消费”分离开：

- 依赖解析交给 Conan
- 主构建图仍然保留在 CMake

### Conan 在这个仓库里的实际工作方式

这里不是“在 `CMakeLists.txt` 里直接调用一堆 `find_package` 看系统装了什么”，而是：

1. 在 `conanfile.py` 里声明依赖和 feature 对应关系。
2. 在 `BootstrapConan.cmake` 里把当前 CMake feature 规范化。
3. 调用 `conan install`。
4. 由 Conan 生成：
   - `conan_toolchain.cmake`
   - 一组 `*-config.cmake`
   - runtime activation script
5. 根 `CMakeLists.txt` 再通过 `LibpgoConanDeps.cmake` 把这些包找出来并做 target 别名兼容。

也就是说，`find_package()` 仍然在用，但它找的不是 `/usr/local/lib/cmake` 下的随机安装，而是 Conan 刚刚生成、并指向 Conan cache 的那一套配置文件。

### profile 为什么重要

在这个仓库里，profile 不是可有可无的附属品，而是 ABI 契约。

一个 Conan profile 至少会描述：

- 操作系统
- 架构
- 编译器
- 编译器版本
- `compiler.cppstd`
- C++ 标准库 ABI，例如 `libstdc++11`
- `build_type`

例如仓库里的 Linux profile `conan/profiles/linux-gcc-release` 明确写了：

```ini
[settings]
os=Linux
arch=x86_64
compiler=gcc
compiler.version=13
compiler.cppstd=20
compiler.libcxx=libstdc++11
build_type=Release
```

这意味着 Conan 不只是“帮你下载包”，还会用 profile 决定：

- 当前环境应该匹配哪个二进制包
- 拉取不到预编译包时，本地要按什么 ABI 重新编译

本地开发通常至少要先运行：

```bash
conan profile detect --force
```

而 CI 不靠自动探测，而是显式选择仓库里的 profile 文件。

## scikit-build-core 在这个仓库里负责什么

### 它是什么

`scikit-build-core` 是现代 Python 构建后端，用来把 Python 打包流程和 CMake 对接起来。

在 `libpgo` 里，它的角色不是替代 CMake，而是让 Python 构建不要绕开 CMake。

### 为什么这里不用传统 `setup.py` 扩展编译逻辑

如果 Python 扩展直接靠 `setup.py` 手写编译参数，通常会慢慢演变成：

- 大量平台分支
- 大量硬编码 `extra_compile_args`
- 手写 `include_dirs` / `library_dirs`
- 头文件、Python ABI、编译器、依赖库路径之间耦合越来越重

而 `libpgo` 已经有一套非平凡的 C++ 构建系统：

- 顶层 feature 开关
- Conan 依赖 bootstrap
- 原生工具和静态库
- Python module 依赖同一批内部 target

因此更合理的做法是：

- Python 层只描述“我要构建一个包”
- 真正的 C++ 构建继续交给 CMake

### 这里它是怎么接到 CMake 的

`pyproject.toml` 里把 `scikit-build-core` 声明为 build backend：

```toml
[build-system]
build-backend = "scikit_build_core.build"
```

同时在 `tool.scikit-build` 下传入了核心 CMake 参数：

```toml
cmake.args = [
    "-DCMAKE_TOOLCHAIN_FILE=../../cmake/BootstrapConan.cmake",
    "-DPGO_ENABLE_TESTS=OFF",
    "-DPGO_FEATURE_PYTHON=ON",
    "-DPGO_FEATURE_ANIMATION_IO=ON",
    "-DPGO_FEATURE_GEOMETRY_STACK=ON",
    "-DPYPGO_VERSION_INFO=0.0.3",
]
```

这几行非常关键，因为它说明 Python 构建并没有走“另一套依赖管理逻辑”，而是又回到了：

- `BootstrapConan.cmake`
- `conanfile.py`
- 根 `CMakeLists.txt`

## 三条用户入口的真实落点

### 1. 原生 C++ 开发

最典型的入口是：

```bash
cmake --preset core-release
cmake --build --preset core-release
```

控制流是：

1. `CMakePresets.json` 选择 `generator=Ninja`、`binaryDir`、`CMAKE_BUILD_TYPE`
2. preset 把 `CMAKE_TOOLCHAIN_FILE` 指到 `cmake/BootstrapConan.cmake`
3. `BootstrapConan.cmake` 运行 Conan，并注入生成结果
4. 根 `CMakeLists.txt` 生成主构建图
5. `Ninja` 负责编译

### 2. Python 开发环境

典型入口是：

```bash
uv sync
```

控制流是：

1. `uv` 解析 Python 环境
2. `pyproject.toml` 指定 `scikit-build-core`
3. `scikit-build-core` 调 CMake
4. CMake 仍然从 `BootstrapConan.cmake` 开始
5. 生成并编译 `pypgo` module

### 3. Python 发包

典型入口是：

```bash
uv build
```

这条链和 `uv sync` 使用同一个 Python build backend，只是目标从“本地开发安装”变成了“构建分发产物”。

## preset 在这个仓库里是什么角色

`CMakePresets.json` 不是第二套构建系统，而是原生 C++ 入口的参数模板。

例如：

- `core-release` 代表核心原生构建
- `geometry-release` 开启几何栈
- `python-release` 开启 Python + animation I/O + geometry stack
- `full-release` 通过 `PGO_PROFILE_FULL=ON` 打开完整特性组合

这让开发者不必每次手写一串 `-D...=ON`：

```bash
cmake --preset full-release
```

而不是：

```bash
cmake -S . -B build/full-release \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_TOOLCHAIN_FILE=cmake/BootstrapConan.cmake \
  -DPGO_PROFILE_FULL=ON
```

## 这一层的最终心智模型

把这一层压缩成一句话：

> `libpgo` 不是“CMake + 一点点 Python”，而是“一套 Conan + CMake 的 C++ 后端，再加两个入口壳：CMake presets 和 Python packaging”。

后面三篇就是沿着这条主线继续往下拆：

- [Phase 02: Conan Bootstrap 控制流](phase02_conan_bootstrap.md)
- [Phase 03: 根 CMakeLists.txt 构建图](phase03_root_cmake_graph.md)
- [Phase 04: Recipes、Python 打包与 CI/CD](phase04_recipes_python_ci.md)
