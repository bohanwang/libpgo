# Phase 03: 根 CMakeLists.txt 构建图

这一篇看的是“Conan bootstrap 之后，主 CMake 构建图是怎么展开的”。

入口文件只有一个：

- `CMakeLists.txt`

但它内部实际上分成了几层很明确的职责：

1. 项目初始化
2. feature 开关与级联推导
3. 编译参数和 interface target
4. Conan 依赖暴露与系统库探测
5. target 工厂函数
6. 子目录挂载
7. 安装 / `SKBUILD` 分支

## 1. 项目初始化

顶层一开始先做标准 CMake 工程初始化：

```cmake
cmake_minimum_required(VERSION 3.28.0 FATAL_ERROR)
project(libpgo LANGUAGES CXX C)
include(GNUInstallDirs)
```

然后设定全局语言和 ABI 基线：

```cmake
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
```

这里有两个值得注意的点：

1. `C++20` 是整个仓库的统一标准。
2. `CMAKE_POSITION_INDEPENDENT_CODE ON` 让静态库也按 PIC 编译，这对后面 Python module 链接内部静态库很重要。

## 2. feature 开关和级联关系

根 `CMakeLists.txt` 定义了一批 `option(...)`，例如：

- `PGO_PROFILE_FULL`
- `PGO_FEATURE_PYTHON`
- `PGO_FEATURE_ANIMATION_IO`
- `PGO_FEATURE_GEOMETRY_STACK`
- `PGO_FEATURE_LIBIGL`
- `PGO_FEATURE_GEOGRAM`
- `PGO_FEATURE_GMSH`
- `PGO_FEATURE_MKL`
- `PGO_FEATURE_ARPACK`
- `PGO_ENABLE_TESTS`

### 根 CMake 也会再做一次归一化

虽然 `BootstrapConan.cmake` 已经做过一次 feature 计算，根 `CMakeLists.txt` 仍然会重新计算一遍 `_LIBPGO_EFFECTIVE_*`，然后 `FORCE` 回写到 cache。

这样做的意义是：即便有人没有通过标准入口，而是直接以某种方式触发了根 configure，主构建图也仍然能保证 feature 关系自洽。

### 当前实现里的级联关系

1. `PGO_PROFILE_FULL=ON` 会强制打开：

- `PGO_FEATURE_PYTHON`
- `PGO_FEATURE_ANIMATION_IO`
- `PGO_FEATURE_GEOMETRY_STACK`
- `PGO_FEATURE_LIBIGL`
- `PGO_FEATURE_GEOGRAM`
- `PGO_FEATURE_GMSH`

2. 以下任一开启时，会强制打开 `PGO_FEATURE_GEOMETRY_STACK`：

- `PGO_FEATURE_PYTHON`
- `PGO_FEATURE_LIBIGL`
- `PGO_FEATURE_GEOGRAM`

3. `PGO_FEATURE_MKL` 和 `PGO_FEATURE_ARPACK` 不属于 full profile 的隐式展开内容；它们保持独立开关。

这个设计的目的不是“多做一次重复计算”，而是把 feature 逻辑集中写死在顶层，而不是散落到各个子目录里。

## 3. 编译参数如何被统一注入

根 `CMakeLists.txt` 没有把大堆编译参数直接散在每个 target 上，而是通过 interface target 聚合：

- `compilation_flag`
- `compilation_flag_for_debug`
- `cuda_compilation_flag`

### 为什么用 interface target

这样做有两个直接好处：

1. 目标函数 `add_libpgo_lib()` / `add_libpgo_tools()` 可以统一链接这些 flag target。
2. 平台分支和优化策略集中在顶层，而不是复制到几十个子 `CMakeLists.txt`。

### 这个阶段会处理哪些编译环境探测

#### OpenMP

```cmake
find_package(OpenMP)
```

如果成功：

- 记录 `OpenMP_CXX_FLAGS`
- 定义 `USE_OPENMP`
- 把 OpenMP 编译/链接参数注入 `compilation_flag`

#### AVX / ISA

- 非 MSVC：调用 `Find_AVX` 做自动探测
- MSVC：不自动猜测，而是通过缓存变量 `PGO_MSVC_ARCH` 显式选择 `default` / `AVX` / `AVX2` / `AVX512`

这说明当前仓库对 Windows 的策略是“显式配置优于自动猜测”。

#### Debug / Release 分支

当前实现中：

- `Debug` 路径会注入 `-O0`、`-ggdb3`
- POSIX 风格的 `compilation_flag_for_debug` 还包含 AddressSanitizer / UBSan / LeakSanitizer
- 非 Debug 路径会走 `-O3`

#### GNU/Clang 额外策略

顶层还会注入：

- `-Wall`
- `-Wextra`
- `-frounding-math`
- `-fvisibility=hidden`
- `-march=native`
- `-mtune=native`

对 GNU 还启用了：

- `-static-libstdc++`
- `-static-libgcc`

这些都是顶层“统一政策”，而不是单个模块的局部决定。

## 4. Conan 生成依赖如何暴露给主构建图

顶层先在必要时把 Conan generators 目录再一次放进搜索路径：

```cmake
if(DEFINED PGO_CONAN_GENERATORS_DIR AND EXISTS "${PGO_CONAN_GENERATORS_DIR}")
  list(PREPEND CMAKE_PREFIX_PATH "${PGO_CONAN_GENERATORS_DIR}")
  list(PREPEND CMAKE_MODULE_PATH "${PGO_CONAN_GENERATORS_DIR}")
endif()
```

然后包含：

```cmake
include(${CMAKE_SOURCE_DIR}/cmake/LibpgoConanDeps.cmake)
```

`LibpgoConanDeps.cmake` 负责两件事：

1. `find_package()` Conan 管理的包
2. 兼容旧 target 名称

### `find_package()` 的分层

#### Core dependencies

无条件尝试：

- `TBB` / `onetbb`
- `Eigen3`
- `fmt`
- `nlohmann_json`
- `spdlog`
- `argparse`
- `tinyobjloader`
- `stb`
- `autodiff`

#### feature-gated dependencies

按 feature 条件尝试：

- geometry stack: `Boost`、`CGAL`、`Ceres`、`NLopt`、`SuiteSparse`
- Python: `pybind11`
- libigl: `libigl`
- geogram: `geogram`
- animation I/O: `Alembic`、`Imath`
- gmsh: `gmsh`
- arpack: `arpackng`

#### repo-local packages

无条件尝试：

- `tetgen`
- `CCD_SafeCCD`
- `CCD_ExactCCD`
- `ASA`

#### tests

如果 `PGO_ENABLE_TESTS=ON`，尝试：

- `GTest`

### 兼容旧 target 名称的策略

这是 `LibpgoConanDeps.cmake` 最关键也最容易被忽略的一层。

当前源码树里已经存在很多历史 target 名称约定，而 Conan 包公开的 target 名不一定完全一致，所以这里做了一层别名桥接。

典型例子包括：

| 兼容目标 | 作用 |
| --- | --- |
| `TBB::tbb` <-> `onetbb::onetbb` | 兼容 oneTBB 命名差异 |
| `Boost::boost` <-> `boost::boost` | 兼容 Boost 命名差异 |
| `NLopt::nlopt` <-> `nlopt::nlopt` | 兼容大小写差异 |
| `Ceres::Ceres` <-> `Ceres::ceres` | 兼容历史 target 名 |
| `fmt::fmt-header-only` -> `fmt::fmt` | 对 header-only 用法做桥接 |
| `spdlog::spdlog_header_only` -> `spdlog::spdlog` | 对 header-only 用法做桥接 |
| `tiny_obj_loader` -> `tinyobjloader::*` | 保持旧源码里的 bare target 名 |
| `tetgen` -> `tetgen::tetgen` | 保持旧源码里的 bare target 名 |
| `ASA_lib` -> `ASA::ASA` | 保持旧源码里的 bare target 名 |
| `geogram` -> `geogram::geogram` | 保持旧源码里的 bare target 名 |

还有几个更偏修补性质的点：

- 把 `boost::boost` 归一回 headers-only 语义，避免不必要的链接库泄漏给下游
- 为 `CGAL` 补一个 `CGAL::TBB_support` shim
- 把 `CCD::SafeCCD` / `CCD::Exact` 映射到源码当前使用的 `CCD_SafeCCD` / `CCD_ExactCCD`

一句话总结：

> `LibpgoConanDeps.cmake` 的目标不是“优雅重命名”，而是尽量不改现有 `src/` 里大量 `target_link_libraries(...)` 的历史写法。

## 5. 系统库和非 Conan 依赖

当前顶层不是所有依赖都交给 Conan。

例如：

- `Threads` 由系统探测：
  `find_package(Threads REQUIRED)`
- `MKL` 通过 `cmake/third-party/mkl.cmake` 处理，而不是 ConanCenter
- `Knitro` 通过 `cmake/third-party/knitro.cmake` 处理

一旦对应 target 存在，还会注入：

- `PGO_HAS_MKL`
- `PGO_HAS_CGAL`
- `PGO_HAS_CERES`
- `PGO_HAS_KNITRO`

这些宏用于把“构建期探测结果”传递给源码编译条件。

## 6. target 工厂函数：为什么顶层封装了 `add_libpgo_lib()` 和 `add_libpgo_tools()`

这两个函数是根 `CMakeLists.txt` 里很重要的“内部 API”。

### `add_libpgo_lib()`

它统一处理：

- 创建静态库
- 输出目录
- `target_include_directories(... PUBLIC ./)`
- 统一链接 `compilation_flag`
- 如果开启 `PGO_RELEASE_MODE_DEBUG`，再链接 `compilation_flag_for_debug`
- 给 target 设置 `libraries` folder

### `add_libpgo_tools()`

它统一处理：

- 创建可执行文件
- 链接依赖 target
- 链接 `compilation_flag`
- 输出目录
- Windows runtime DLL 拷贝逻辑

这套封装让子目录不需要每次重新写：

- 平台分支
- 输出路径
- 通用编译参数
- Windows DLL 复制逻辑

## 7. 源码树如何挂到主构建图上

顶层实际挂载顺序很直接：

```cmake
add_subdirectory(src)
```

然后按条件决定：

- `PGO_ENABLE_TESTS=ON` 时：
  `add_subdirectory(tests)`
- `PGO_BUILD_SUBPROJECTS=ON` 时：
  `add_subdirectory(projects)`
- `PGO_FEATURE_PYTHON=ON` 时：
  `include(cmake/python_list.cmake)`

这意味着：

- `src/` 永远是主图的一部分
- tests、subprojects、Python 绑定都属于 feature-gated 或 option-gated 扩展层

## 8. `SKBUILD` 分支为什么跳过 `libpgoConfig.cmake`

顶层最后有一个很关键的分支：

```cmake
if(NOT SKBUILD)
  configure_package_config_file(...)
  install(FILES ... DESTINATION ...)
else()
  message(STATUS "SKBUILD detected, skipping generation of libpgoConfig.cmake")
endif()
```

这个分支的含义是：

- 原生 C++ 构建希望生成并安装 `libpgoConfig.cmake`
- `scikit-build-core` 驱动的 Python 构建不需要走这套“传统 CMake package install”路径

也就是说：

- `cmake --preset ...` 更偏“原生 C++ 项目 / 安装包”语境
- `uv sync` / `uv build` 更偏“我要产出 Python module”语境

两者共享同一个源码构建图，但发布终点不一样。

## 9. 输出布局

### 原生 C++ 构建

顶层把 `PGO_BINARY_ROOT` 设为当前二进制目录，也就是 preset 的 `binaryDir`。

对于 GNU/Clang/Apple Clang：

| 产物 | 目录 |
| --- | --- |
| 静态库 | `build/<preset>/lib` |
| 可执行文件 | `build/<preset>/bin` |

对于 MSVC：

| 产物 | 目录 |
| --- | --- |
| 静态库 | `build/<preset>/lib/Debug` 或 `build/<preset>/lib/Release` |
| 可执行文件 | `build/<preset>/bin/Debug` 或 `build/<preset>/bin/Release` |

### Python 构建

Python 路径由 `scikit-build-core` 把 build tree 固定在：

- `build/scikit-build`

此时：

- Conan 生成物在 `build/scikit-build/.conan-build/...`
- CMake 构建图也在这个 build tree 内
- `pypgo` module 通过 `src/api/python/pypgo/CMakeLists.txt` 安装到 wheel / editable 环境需要的位置

## 10. 这一层的最终结论

根 `CMakeLists.txt` 做的事可以压缩成一句话：

> 它把“feature、编译策略、依赖暴露、目标工厂、发布分支”集中在顶层统一定义，然后再把真正的源码实现下发到 `src/`、`tests/` 和 Python 绑定子树。

如果 Phase 02 解决的是“依赖环境从哪里来”，那么这一篇解决的是“依赖环境可用之后，主构建图怎样长出来”。

下一篇再看这套系统最外层的补充部分：

- [Phase 04: Recipes、Python 打包与 CI/CD](phase04_recipes_python_ci.md)
