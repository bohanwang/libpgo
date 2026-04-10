# Phase 02: 核心 CMake 文件

本阶段讲解三个核心 CMake 文件的完整控制流：

1. `cmake/BootstrapConan.cmake` — configure 入口，负责依赖环境就绪
2. `cmake/LibpgoConanDeps.cmake` — 导入依赖库 + 兼容性别名
3. `CMakeLists.txt` — 主构建图

## 1. `BootstrapConan.cmake`

所有入口（preset / uv sync / uv build）共享的 configure 入口。它负责让依赖环境和 feature 环境可用，不负责定义源码 target。

### 完整控制流

```text
[开始解析 BootstrapConan.cmake]
        |
        v
[1. 初始化与变量设置]
  * 设置默认 CMAKE_BUILD_TYPE = Release（如果未定义）
  * 定义需要传递给 Conan 的 CMake 变量列表
        |
        v
[2. 特性推导与规范化 (_libpgo_compute_effective_features)]
  * 读取传入的特性开关（如 PGO_FEATURE_PYTHON）
  * 执行依赖推导：
      - PGO_PROFILE_FULL=ON → 强制开启所有子特性
      - Python / libigl / geogram 任一开启 → 强制开启 GEOMETRY_STACK
  * 将 CMake 的 ON/OFF 转换为 Conan 的 True/False 字符串
  * 结果 FORCE 回写 CMake cache
        |
        v
[3. 环境与路径探测]
  * 获取 source dir、cppstd（默认 20）
  * 解析 host/build profile 的绝对路径
  * 设定 Conan 输出目录和 generators 目录
        |
        v
[4. 构建指纹签名 (Feature Signature)]
  * 将 build_type、cppstd、profile、所有 feature 开关
    拼接成签名字符串
        |
        v
[5. 指纹比对]
  * 检查磁盘上 feature-signature.txt 和 conan_toolchain.cmake
        |
   +---------+---------+
   |                   |
  匹配               不匹配/不存在
   |                   |
   v                   v
[跳过 Conan]    [6. 查找 Python 环境]
   |              * _libpgo_find_host_python()
   |              * 优先 UV_PYTHON，fallback python3
   |                   |
   |                   v
   |            [7. 导出本地 Recipe]
   |              * 运行 export_recipes.py
   |              * 把仓库自带 recipe 写入 Conan cache
   |                   |
   |                   v
   |            [8. 执行 conan install]
   |              * 拼接命令，传入所有 feature 开关和 profile
   |              * Conan 下载/编译依赖
   |              * 在 generators 目录产生 CMake 配置文件
   |                   |
   |                   v
   |            [9. 更新签名]
   |              * 写入新的 feature-signature.txt
   |                   |
   +---------+---------+
             |
             v
[10. 注入 Conan 环境到 CMake]
  * include(conan_toolchain.cmake)
  * Prepend generators 目录到 CMAKE_PREFIX_PATH / CMAKE_MODULE_PATH
  * 控制权交还给根 CMakeLists.txt
```

### `conan install` 的命令形态

```bash
conan install <source_dir> \
  --output-folder <conan_root> \
  --build missing \
  -s build_type=<CMAKE_BUILD_TYPE> \
  -s compiler.cppstd=<PGO_CONAN_CPPSTD> \
  -o &:with_python=<True|False> \
  -o &:with_animation_io=<True|False> \
  -o &:with_geometry_stack=<True|False> \
  ...
```

- `--build missing`：cache 中缺失时允许本地构建
- `-o &:...`：option 传给根 recipe（`conanfile.py`）
- 如果设置了 profile，附带 `-pr:h` / `-pr:b`

### Feature Signature 避免重复跑 Conan

跳过 `conan install` 的条件：`conan_toolchain.cmake` 存在 **且** `feature-signature.txt` 与当前完全一致。任何一项变化都会触发重新安装。

### 输出目录布局

以 `core-release` 为例：

```text
build/core-release/.conan-build/Release/
build/core-release/.conan-build/Release/build/Release/generators/
```

不同 preset / build type 彼此隔离。

### Conan 生成的两类产物

| 产物 | 由谁生成 | 解决什么 |
| --- | --- | --- |
| `conan_toolchain.cmake` | `CMakeToolchain` | "当前怎么编" — 编译器/ABI 设置、feature 变量回写 |
| `*-config.cmake` | `CMakeDeps` | "依赖在哪里" — `find_package()` 消费的 package config |

### 结尾桥接

```cmake
include("${_LIBPGO_TOOLCHAIN_FILE}")
list(PREPEND CMAKE_PREFIX_PATH "${_LIBPGO_GENERATORS_DIR}")
list(PREPEND CMAKE_MODULE_PATH "${_LIBPGO_GENERATORS_DIR}")
```

Prepend 而非 append — 优先消费当前 Conan install 生成的配置，避免系统目录漂移。

## 2. `LibpgoConanDeps.cmake`

在 `BootstrapConan.cmake` 注入 Conan 环境之后，由根 `CMakeLists.txt` 加载。职责：导入所有依赖库 + 建立兼容性别名。

### 完整控制流

```text
[开始解析 LibpgoConanDeps.cmake]
        |
        v
[1. 防重复包含守卫]
  * 检查 LIBPGO_CONAN_DEPS_INCLUDED 是否已定义
  * 若已定义 → return
        |
        v
[2. 定义辅助工具函数]
  * _libpgo_try_find_package()：包装 find_package，
    静默查找 CONFIG 模式的库，打印成功/失败日志
  * _libpgo_add_ns_alias()：创建 INTERFACE "幽灵目标"
    作为别名，链接到真实 Conan 目标
        |
        v
[3. 加载核心必选依赖]
  * TBB/onetbb、Eigen3、fmt、nlohmann_json、spdlog、
    argparse、tinyobjloader、stb、autodiff
        |
        v
[4. 加载条件可选依赖（由 feature 开关控制）]
  * GEOMETRY_STACK → Boost、CGAL、Ceres、NLopt、SuiteSparse
  * PYTHON → pybind11
  * LIBIGL → libigl
  * GEOGRAM → geogram
  * ANIMATION_IO → Alembic、Imath
  * GMSH → gmsh
  * ARPACK → arpackng
        |
        v
[5. 加载内部私有库与测试库]
  * tetgen、CCD_SafeCCD、CCD_ExactCCD、ASA
  * PGO_ENABLE_TESTS=ON → GTest
        |
        v
[6. 兼容性别名映射]
  * 目的：让 src/ 目录下旧代码的 target_link_libraries
    名字不用改
        |
        v
[脚本结束]
```

### 别名映射详表

| 映射类型 | 具体操作 | 原因 |
| --- | --- | --- |
| 双向绑定 | `TBB::tbb` ↔ `onetbb::onetbb` | oneTBB 命名差异 |
| 大小写桥接 | `Boost::boost` ↔ `boost::boost`、NLopt/nlopt、Ceres/ceres | 不同 CMake 版本 target 名大小写不一致 |
| Header-only 桥接 | `fmt::fmt-header-only` → `fmt::fmt`、`spdlog::spdlog_header_only` → `spdlog::spdlog` | Conan 导出的 target 名 vs 源码中 link 的名字 |
| 去命名空间 | `tinyobjloader::*` → `tiny_obj_loader`、`tetgen::tetgen` → `tetgen`、`geogram::geogram` → `geogram`、`ASA::ASA` → `ASA_lib` | 源码中使用 bare target 名 |
| 补齐特性 | 创建 `CGAL::TBB_support` 幽灵目标 | CGAL 内部 TBB shim |
| 降级处理 | `boost::boost` 强制降级为只依赖 `Boost::headers` | 防止 Boost.Test 等不需要的库被强行拉进下游（如 CGAL） |
| 私有库重命名 | `CCD::SafeCCD` → `CCD_SafeCCD` 等 | 历史代码中的 target 名 |

### 系统库和非 Conan 依赖

- `Threads`：系统探测
- MKL：`cmake/third-party/mkl.cmake`
- Knitro：`cmake/third-party/knitro.cmake`

探测成功后注入宏（`PGO_HAS_MKL`、`PGO_HAS_CGAL` 等）供源码条件编译。

## 3. 根 `CMakeLists.txt`

Conan bootstrap 之后，根 `CMakeLists.txt` 展开主构建图。

### 完整控制流

```text
[开始解析 CMakeLists.txt]
        |
        v
[1. 全局基础初始化]
  * cmake_minimum_required(VERSION 3.28.0)
  * project(libpgo LANGUAGES CXX C)
  * CMAKE_CXX_STANDARD = 20，PIC = ON
        |
        v
[2. Feature 开关定义与级联推导]
  * 定义 PGO_FEATURE_* 系列 option
  * 依赖推导（与 BootstrapConan 同一套规则）：
      - PROFILE_FULL=ON → 强制开启所有子特性
      - Python/libigl/geogram → 强制开启 GEOMETRY_STACK
  * 推导结果 FORCE 回写 CMake cache
        |
        v
[3. 配置编译参数（Interface Target 模式）]
  * 探测 OpenMP 和 AVX/AVX2/AVX512
  * 创建 interface target：
      - compilation_flag：平台编译参数
      - compilation_flag_for_debug：Debug 专用参数
  * 平台分支：
      - GNU/Clang：-O3, -march=native, -Wall -Wextra,
        -fvisibility=hidden, -frounding-math,
        -static-libstdc++ -static-libgcc, Sanitizers
      - MSVC：/O2, /bigobj, /MP (多核编译),
        AVX 通过 PGO_MSVC_ARCH 显式选择
  * 所有 flag 统一绑定到 interface target
        |
        v
[4. 接入第三方依赖库]
  * include(LibpgoConanDeps.cmake) → 上一节的流程
  * 探测 MKL、CGAL、Ceres、Knitro
  * 找到的库注入 C++ 宏（-DPGO_HAS_MKL 等）
        |
        v
[5. 定义 Target 工厂函数]
  * add_libpgo_lib()：
      - 创建静态库
      - 统一输出目录（lib/Debug 或 lib/Release）
      - 自动继承 compilation_flag
      - PUBLIC include 设置
  * add_libpgo_tools()：
      - 创建可执行文件
      - 统一输出目录（bin/）
      - Windows：powershell 脚本自动拷贝依赖 DLL
        到 exe 旁边
        |
        v
[6. 下发构建任务（挂载子目录）]
  * add_subdirectory(src)       — 永远挂载
  * add_subdirectory(tests)     — PGO_ENABLE_TESTS=ON
  * add_subdirectory(projects)  — PGO_BUILD_SUBPROJECTS=ON
  * include(cmake/python_list.cmake) — PGO_FEATURE_PYTHON=ON
        |
        v
[7. 发布与安装策略]
  * 检查 SKBUILD 变量：
      |
   +--+--+
   |     |
  Yes    No（传统 C++ 开发者）
   |     |
   |     v
   |  [生成 libpgoConfig.cmake]
   |  [install() 到 /usr/local/lib 等]
   |     |
   v     v
[结束：生成底层构建图（Ninja / MSBuild）]
```

### Interface Target 编译参数详表

| Target | 内容 |
| --- | --- |
| `compilation_flag` | `-Wall -Wextra -frounding-math -fvisibility=hidden -march=native`、OpenMP、AVX |
| `compilation_flag_for_debug` | `-O0 -ggdb3`、ASan/UBSan/LeakSan |
| `cuda_compilation_flag` | CUDA 编译参数 |

子目录不需要处理平台分支、输出路径、通用编译参数 — 全部由 interface target 和工厂函数封装。

### 输出布局

**GNU/Clang/Apple Clang：**

| 产物 | 目录 |
| --- | --- |
| 静态库 | `build/<preset>/lib` |
| 可执行文件 | `build/<preset>/bin` |

**MSVC：**

| 产物 | 目录 |
| --- | --- |
| 静态库 | `build/<preset>/lib/{Debug,Release}` |
| 可执行文件 | `build/<preset>/bin/{Debug,Release}` |

**Python：** build tree 在 `build/scikit-build`，Conan 生成物在其 `.conan-build/` 子目录。

---

上一阶段：[Phase 01](phase01_stack_and_entrypoints.md)
下一阶段：[Phase 03: Recipes 与 CI/CD](phase03_recipes_and_ci.md)
