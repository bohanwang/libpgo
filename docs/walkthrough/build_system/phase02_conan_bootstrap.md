# Phase 02: Conan Bootstrap 控制流

这一篇只看一个文件：`cmake/BootstrapConan.cmake`。

它是 `libpgo` 构建系统的真正 configure 入口，因为无论你从哪条路径进入：

- `cmake --preset ...`
- `uv sync`
- `uv build`

最后都会先经过这个文件。

## 它到底负责什么

可以把 `BootstrapConan.cmake` 理解成 configure 阶段的依赖总调度器。它负责：

1. 规范化 feature 开关
2. 规范化 profile 和 Conan 输出目录
3. 判断当前配置是否需要重新跑 Conan
4. 在需要时导出本地 recipe 并执行 `conan install`
5. 把 Conan 生成的 toolchain 和 package config 接回 CMake

它不负责编译源码，也不直接声明 `src/` 里的 target；那是根 `CMakeLists.txt` 的工作。

## 输入和输出

### 主要输入

`BootstrapConan.cmake` 消费的核心输入包括：

- `CMAKE_BUILD_TYPE`
- `PGO_CONAN_CPPSTD`
- `PGO_CONAN_PROFILE_HOST`
- `PGO_CONAN_PROFILE_BUILD`
- `PGO_PROFILE_FULL`
- `PGO_FEATURE_PYTHON`
- `PGO_FEATURE_ANIMATION_IO`
- `PGO_FEATURE_GEOMETRY_STACK`
- `PGO_FEATURE_LIBIGL`
- `PGO_FEATURE_GEOGRAM`
- `PGO_FEATURE_GMSH`
- `PGO_FEATURE_MKL`
- `PGO_FEATURE_ARPACK`

这些值有些来自 preset，有些来自 `pyproject.toml`，有些来自命令行 `-D...`。

### 主要输出

它最终会准备出三类结果：

1. Conan 输出根目录

- 原生 preset 常见位置：
  `build/<preset>/.conan-build/<BuildType>/`
- Python 路径常见位置：
  `build/scikit-build/.conan-build/<BuildType>/`

2. Conan 生成器目录

- `.../build/<BuildType>/generators/`

3. CMake 侧的可消费对象

- `conan_toolchain.cmake`
- 各个依赖包的 `*-config.cmake`
- runtime activation script，例如 `conanrun.sh` / `conanrun.bat`

## 控制流总览

下面是按当前实现整理过的控制流：

```text
1. 规范化 CMAKE_BUILD_TYPE
2. 把关键变量加入 CMAKE_TRY_COMPILE_PLATFORM_VARIABLES
3. 计算 effective features
4. 解析 source dir、cppstd、host/build profile
5. 计算 Conan 输出目录与 generators 目录
6. 生成 configuration signature
7. 检查 signature 和 toolchain 是否已存在
8. 若需要：
   8.1 找到 conan 命令
   8.2 找到 host Python
   8.3 运行 conan/recipes/export_recipes.py
   8.4 执行 conan install
   8.5 写回新的 feature-signature.txt
9. include(conan_toolchain.cmake)
10. prepend generators 目录到 CMAKE_PREFIX_PATH / CMAKE_MODULE_PATH
```

## 1. 默认 `CMAKE_BUILD_TYPE`

文件一开始就做了一个兜底：

```cmake
if(NOT DEFINED CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
```

这意味着：

- preset 没给时，默认是 `Release`
- Python 路径也会继承这个默认值，除非别处显式改掉

这一步看起来很普通，但它影响后面几乎所有路径计算：

- Conan 输出根目录
- signature 内容
- `conan install -s build_type=...`

## 2. 为什么要写 `CMAKE_TRY_COMPILE_PLATFORM_VARIABLES`

这个文件显式把一批变量加入 `CMAKE_TRY_COMPILE_PLATFORM_VARIABLES`。

原因是：CMake 在 configure 阶段会做一些 `try_compile`，而这些子 configure 也需要看到当前构建环境的核心变量。否则会出现：

- 主 configure 知道自己用了哪个 build type / 哪组 feature
- `try_compile` 子 configure 却看不到这些值

`BootstrapConan.cmake` 现在传进去的包括：

- `CMAKE_BUILD_TYPE`
- `PGO_CONAN_OUTPUT_DIR`
- `PGO_CONAN_GENERATORS_DIR`
- `PGO_CONAN_CPPSTD`
- 所有主要 `PGO_FEATURE_*`

这一步的目标不是“多传一点变量”，而是让 configure 子上下文和主上下文保持一致。

## 3. feature 归一化与级联推导

`_libpgo_compute_effective_features()` 做了第一轮 feature 规范化。

核心规则是：

1. `PGO_PROFILE_FULL=ON` 会强制打开：

- `PGO_FEATURE_PYTHON`
- `PGO_FEATURE_ANIMATION_IO`
- `PGO_FEATURE_GEOMETRY_STACK`
- `PGO_FEATURE_LIBIGL`
- `PGO_FEATURE_GEOGRAM`
- `PGO_FEATURE_GMSH`

2. 以下任一开启，都会隐式要求几何栈：

- `PGO_FEATURE_PYTHON`
- `PGO_FEATURE_LIBIGL`
- `PGO_FEATURE_GEOGRAM`

最后这些结果会被 `FORCE` 回写到 CMake cache。

这很重要，因为后面 Conan options 不能根据“半展开状态”的 feature 去生成。它必须拿到一套已经补齐依赖关系的最终值。

## 4. host/build profile 和路径规范化

这个阶段主要做三件事：

### 4.1 计算 source dir

```cmake
get_filename_component(_LIBPGO_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)
```

后面执行：

- `export_recipes.py`
- `conan install "${_LIBPGO_SOURCE_DIR}"`

都依赖这个绝对路径。

### 4.2 规范化 `PGO_CONAN_CPPSTD`

如果外部没传，默认使用：

```cmake
set(PGO_CONAN_CPPSTD "20" CACHE STRING ... FORCE)
```

这和根 `CMakeLists.txt` 的 `CMAKE_CXX_STANDARD 20` 对齐。

### 4.3 把 profile 解析成绝对路径

如果 `PGO_CONAN_PROFILE_HOST` / `PGO_CONAN_PROFILE_BUILD` 是相对路径，bootstrap 会把它们解释为“相对仓库根目录”的路径，再转成绝对路径。

这样做的好处是：

- 本地命令行、preset、CI 都可以稳定地传 profile
- 后续写 signature 时不会因为相对路径写法不同而产生伪差异

## 5. Conan 输出目录与 generators 目录

如果外部没有手动指定，bootstrap 会自动推导默认目录：

- Conan 输出根目录：
  `CMAKE_BINARY_DIR/.conan-build/${CMAKE_BUILD_TYPE}`
- generators 目录：
  `${conan_root}/build/${CMAKE_BUILD_TYPE}/generators`

以 `core-release` 为例，常见路径是：

```text
build/core-release/.conan-build/Release/
build/core-release/.conan-build/Release/build/Release/generators/
```

这套布局的好处是：

- Conan 生成物跟当前 build tree 放在一起
- 不同 preset、不同 build type 彼此隔离
- 删除某个 build tree 时，不会把别的入口的生成物一起删掉

## 6. feature signature：为什么它能避免重复跑 Conan

这个文件会把当前构建环境压缩成一个字符串签名，写到：

- `feature-signature.txt`

签名内容包括：

- `build_type`
- `cppstd`
- host/build profile
- `profile_full`
- 所有主要 feature 对应的 `with_*`

一个典型签名长这样：

```text
build_type=Release
cppstd=20
profile_host=/abs/path/to/profile
profile_build=/abs/path/to/profile
profile_full=False
with_python=True
with_animation_io=True
...
```

然后 bootstrap 会做两件检查：

1. `conan_toolchain.cmake` 是否存在
2. `feature-signature.txt` 是否存在且与当前完全一致

只有在两者都满足时，才会跳过 `conan install`。

这意味着重复 configure 的常见快路径是：

- feature 没变
- profile 没变
- build type 没变
- toolchain 还在

如果上述任何一个变化，bootstrap 就会重新跑 Conan。

## 7. 找 host Python：为什么不是直接写死 `python3`

bootstrap 里有一个 `_libpgo_find_host_python()`：

1. 如果环境变量 `UV_PYTHON` 存在且指向可执行文件，优先使用它
2. 否则再 `find_program(python3 python REQUIRED)`

这一步的直接用途是运行：

- `conan/recipes/export_recipes.py`

优先认 `UV_PYTHON` 的意义在于：Python 打包路径里，当前活跃 Python 环境往往就是由 `uv` 管理的。

## 8. 为什么要先跑 `export_recipes.py`

在真正执行 `conan install` 之前，bootstrap 会先运行：

```bash
python conan/recipes/export_recipes.py
```

这是因为 repo 里有一批本地 recipe：

- `asa`
- `ccd-safe`
- `ccd-exact`
- `autodiff`
- `tetgen`
- `suitesparse`
- `arpack-ng`
- `geogram`
- `gmsh`
- `alembic`

这些 recipe 先 `conan export` 到本地 Conan cache，后面的 `conan install` 才能把它们当成正常依赖解析。

这里要区分两件事：

- `conan export`：只把 recipe 写进 cache
- `conan install`：解析依赖、拉包、必要时编译、生成 CMake 配置、准备运行环境

所以 `export_recipes.py` 不是可选优化，而是 repo 自带 recipe 生效的前提。

## 9. `conan install` 实际传了什么

bootstrap 构造出来的命令核心形态是：

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

如果 host/build profile 被设置了，还会附带：

```bash
-pr:h <host_profile>
-pr:b <build_profile>
```

其中最关键的几件事是：

1. `--build missing`

- 本地 cache 和远端都没有匹配二进制时，允许 Conan 本地构建缺失包。

2. `-o &:...`

- 这些 option 传给当前根 recipe，也就是 `conanfile.py`。
- `&:` 这里表示“当前消费中的根包”。

3. profile + `compiler.cppstd`

- 共同决定当前 ABI / toolchain 约束。

## 10. `CMakeToolchain` 和 `CMakeDeps` 分别解决什么问题

这两个东西都由 `conanfile.py` 的 `generate()` 生成，但职责不同。

### `CMakeToolchain`

它的主要作用是把“Conan 解析出来的构建环境”注入给 CMake，包括：

- 编译器相关设置
- ABI 相关设置
- 构建类型
- Conan 侧推导出的 CMake 变量

在 `libpgo` 里，`generate()` 把以下变量写回 toolchain：

- `PGO_PROFILE_FULL`
- `PGO_FEATURE_PYTHON`
- `PGO_FEATURE_ANIMATION_IO`
- `PGO_FEATURE_GEOMETRY_STACK`
- `PGO_FEATURE_LIBIGL`
- `PGO_FEATURE_GEOGRAM`
- `PGO_FEATURE_GMSH`
- `PGO_FEATURE_MKL`
- `PGO_FEATURE_ARPACK`

这让 CMake configure 后续看到的是一套已经和 Conan 对齐过的 feature 状态。

### `CMakeDeps`

它的作用是生成给 `find_package()` 使用的 package config。

也就是说，根 `CMakeLists.txt` 和 `LibpgoConanDeps.cmake` 后续调用：

```cmake
find_package(TBB CONFIG)
find_package(Eigen3 CONFIG)
find_package(pybind11 CONFIG)
```

能找到的那套 `*-config.cmake`，很多就是 `CMakeDeps` 生成在 generators 目录里的。

一句话区分：

- `CMakeToolchain` 解决“当前怎么编”
- `CMakeDeps` 解决“依赖在哪里，`find_package()` 怎么找”

## 11. 为什么要把 generators 目录 prepend 到 `CMAKE_PREFIX_PATH`

bootstrap 结尾会：

```cmake
list(PREPEND CMAKE_PREFIX_PATH "${_LIBPGO_GENERATORS_DIR}")
list(PREPEND CMAKE_MODULE_PATH "${_LIBPGO_GENERATORS_DIR}")
```

这一步是整个 Conan -> CMake handoff 最关键的桥接之一。

为什么是 prepend，而不是 append？

因为 `libpgo` 希望当前 configure 尽量优先消费“这次 Conan install 刚生成的那套配置”，而不是：

- 系统目录里的旧版本
- 用户机器上其他项目留下的全局安装

这能显著减少“`find_package()` 找到了包，但不是当前 Conan 解析出来的那份”的漂移问题。

## 12. bootstrap 结束后，控制权交给谁

bootstrap 最后只做两件事：

1. `include("${_LIBPGO_TOOLCHAIN_FILE}")`
2. 把 generators 目录塞进 CMake 搜索路径

之后它就结束了，控制权回到根 `CMakeLists.txt`。

也就是说，这个文件的职责边界非常清楚：

- 它负责让依赖环境和 feature 环境可用
- 它不负责定义源码 target

源码 target、编译 flag、输出目录、安装逻辑，都在下一层完成。

下一篇继续看这个下一层：

- [Phase 03: 根 CMakeLists.txt 构建图](phase03_root_cmake_graph.md)
