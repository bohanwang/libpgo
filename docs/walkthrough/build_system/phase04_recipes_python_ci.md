# Phase 04: Recipes、Python 打包与 CI/CD

前两篇已经解释了：

- `BootstrapConan.cmake` 怎样把依赖环境准备好
- 根 `CMakeLists.txt` 怎样生成主构建图

这一篇看三块外围但很关键的内容：

1. repo-local Conan recipes
2. `conanfile.py` 和 `pyproject.toml`
3. GitHub Actions 里的 CI/CD 是如何复用这套模型的

## 1. repo-local Conan recipes 的角色

`libpgo` 并不完全依赖 ConanCenter 上现成的包形态。仓库里有一批自带 recipe，位于：

- `conan/recipes/asa/all/conanfile.py`
- `conan/recipes/ccd-safe/all/conanfile.py`
- `conan/recipes/ccd-exact/all/conanfile.py`
- `conan/recipes/autodiff/all/conanfile.py`
- `conan/recipes/tetgen/all/conanfile.py`
- `conan/recipes/suitesparse/all/conanfile.py`
- `conan/recipes/arpack-ng/all/conanfile.py`
- `conan/recipes/geogram/all/conanfile.py`
- `conan/recipes/gmsh/all/conanfile.py`
- `conan/recipes/alembic/all/conanfile.py`

这些 recipe 的存在说明两件事：

1. 这个仓库希望把依赖解析控制权尽量放在 repo 内，而不是完全依赖外部包仓库的默认形态。
2. 某些包即使在 ConanCenter 上存在，repo 也可能需要自己的 recipe 来固定版本、补充 patch、修正 package layout 或统一 target 暴露方式。

## 2. `export_recipes.py` 具体做了什么

`conan/recipes/export_recipes.py` 的职责很单纯：

- 维护一张 `{name, version}` 列表
- 对每个 recipe 运行：
  `conan export <recipe_dir> --name <name> --version <version>`

它只负责把 recipe 写进本地 Conan cache，并不直接触发依赖解析和构建。

这就是为什么整个流程要分成两步：

1. `conan export`
2. `conan install`

### 什么时候只需要 `conan export`

如果你只是修改了某个本地 recipe，自然需要先把新 recipe 导出到 cache。

### 什么时候会进入完整链路

一旦触发 `BootstrapConan.cmake` 里的 `conan install`，Conan 才会继续做：

- 依赖图解析
- 下载二进制或源码
- 必要时本地构建缺失包
- 生成 `CMakeToolchain`
- 生成 `CMakeDeps`

## 3. `conanfile.py` 是这套依赖系统的声明中心

`conanfile.py` 当前实现里最重要的部分有五个。

### 3.1 `build_requirements()`

```python
def build_requirements(self):
    self.tool_requires("cmake/[>=3.28 <4.0.0]")
    self.tool_requires("ninja/[>=1.11 <2.0]")
```

这表明：

- Conan 不只管运行时 / 链接时依赖
- 它也可以把构建工具链的一部分纳入依赖模型

### 3.2 `requirements()`

这里定义了真正的 C++ 依赖图。

#### 核心无条件依赖

- `onetbb`
- `eigen`
- `fmt`
- `nlohmann_json`
- `spdlog`
- `argparse`
- `tinyobjloader`
- `stb`
- `autodiff`
- `tetgen`
- `ccd-safe`
- `ccd-exact`
- `asa`

#### 条件依赖

由 `_effective_features()` 控制：

- geometry stack: `boost`、`cgal`、`ceres-solver`、`nlopt`、`suitesparse`
- Python: `pybind11`
- animation I/O: `alembic`、`imath`
- libigl: `libigl`
- geogram: `geogram`
- gmsh: `gmsh`
- arpack: `arpack-ng`

此外还固定有：

```python
self.test_requires("gtest/1.14.0")
```

也就是说，测试依赖的声明也集中在这里，而不是散落在 CMake 中。

### 3.3 `_effective_features()`

这是 Conan 世界里的 feature 归一化逻辑，对应 CMake 世界里的 effective feature 推导。

它会做三类事：

1. 展开 `profile_full`
2. 让 `with_python` / `with_libigl` / `with_geogram` 隐式拉起 `with_geometry_stack`
3. 产出一组标准化后的布尔字典，供 `requirements()` 和 `generate()` 共用

这样 Conan 依赖图和 CMake feature 图不会各算各的。

### 3.4 `layout()`

```python
def layout(self):
    cmake_layout(self)
```

这个调用很重要，因为它决定了 Conan 为 CMake 生成文件时采用的目录布局。

配合 `BootstrapConan.cmake` 里的 `--output-folder`，最终你会看到类似这样的结构：

```text
build/core-release/.conan-build/Release/
build/core-release/.conan-build/Release/build/Release/generators/
```

### 3.5 `generate()`

这里是 Conan 到 CMake 的桥梁。

```python
tc = CMakeToolchain(self)
deps = CMakeDeps(self)
```

它做了两件关键事：

1. 生成 `CMakeToolchain`

- 把 feature 状态回写成 CMake 变量
- 让后续 CMake configure 和 Conan 解析结果对齐

2. 生成 `CMakeDeps`

- 产生后续 `find_package()` 消费的 package config 文件

当前实现还做了一个细节：

```python
tc.user_presets_path = False
```

这表示 Conan 不去额外生成用户 preset 文件，仓库自己维护 `CMakePresets.json` 作为主入口。

## 4. `pyproject.toml` 如何把 Python 构建接回 C++

`pyproject.toml` 在这个仓库里做的是“Python 包的入口描述”，不是“手写一套替代 CMake 的扩展编译逻辑”。

### build backend

```toml
[build-system]
build-backend = "scikit_build_core.build"
```

这一步决定：

- `uv sync`
- `uv build`

最终都会交给 `scikit-build-core`。

### build 目录

```toml
[tool.scikit-build]
build-dir = "build/scikit-build"
```

这样 Python 构建有自己稳定的 build tree，不会和 preset 目录混在一起。

### 传给 CMake 的参数

当前默认参数是：

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

这几行可以直接读成一句话：

> Python 构建不是另一套后端，它只是用 `scikit-build-core` 帮你把 CMake 调起来，然后继续走 `BootstrapConan.cmake`。

### `src/api/python/pypgo/CMakeLists.txt` 在做什么

这个子目录负责定义 Python module target：

- `add_library(pypgo MODULE ...)`
- 链接 `pybind11::module`
- 链接内部静态库 `pgo_c_static`
- `install(TARGETS pypgo LIBRARY DESTINATION . RUNTIME DESTINATION .)`

所以 Python wheel / editable install 里的二进制核心，实际上仍然是这套主 CMake 构建图的一部分。

## 5. 当前 CI/CD 的真实覆盖范围

这里必须以 `.github/workflows/ci-cd.yml` 的当前实现为准，而不是凭印象总结。

### 5.1 C++ job

当前 C++ matrix 是：

| OS | preset |
| --- | --- |
| `ubuntu-latest` | `core-release` |
| `ubuntu-latest` | `geometry-release` |
| `ubuntu-latest` | `full-release` |
| `macos-latest` | `core-release` |
| `windows-latest` | `core-release` |

这说明当前公开 workflow 的事实是：

- Linux 的 C++ matrix 只覆盖 `core-release`、`geometry-release`、`full-release`
- macOS 和 Windows 当前只覆盖 `core-release`
- 并没有在这个 workflow 里把所有 preset 全量跑一遍

这也是本系列需要显式修正的一点：不要把文档写成“Ubuntu 构建了所有 profile”。

### C++ job 的控制流

每个平台都遵循同一套形态：

1. 安装 CMake
2. 安装 Python
3. 安装 `conan==2.21.0`
4. `conan profile detect --force`
5. 选择仓库里的 profile 文件并写入：
   - `PGO_CONAN_PROFILE_HOST`
   - `PGO_CONAN_PROFILE_BUILD`
6. 安装平台额外系统依赖
7. `cmake --preset <preset> -DPGO_CONAN_PROFILE_HOST=... -DPGO_CONAN_PROFILE_BUILD=...`
8. `cmake --build --preset <preset>`
9. `ctest --preset <preset>`

Windows 额外处理了：

- Git long paths
- `msvc-dev-cmd`
- `conanrun.bat` 激活后再跑 `ctest`

### 5.2 Python job

当前 Python job 只跑在 `ubuntu-latest`，流程是：

1. 安装 CMake
2. 安装 `uv`
3. 安装 Python 3.11
4. 安装 Conan
5. 选择 `conan/profiles/linux-gcc-release`
6. 通过环境变量 `CMAKE_ARGS` 把 host/build profile 传给 `uv`
7. 执行：

```bash
uv sync --python "${UV_PYTHON_VERSION}"
```

8. 找到 Conan 生成的 `conanrun.sh`
9. `source` 之后执行：

```bash
uv run pytest
```

这条链验证的是：

- Python build backend
- 同一套 Conan bootstrap
- Python tests

### 5.3 CD release job

tag 发布时，`cd-release` job 会在 Ubuntu 上执行：

```bash
uv build
```

同样，它也会先安装 Conan、选定 repository profile，并通过 `CMAKE_ARGS` 把 profile 传进 Python build backend。

因此 release job 并没有引入第二套“发布专用构建系统”；它只是换了一个前门。

## 6. 当前系统的最终理解方式

把整个构建系统压缩成最短的一张关系图，可以写成：

```text
原生入口: cmake --preset ...
Python 入口: uv sync / uv build
        |
        v
BootstrapConan.cmake
        |
        v
conanfile.py + conan/recipes/
        |
        v
CMakeToolchain + CMakeDeps
        |
        v
根 CMakeLists.txt + src/tests/api/python
        |
        +--> 原生库 / 工具 / tests
        |
        +--> pypgo module / wheel
```

也就是说，`libpgo` 当前的工程化重点不是“支持很多入口”，而是“无论从哪个入口进来，最终都收敛到同一个构建后端”。

回到系列首页：

- [构建系统 Walkthrough](index.md)
