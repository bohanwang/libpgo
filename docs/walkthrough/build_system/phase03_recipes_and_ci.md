# Phase 03: Repo-Local Recipes 与 CI/CD

## 1. Repo-Local Conan Recipes

`conan/recipes/` 下的自带 recipe：asa、ccd-safe、ccd-exact、autodiff、tetgen、suitesparse、arpack-ng、geogram、gmsh、alembic。

这些包在 ConanCenter 上缺失或不完全匹配 repo 需要（版本、patch、target 暴露方式），因此 repo 自带 recipe 来固定控制。

`export_recipes.py` 对每个 recipe 运行 `conan export`，只写 recipe 到 cache，不触发依赖解析和构建。

## 2. `conanfile.py` — 依赖声明中心

### `build_requirements()`

```python
self.tool_requires("cmake/[>=3.28 <4.0.0]")
self.tool_requires("ninja/[>=1.11 <2.0]")
```

构建工具链也纳入依赖模型。

### `requirements()`

**Core 无条件依赖：** onetbb、eigen、fmt、nlohmann_json、spdlog、argparse、tinyobjloader、stb、autodiff、tetgen、ccd-safe、ccd-exact、asa

**条件依赖（由 `_effective_features()` 控制）：**

- Geometry stack: boost、cgal、ceres-solver、nlopt、suitesparse
- Python: pybind11
- Animation I/O: alembic、imath
- libigl / geogram / gmsh / arpack: 各自独立

**测试：** `gtest/1.14.0`（`test_requires`）

### `_effective_features()`

Conan 侧的 feature 归一化，与 CMake 侧规则一致（展开 profile_full、隐式拉起 geometry_stack）。产出标准化布尔字典供 `requirements()` 和 `generate()` 共用。

### `generate()`

```python
tc = CMakeToolchain(self)    # → conan_toolchain.cmake（feature 变量回写）
deps = CMakeDeps(self)        # → *-config.cmake（find_package 消费）
tc.user_presets_path = False  # 仓库自己维护 CMakePresets.json
```

### `layout()`

`cmake_layout(self)` — 决定 Conan 生成文件的目录布局，配合 bootstrap 的 `--output-folder`。

## 3. `pyproject.toml` — Python 构建入口

```toml
[build-system]
build-backend = "scikit_build_core.build"

[tool.scikit-build]
build-dir = "build/scikit-build"
cmake.args = [
    "-DCMAKE_TOOLCHAIN_FILE=../../cmake/BootstrapConan.cmake",
    "-DPGO_ENABLE_TESTS=OFF",
    "-DPGO_FEATURE_PYTHON=ON",
    "-DPGO_FEATURE_ANIMATION_IO=ON",
    "-DPGO_FEATURE_GEOMETRY_STACK=ON",
]
```

Python 构建有独立 build tree（不与 preset 目录混），但回到同一个 C++ 后端。

`src/api/python/pypgo/CMakeLists.txt` 定义 Python module target（`pybind11::module` + 内部静态库 `pgo_c_static`），install 到 wheel / editable 环境。

## 4. CI/CD

以 `.github/workflows/ci-cd.yml` 当前实现为准。

### C++ Job Matrix

在三个平台的 GitHub runner 上进行构建测试，并在 Ubuntu 上构建了所有 profile：

| OS | Preset |
| --- | --- |
| `ubuntu-latest` | `core-release`, `geometry-release`, `full-release` |
| `macos-latest` | `core-release` |
| `windows-latest` | `core-release` |

每个平台的流程：

1. 安装 CMake + Python + `conan==2.21.0`
2. `conan profile detect --force`
3. 选择仓库 profile → `PGO_CONAN_PROFILE_HOST` / `PGO_CONAN_PROFILE_BUILD`
4. 安装平台系统依赖
5. `cmake --preset <preset> -DPGO_CONAN_PROFILE_HOST=... -DPGO_CONAN_PROFILE_BUILD=...`
6. `cmake --build --preset <preset>`
7. `ctest --preset <preset>`

Windows 额外：Git long paths、`msvc-dev-cmd`、`conanrun.bat` 激活后再 ctest。

### Python Job

只跑 `ubuntu-latest`：

1. 安装 CMake + uv + Python 3.11 + Conan
2. 选择 `conan/profiles/linux-gcc-release`
3. 通过 `CMAKE_ARGS` 环境变量传 profile
4. `uv sync`
5. `source conanrun.sh && uv run pytest`

### CD Release Job

Tag 发布时在 Ubuntu 执行 `uv build`，同一套 Conan bootstrap，没有第二套发布专用构建系统。

---

上一阶段：[Phase 02: 核心 CMake 文件](phase02_cmake_files.md)
返回总览：[构建系统 Walkthrough](index.md)
