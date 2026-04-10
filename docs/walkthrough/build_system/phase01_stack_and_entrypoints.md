# Phase 01: 技术栈

## Conan

### 是什么

Conan 是 C/C++ 的去中心化包管理器 — 可以类比为 C++ 的 uv / pip。它解决的核心问题是：给定一组依赖声明和目标平台的 ABI 约束（OS、compiler、C++ standard、build type），自动获取匹配的预编译二进制或在本地从源码构建。

### 不用 Conan 会怎样

**`find_package()` + 系统包管理器（apt/brew）：**

- 版本固定在发行版，不同机器上版本不一致
- ABI 不可控（libstdc++ 版本、C++ standard 不匹配会导致链接期 crash）
- `find_package()` 经常找不到库，需要手动设 `CMAKE_PREFIX_PATH`

**git submodule / FetchContent：**

- build 目录体积膨胀（每个项目重新编译所有依赖）
- 无法跨项目复用已编译的 package
- `add_subdirectory()` 机制容易污染主 CMake 的变量设置（反之亦然）

### 有什么好处

- **跨平台**：同一份 `conanfile.py` 在 Linux/macOS/Windows 上都能工作
- **构建系统无关**：Conan 本身不绑定 CMake，通过 generator 适配不同构建系统
- **ABI 兼容**：通过 profile 精确描述目标环境，自动匹配或构建 ABI 兼容的二进制
- **下载方便**：[ConanCenter](https://conan.io/center) 提供大量公共包；项目也可以自带 recipe 补齐缺失的包

### 核心概念

**Recipe**：描述一个包"怎么构建、怎么消费"的 Python 脚本（`conanfile.py`）。包含依赖声明、构建步骤、导出的 target 信息。ConanCenter 是公共 recipe 仓库，项目也可以自带 recipe（见 [Phase 03](phase03_recipes_and_ci.md)）。

**Profile**：描述目标平台的 ABI 契约。例如 `conan/profiles/linux-gcc-release`：

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

Conan 用 profile 做三件事：

1. `conan profile detect --force` — 自动检测当前编译器环境生成 profile
2. 根据 profile 匹配或构建 ABI 兼容的二进制包
3. 通过 generator 把 ABI 信息注入构建系统

**Package binary**：一个 recipe 在特定 profile 下的编译产物。Conan 用 package ID（recipe + settings + options 的 hash）做精确匹配。

**Local cache**：`~/.conan2/` 下的本地仓库，存放已下载/已编译的 package binary。多个项目共享同一份 cache，避免重复编译。

**Generator**：Conan 安装完依赖后，需要告诉构建系统"依赖在哪里、怎么链接"。本项目使用两个 CMake generator：

| Generator | 产物 | 作用 |
| --- | --- | --- |
| `CMakeToolchain` | `conan_toolchain.cmake` | 注入 ABI 信息：C++ standard、libstdc++ 版本、build type、compiler 路径、`CMAKE_PREFIX_PATH` |
| `CMakeDeps` | `*-config.cmake`（在 `build/<preset>/generators/`） | 供 `find_package()` 消费，内部硬编码绝对路径指向 `~/.conan2/p/` 中的头文件和库文件 |

### 如何工作

#### 1. 在 `conanfile.py` 中声明依赖

```python
class MyPkg(ConanFile):
    settings = "os", "compiler", "build_type", "arch"
    options  = {"with_feature": [True, False]}

    def build_requirements(self):     # 构建工具依赖
        self.tool_requires("cmake/[>=3.28]")

    def requirements(self):           # 运行时依赖
        self.requires("eigen/3.4.0")
        if self.options.with_feature:
            self.requires("boost/1.83.0")

    def generate(self):               # 配置 generator
        tc = CMakeToolchain(self)     # → conan_toolchain.cmake
        deps = CMakeDeps(self)        # → *-config.cmake
        tc.generate()
        deps.generate()

    def layout(self):
        cmake_layout(self)            # 设定目录布局
```

`requirements()` 声明"需要什么"，`generate()` 决定"怎么告诉 CMake"，两者分离。

#### 2. 安装依赖

```bash
conan install . --build=missing
```

这条命令的内部流程：

1. 解析 `conanfile.py` 的依赖树，处理不同库之间的版本冲突
2. 检查本地缓存中是否已有符合当前 profile 的二进制包
3. 如果没有，去远程仓库（ConanCenter）下载
4. 如果远程只有源码没有对应的二进制包（`--build=missing` 触发），在本地自动编译并缓存
5. 运行 generator，在 build 目录下生成 CMake 消费文件

#### 3. 告诉构建系统依赖在哪里

传统 `find_package()` 去 `/usr/local/lib/cmake` 找依赖 — 这依赖系统安装，不可靠。

Conan 的 `CMakeDeps` generator 在 `build/Release/generators/` 下生成 `*-config.cmake`，里面硬编码绝对路径指向 Conan 全局缓存（`~/.conan2/p/`）中的头文件和库文件。`CMakeToolchain` 同时把 generators 目录注入 `CMAKE_PREFIX_PATH`，这样 `find_package()` 就能找到。

#### 4. 告诉构建系统 ABI 信息

`CMakeToolchain` 生成的 `conan_toolchain.cmake` 注入：

- C++ standard（`CMAKE_CXX_STANDARD`）
- libstdc++ 版本选择
- Build type（Release/Debug）
- Compiler 路径
- `CMAKE_PREFIX_PATH`（指向 generators 目录）

#### 5. 使用生成的 toolchain 配置 CMake

```bash
cmake -B build -DCMAKE_TOOLCHAIN_FILE=build/conan_toolchain.cmake
cmake --build build
```

### 如何创建 Recipe

Recipe 有两种典型模式：

**编译本地源码：**

```python
class MyLib(ConanFile):
    exports_sources = "src/*", "CMakeLists.txt"

    def layout(self):
        cmake_layout(self)

    def generate(self):
        tc = CMakeToolchain(self)
        tc.generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        cmake = CMake(self)
        cmake.install()

    def package_info(self):
        self.cpp_info.libs = ["mylib"]      # 库文件名
        # 自动导出 target 名供下游 find_package 消费
```

**下载远程源码：**

```yaml
# conandata.yml
sources:
  "1.0":
    url: "https://github.com/foo/bar/archive/v1.0.tar.gz"
    sha256: "abc123..."
```

```python
class RemoteLib(ConanFile):
    def source(self):
        get(self, **self.conan_data["sources"][self.version])

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()
    # ...
```

**安装到 cache：**

- `conan export .` — 只把 recipe 写入 cache，不触发构建
- `conan install .` — 写入 recipe + 解析依赖 + 构建 + 打包

## scikit-build-core

### 是什么

scikit-build-core 是现代的 Python C++ 构建后端（PEP 517），把 Python 打包流程接回 CMake。

### 不用 scikit-build-core 会怎样

传统方式是 `setup.py` + `setuptools`：

- 硬编码大量 `extra_compile_args`（`-std=c++17`, `/O2`, `-fPIC`）
- `if sys.platform == "win32"` 到处散落
- 手写 `include_dirs` / `library_dirs`
- 改一处 C++ 代码要全量重编译（没有增量编译）

### 有什么好处

- **关注点分离**：C++ 编译交给 CMake，Python 打包配置放在 `pyproject.toml`
- **跨平台**：自动检测 Python 环境（路径、头文件、ABI 版本如 `cp310-cp310`），静默注入给 CMake
- **增量编译**：底层是 CMake + Ninja，改一个文件只重编译受影响的 target
- **静态配置**：所有配置在 `pyproject.toml` 中声明，不需要可执行的 `setup.py`

### 如何工作

当 `uv` 或 `pip` 启动构建时，scikit-build-core 接管并调用 CMake：

```text
uv sync / uv build / pip install
  → 读取 pyproject.toml
  → 发现 build-backend = scikit_build_core.build
  → scikit-build-core 接管：
      → cmake configure（传入 cmake.args）
      → cmake build
      → 把 CMake install 的产物打进 wheel
```

`pyproject.toml` 中的配置：

```toml
[build-system]
requires = ["scikit-build-core"]
build-backend = "scikit_build_core.build"

[tool.scikit-build]
build-dir = "build/scikit-build"          # 独立 build tree
cmake.args = [                            # 传给 cmake configure 的参数
    "-DCMAKE_TOOLCHAIN_FILE=../../cmake/BootstrapConan.cmake",
    "-DPGO_ENABLE_TESTS=OFF",
    "-DPGO_FEATURE_PYTHON=ON",
    "-DPGO_FEATURE_ANIMATION_IO=ON",
    "-DPGO_FEATURE_GEOMETRY_STACK=ON",
]
```

关键点：

- **同一个 CMake 后端**：Python 构建回到同一个 `BootstrapConan.cmake` + `conanfile.py` + 根 `CMakeLists.txt`，不是另一套编译系统
- **独立 build tree**：`build-dir` 与 preset 的 build 目录分开，互不干扰
- **SKBUILD 变量**：CMake 端通过 `if(SKBUILD)` 区分 Python 打包 vs 原生 C++ 安装，两条路共享源码构建图但安装目标不同

### 与 uv 的关系

`uv` 是 Python 环境和包管理工具（类似 pip + venv 的组合）。`uv sync` 读取 `pyproject.toml`，发现 build-backend 是 scikit-build-core，就调用它执行 CMake 构建并安装 Python module。`uv build` 做同样的事但产出 wheel 分发包。

## 三条构建入口

### 1. 原生 C++

```bash
cmake --preset core-release
cmake --build --preset core-release
```

`CMakePresets.json` → `CMAKE_TOOLCHAIN_FILE=BootstrapConan.cmake` → Conan install → 根 `CMakeLists.txt` → Ninja 编译。

常用 preset：

| Preset | 内容 |
| --- | --- |
| `core-release` | 核心构建 |
| `geometry-release` | 开启几何栈 |
| `python-release` | Python + animation I/O + geometry stack |
| `full-release` | `PGO_PROFILE_FULL=ON`，完整特性 |

Preset 可以通过 CMake 的 feature 开关（`PGO_FEATURE_*`）灵活组合构建内容。

### 2. Python 开发

```bash
uv sync
```

`pyproject.toml` → `scikit-build-core` → CMake → `BootstrapConan.cmake` → 同一套后端 → `pypgo` module。

### 3. Python 发包

```bash
uv build
```

与 `uv sync` 同一个 build backend，目标从本地开发安装变成分发产物。

---

下一阶段：[Phase 02: 核心 CMake 文件](phase02_cmake_files.md)
