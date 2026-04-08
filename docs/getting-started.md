# Getting Started

## Prerequisites

1. **Conan 2.x** — C++ dependency manager

    ```bash
    pip install conan
    conan profile detect
    ```

2. **uv** — Python build and development tool

    ```bash
    pip install uv
    ```

3. **Compilers**
    - GCC 11–13 (Ubuntu)
    - Apple Clang >= 15 (macOS)
    - Visual Studio 2022 (Windows)

## Quick Start (C++ Core)

```bash
cmake --preset core-release
cmake --build --preset core-release
ctest --preset core-release
```

## Quick Start (Python)

```bash
uv sync
uv run python -c "import pypgo; print(pypgo.__version__)"
```

See the [Build System Documentation](build-system.md) for the full toolchain
overview, preset-to-feature mapping, and Conan bootstrap flow.

See the [C++ Library Guide](guide/cpp-library.md) and [Python Bindings Guide](guide/python-bindings.md) for task-oriented build instructions.

## Building the Docs Locally

```bash
pip install -r docs/requirements.txt
mkdocs build --strict
mkdocs serve
```

The site will be available at `http://127.0.0.1:8000/`.
