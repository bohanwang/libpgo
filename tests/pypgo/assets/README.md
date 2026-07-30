These lightweight fixtures were generated from
`examples/assets/volume/box/box.obj` with:

```text
python examples/scripts/generate_lite_tet_cubic_asset.py \
  --build-dir build/pypgo \
  --cubic-resolution 2 \
  --output-dir <temporary-directory>
```

They are checked in so both build-tree and installed-wheel tests can exercise
the Python simulation entry point without requiring mesher executables at test
runtime.
