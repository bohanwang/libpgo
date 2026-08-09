These compact test fixtures are checked in so both build-tree and
installed-wheel tests can exercise the Python simulation entry point without
requiring mesher executables at test runtime. They are test-only assets, not
alternate versions of the retained examples.

The `installed-cli-*` configs also verify the packaged `pgo-run-sim`,
`pgo-dump-abc`, and `pgo-run-cases` console entry points.
