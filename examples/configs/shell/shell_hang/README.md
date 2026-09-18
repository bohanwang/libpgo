# Shell hang

This scene fixes one side of a shell and lets gravity deform the remaining
surface. It covers both solve modes supported by the shell IPC runner:

| Configuration | Mode | Timesteps | Purpose |
| --- | --- | ---: | --- |
| `shell-dynamic-ipc.json` | dynamic | 10000 | Time-dependent hanging-shell motion |
| `shell-static-ipc.json` | static | 1 | Direct static equilibrium under standard gravity |

Run either case directly or run the complete scene group:

```bash
uv run pgo-run-cases examples/cases.json --group shell
```

The simulation outputs are written to
`examples/generated/output/shell-dynamic-ipc/` and
`examples/generated/output/shell-static-ipc/`. Their sibling animation configs
can be exported through the animation stage of the same batch command.

Both configurations use standard gravity, `g = [0, -9.81, 0]`.
The reference static run converged to `solver-eps = 1e-4` after 1249 Newton
iterations; its configured upper bound is 5000 iterations.
