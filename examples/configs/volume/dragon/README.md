# Static dragon

This scene covers the static sampled-runner path for both volume mesh types:

| Configuration | Simulation mesh | Mode |
| --- | --- | --- |
| `dragon-static-tet-sampled.json` | tet | static |
| `dragon-static-cubic-sampled.json` | cubic | static |

Both configurations use the same dragon surface, gravity load, and fixed
geometric region. Because sampled contact is not assembled into the static
Newton energy, these cases exercise elasticity, gravity, and attachments only;
the runner prints an explicit warning.

Run both through their scene group:

```bash
uv run pgo-run-cases examples/cases.json --group dragon
```
