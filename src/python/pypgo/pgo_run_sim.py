import argparse
from collections.abc import Sequence

import pypgo


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="pgo-run-sim",
        description="Run a libpgo simulation from a JSON configuration file.",
    )
    parser.add_argument("config", help="simulation JSON configuration file")
    arguments = parser.parse_args(argv)
    return pypgo.run_sim_from_config(arguments.config)


if __name__ == "__main__":
    raise SystemExit(main())
