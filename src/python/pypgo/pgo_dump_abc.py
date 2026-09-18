import argparse
from collections.abc import Sequence

import pypgo


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="pgo-dump-abc",
        description="Convert libpgo OBJ animation output to Alembic files.",
    )
    parser.add_argument("config", help="animation JSON configuration file")
    parser.add_argument("output_folder", help="folder for generated Alembic files")
    arguments = parser.parse_args(argv)
    return pypgo.convert_animation_to_abc(
        arguments.config,
        arguments.output_folder,
    )


if __name__ == "__main__":
    raise SystemExit(main())
