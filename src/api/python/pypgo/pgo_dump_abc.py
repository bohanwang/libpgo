import argparse
from pathlib import Path

import pypgo


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert an anim.json sequence config into an Alembic file."
    )
    parser.add_argument("anim_json", help="Path to the anim.json config file.")
    parser.add_argument(
        "output_dir",
        nargs="?",
        help="Optional output directory. Defaults to the directory containing anim.json.",
    )
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    output_dir = args.output_dir
    if output_dir is None:
        output_dir = str(Path(args.anim_json).parent or Path("."))

    print(pypgo.convert_animation_to_abc(args.anim_json, output_dir))
