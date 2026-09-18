import importlib
import json
from pathlib import Path

import pytest


def _write_manifest(
    directory: Path,
    cases: list[dict],
    groups: dict[str, list[str]] | None = None,
) -> Path:
    manifest = directory / "cases.json"
    manifest.write_text(
        json.dumps(
            {
                "version": 1,
                "abc-output-root": "abc",
                "groups": groups or {},
                "cases": cases,
            }
        )
    )
    return manifest


def _case(name: str) -> dict:
    return {
        "name": name,
        "simulation-config": f"configs/{name}.json",
        "animation-config": f"configs/{name}-animation.json",
    }


def _write_case_configs(directory: Path, name: str) -> None:
    config_directory = directory / "configs"
    config_directory.mkdir(exist_ok=True)
    (config_directory / f"{name}.json").write_text("{}\n")
    (config_directory / f"{name}-animation.json").write_text("{}\n")


def test_batch_manifest_paths_are_relative_to_manifest(tmp_path: Path):
    cli = importlib.import_module("pypgo.pgo_run_cases")
    manifest = _write_manifest(tmp_path, [_case("first")])

    loaded_manifest = cli.load_manifest(manifest)

    assert loaded_manifest.cases == (
        cli.Case(
            name="first",
            simulation_config=(tmp_path / "configs/first.json").resolve(),
            animation_config=(
                tmp_path / "configs/first-animation.json"
            ).resolve(),
            abc_output=(tmp_path / "abc/first").resolve(),
        ),
    )
    assert loaded_manifest.groups == {}


def test_batch_cli_runs_only_requested_case(tmp_path: Path, monkeypatch):
    cli = importlib.import_module("pypgo.pgo_run_cases")
    for name in ("first", "second"):
        _write_case_configs(tmp_path, name)
    manifest = _write_manifest(tmp_path, [_case("first"), _case("second")])
    calls = []
    monkeypatch.setattr(
        cli.pypgo,
        "run_sim_from_config",
        lambda config: calls.append(("simulation", Path(config).name)) or 0,
    )
    monkeypatch.setattr(
        cli.pypgo,
        "convert_animation_to_abc",
        lambda config, output: calls.append(
            ("animation", Path(config).name, Path(output).name)
        )
        or 0,
    )

    assert cli.main([str(manifest), "second"]) == 0
    assert calls == [
        ("simulation", "second.json"),
        ("animation", "second-animation.json", "second"),
    ]
    assert (tmp_path / "abc/second").is_dir()
    assert not (tmp_path / "abc/first").exists()


def test_batch_cli_dry_run_prints_individual_pypgo_commands(
    tmp_path: Path, monkeypatch, capsys
):
    cli = importlib.import_module("pypgo.pgo_run_cases")
    _write_case_configs(tmp_path, "first")
    manifest = _write_manifest(tmp_path, [_case("first")])
    monkeypatch.setattr(
        cli.pypgo,
        "run_sim_from_config",
        lambda config: pytest.fail("dry-run executed a simulation"),
    )

    assert cli.main([str(manifest), "--dry-run"]) == 0
    output = capsys.readouterr().out
    assert "pgo-run-sim" in output
    assert "pgo-dump-abc" in output
    assert "mkdir -p" in output
    assert not (tmp_path / "abc").exists()


def test_batch_cli_expands_scene_groups_and_deduplicates_cases(
    tmp_path: Path, capsys
):
    cli = importlib.import_module("pypgo.pgo_run_cases")
    for name in ("first", "second", "third"):
        _write_case_configs(tmp_path, name)
    manifest = _write_manifest(
        tmp_path,
        [_case("first"), _case("second"), _case("third")],
        {"scene-a": ["first", "second"], "scene-b": ["second", "third"]},
    )

    assert cli.main(
        [
            str(manifest),
            "first",
            "--group",
            "scene-a",
            "--group",
            "scene-b",
            "--list",
        ]
    ) == 0
    assert capsys.readouterr().out.splitlines() == ["first", "second", "third"]


def test_batch_cli_lists_groups(tmp_path: Path, capsys):
    cli = importlib.import_module("pypgo.pgo_run_cases")
    manifest = _write_manifest(
        tmp_path,
        [_case("first"), _case("second")],
        {"scene": ["first", "second"]},
    )

    assert cli.main([str(manifest), "--list-groups"]) == 0
    assert capsys.readouterr().out == "scene: first, second\n"


def test_batch_cli_keep_going_reports_failure_and_runs_later_cases(
    tmp_path: Path, monkeypatch
):
    cli = importlib.import_module("pypgo.pgo_run_cases")
    for name in ("first", "second"):
        _write_case_configs(tmp_path, name)
    manifest = _write_manifest(tmp_path, [_case("first"), _case("second")])
    simulations = []
    animations = []

    def run_simulation(config: str) -> int:
        simulations.append(Path(config).name)
        return 4 if Path(config).name == "first.json" else 0

    monkeypatch.setattr(cli.pypgo, "run_sim_from_config", run_simulation)
    monkeypatch.setattr(
        cli.pypgo,
        "convert_animation_to_abc",
        lambda config, output: animations.append(Path(config).name) or 0,
    )

    assert cli.main([str(manifest), "--keep-going"]) == 1
    assert simulations == ["first.json", "second.json"]
    assert animations == ["second-animation.json"]


def test_batch_manifest_rejects_duplicate_names(tmp_path: Path):
    cli = importlib.import_module("pypgo.pgo_run_cases")
    manifest = _write_manifest(tmp_path, [_case("same"), _case("same")])

    with pytest.raises(cli.ManifestError, match="duplicate case name"):
        cli.load_manifest(manifest)


def test_batch_manifest_rejects_unknown_group_member(tmp_path: Path):
    cli = importlib.import_module("pypgo.pgo_run_cases")
    manifest = _write_manifest(
        tmp_path,
        [_case("first")],
        {"scene": ["missing"]},
    )

    with pytest.raises(cli.ManifestError, match="references unknown case"):
        cli.load_manifest(manifest)
