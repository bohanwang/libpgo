import argparse
import json
import sys
import zipfile
from pathlib import Path

import pytest


SCRIPTS_DIR = Path(__file__).resolve().parents[2] / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

import release_wheel_provenance as provenance  # noqa: E402
from release_wheel_contract import get_platform_contract  # noqa: E402


def write_wheel(
    wheel_dir: Path,
    platform: str,
) -> Path:
    contract = get_platform_contract(platform)
    wheel = wheel_dir / f"pypgo-0.0.4-{contract.expected_tag}.whl"
    dist_info = "pypgo-0.0.4.dist-info"
    with zipfile.ZipFile(wheel, "w") as archive:
        archive.writestr(
            f"{dist_info}/METADATA",
            "Metadata-Version: 2.1\nName: pypgo\nVersion: 0.0.4\n",
        )
        archive.writestr(
            f"{dist_info}/WHEEL",
            "Wheel-Version: 1.0\n"
            f"Tag: {contract.expected_tag}\n",
        )
    return wheel


def write_file(directory: Path, name: str) -> Path:
    path = directory / name
    path.write_text(f"evidence for {name}\n", encoding="utf-8")
    return path


def record_args(
    tmp_path: Path,
    platform: str,
) -> tuple[argparse.Namespace, dict[str, str]]:
    source_dir = tmp_path / "source"
    evidence_dir = tmp_path / "evidence"
    wheel_dir = tmp_path / "wheelhouse"
    source_dir.mkdir()
    evidence_dir.mkdir()
    wheel_dir.mkdir()

    checkout = {"commit": "a" * 40, "tree": "b" * 40}
    (evidence_dir / provenance.SOURCE_RECORD).write_text(
        json.dumps({"source": checkout}),
        encoding="utf-8",
    )
    write_wheel(wheel_dir, platform)
    contract = get_platform_contract(platform)
    audit_evidence = [
        write_file(evidence_dir, name)
        for name in reversed(contract.required_audit_evidence)
    ]

    args = argparse.Namespace(
        source_dir=source_dir,
        evidence_dir=evidence_dir,
        wheel_dir=wheel_dir,
        platform=platform,
        cmake_cache=write_file(evidence_dir, "CMakeCache.txt"),
        dependency_evidence=[
            write_file(evidence_dir, "uv-freeze.txt"),
            write_file(evidence_dir, "system-packages.txt"),
        ],
        audit_evidence=audit_evidence,
        test_evidence=[
            write_file(evidence_dir, "installed-pytest.txt"),
            write_file(evidence_dir, "installed-smoke.txt"),
        ],
    )
    return args, checkout


@pytest.mark.parametrize("platform", ["macos", "linux", "windows"])
def test_record_captures_platform_and_sorted_evidence(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    platform: str,
) -> None:
    args, checkout = record_args(tmp_path, platform)
    monkeypatch.setattr(provenance, "read_project_version", lambda _source: "0.0.4")
    monkeypatch.setattr(provenance, "require_clean_checkout", lambda _source: checkout)

    provenance.command_record(args)

    record = json.loads(
        (args.evidence_dir / provenance.FINAL_RECORD).read_text(encoding="utf-8")
    )
    assert record["schema_version"] == 1
    assert record["platform"] == platform
    assert record["wheel"]["tags"] == [
        get_platform_contract(platform).expected_tag
    ]
    assert [item["name"] for item in record["audit_evidence"]] == sorted(
        get_platform_contract(platform).required_audit_evidence
    )
    assert [item["name"] for item in record["dependency_evidence"]] == [
        "system-packages.txt",
        "uv-freeze.txt",
    ]
    assert [item["name"] for item in record["test_evidence"]] == [
        "installed-pytest.txt",
        "installed-smoke.txt",
    ]


def test_record_rejects_missing_evidence_file(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    args, checkout = record_args(tmp_path, "linux")
    args.audit_evidence[0].unlink()
    monkeypatch.setattr(provenance, "read_project_version", lambda _source: "0.0.4")
    monkeypatch.setattr(provenance, "require_clean_checkout", lambda _source: checkout)

    with pytest.raises(provenance.ProvenanceError, match="missing evidence files"):
        provenance.command_record(args)
