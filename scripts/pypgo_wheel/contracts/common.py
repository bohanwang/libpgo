"""Release-wheel requirements shared by every supported platform."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final


RELEASE_DISTRIBUTION: Final = "pypgo"
PYTHON_TAG: Final = "cp312"
ABI_TAG: Final = "cp312"

COMMON_NATIVE_COMPONENTS: Final = ("tbb", "gmp", "gmpxx", "mpfr")


@dataclass(frozen=True)
class WheelPlatformContract:
    """Immutable release requirements for one wheel platform."""

    platform: str
    platform_tag: str
    repair_report: str
    required_native_components: tuple[str, ...]
    required_original_runtime_names: tuple[str, ...] = ()

    @property
    def expected_tag(self) -> str:
        return f"{PYTHON_TAG}-{ABI_TAG}-{self.platform_tag}"

    @property
    def required_audit_evidence(self) -> tuple[str, ...]:
        return (self.repair_report, "wheel-linkage.txt")
