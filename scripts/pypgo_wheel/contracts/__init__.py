"""Registry for selecting a release-wheel contract by platform name."""

from __future__ import annotations

from types import MappingProxyType
from typing import Final, Mapping

from .common import WheelPlatformContract
from .linux import LINUX_CONTRACT
from .macos import MACOS_CONTRACT
from .windows import WINDOWS_CONTRACT


PLATFORM_CONTRACTS: Final[Mapping[str, WheelPlatformContract]] = MappingProxyType(
    {
        contract.platform: contract
        for contract in (LINUX_CONTRACT, MACOS_CONTRACT, WINDOWS_CONTRACT)
    }
)
SUPPORTED_PLATFORMS: Final = tuple(PLATFORM_CONTRACTS)


def get_platform_contract(platform: str) -> WheelPlatformContract:
    """Return the registered contract for a CLI platform name."""

    try:
        return PLATFORM_CONTRACTS[platform]
    except KeyError as error:
        supported = ", ".join(SUPPORTED_PLATFORMS)
        raise ValueError(
            f"unsupported release-wheel platform {platform!r}; expected one of: "
            f"{supported}"
        ) from error
