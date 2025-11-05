"""Domain model and helpers for parsing and reasoning about HLA alleles."""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Optional, Protocol, runtime_checkable

from igen.shared.string_utils import is_blank

from igen.core.enum.hla_locus import HlaLocus

_LOCUS_DISPLAY_MAP: dict[HlaLocus, str] = {
    HlaLocus.DRB3: "B3",
    HlaLocus.DRB4: "B4",
    HlaLocus.DRB5: "B5",
    HlaLocus.DRB345: "",
}

_LOCUS_ALTERNATIVE_MAP: dict[str, HlaLocus] = {
    "B3": HlaLocus.DRB3,
    "B4": HlaLocus.DRB4,
    "B5": HlaLocus.DRB5,
}

_NEGATIVE_SPECIFICITIES = {"NEGATIVO", "NEGATIVE"}
_MISSING_SPECIFICITIES = {"AUSENTE", "MISSING"}
_SPECIFICITY_EXCEPTIONS = _NEGATIVE_SPECIFICITIES | _MISSING_SPECIFICITIES | {"?"}
_SPECIFICITY_REGEX = re.compile(r"^\d{2,}((:[A-Z]{2,})|(:\d{2,}[A-Z]?)*)$", re.IGNORECASE)


@runtime_checkable
class HlaAlleleProtocol(Protocol):
    """Structural typing contract for HLA allele models."""

    locus: HlaLocus
    specificity: str
    field_count: int
    mac_code: Optional[str]
    display_field_count: int

    @classmethod
    def from_string(cls, allele: str) -> HlaAlleleProtocol:
        """Parse a raw allele string into an instance."""
        ...

    @classmethod
    def is_valid_allele(cls, allele: str) -> bool:
        """Return ``True`` when the provided string is a well-formed HLA allele."""
        ...

    def clone(self) -> HlaAlleleProtocol:
        """Return a shallow copy of the allele."""
        ...

    def __str__(self) -> str:
        """Return the canonical allele string."""
        ...

    @property
    def allele(self) -> str:
        """Expose the canonical allele string as a property."""
        ...

    def display(self) -> str:
        """Return the allele string truncated to ``display_field_count`` fields."""
        ...

    def display_specificity(self) -> str:
        """Return the truncated specificity with any locus-specific adjustments applied."""
        ...

    def contains(self, other: "HlaAlleleProtocol | str") -> bool:
        """Return ``True`` when ``other`` is a less specific allele prefix of this one."""
        ...

    @property
    def has_mac_code(self) -> bool:
        """Whether the allele specificity ends with a MAC code suffix."""
        ...

    def with_display_field_count(self, value: int) -> HlaAlleleProtocol:
        """Return a new allele configured to display the specified number of fields."""
        ...

    @property
    def is_drb345(self) -> bool:
        """Whether the allele maps to one of the DRB3/4/5 loci."""
        ...

    @property
    def is_negative(self) -> bool:
        """Whether the specificity indicates a negative typing result."""
        ...

    @property
    def is_missing(self) -> bool:
        """Whether the specificity denotes a missing or absent allele."""
        ...

    @property
    def is_class_i(self) -> bool:
        """Whether the allele belongs to an HLA Class I locus."""
        ...

    @property
    def is_class_ii(self) -> bool:
        """Whether the allele belongs to an HLA Class II locus."""
        ...


def _is_valid_specificity(specificity: str) -> bool:
    if is_blank(specificity):
        return True

    normalized = specificity.upper()

    if normalized in _SPECIFICITY_EXCEPTIONS:
        return True

    return bool(_SPECIFICITY_REGEX.fullmatch(specificity))


def _parse_locus(locus_token: str) -> HlaLocus:
    parsed = HlaLocus.from_value(locus_token)
    if parsed is not None:
        return parsed

    alternative = _LOCUS_ALTERNATIVE_MAP.get(locus_token)
    if alternative is not None:
        return alternative

    raise ValueError(f"Unknown locus: {locus_token}")


def _extract_parts(allele: str) -> tuple[str, str]:
    normalized = allele.strip().upper()

    if normalized in _NEGATIVE_SPECIFICITIES or normalized in _MISSING_SPECIFICITIES:
        return "DRB345", normalized

    if "*" not in allele:
        raise ValueError(f"Allele must contain '*': {allele}")

    locus, specificity = allele.split("*", 1)
    return locus, specificity


@dataclass(frozen=True)
class HlaAllele(HlaAlleleProtocol):
    """Immutable representation of a parsed HLA allele string."""

    locus: HlaLocus
    specificity: str
    field_count: int
    mac_code: Optional[str]
    display_field_count: int = 2

    def __post_init__(self):
        """Validate specificity component immediately after initialization."""
        if not _is_valid_specificity(self.specificity):
            raise ValueError(f"Invalid specificity: {self.specificity}")

    @classmethod
    def from_string(cls, allele: str) -> HlaAllele:
        """Parse a raw allele string (e.g. ``'A*01:01:01'``) into an instance."""
        locus_token, specificity = _extract_parts(allele)
        locus = _parse_locus(locus_token)
        field_count = specificity.count(":") + 1
        mac_code = cls._extract_mac_code(specificity)
        return cls(locus=locus, specificity=specificity, field_count=field_count, mac_code=mac_code)

    @staticmethod
    def _extract_mac_code(specificity: str) -> Optional[str]:
        last_field = specificity.split(":")[-1]
        return last_field if last_field.isalpha() and last_field.isupper() and len(last_field) >= 2 else None

    @classmethod
    def is_valid_allele(cls, allele: str) -> bool:
        """Return ``True`` when the provided string is a well-formed HLA allele."""
        try:
            locus_token, specificity = _extract_parts(allele)
        except ValueError:
            return False

        try:
            _parse_locus(locus_token)
        except ValueError:
            return False

        return _is_valid_specificity(specificity)

    def clone(self) -> HlaAllele:
        """Return a shallow copy of the current allele."""
        return HlaAllele(self.locus, self.specificity, self.field_count, self.mac_code, self.display_field_count)

    def __str__(self) -> str:
        """Return the canonical allele string (e.g. ``'A*01:01'``)."""
        return f"{self.locus.value.value}*{self.specificity}"

    @property
    def allele(self) -> str:
        """Expose the canonical allele string as a property."""
        return str(self)

    def display(self) -> str:
        """Return the allele string truncated to ``display_field_count`` fields."""
        return self._reduce_specificity(str(self))

    def display_specificity(self) -> str:
        """Return the truncated specificity, applying DRB345 aliases when needed."""
        reduced = self._reduce_specificity(self.specificity)

        if not self.is_drb345:
            return reduced

        alias = _LOCUS_DISPLAY_MAP.get(self.locus)
        if is_blank(alias):
            return reduced

        return f"{alias}*{reduced}"

    def _reduce_specificity(self, specificity: str) -> str:
        if re.search(r"\d[A-Z]$", specificity):
            return specificity

        parts = specificity.split(":")
        return ":".join(parts[: self.display_field_count])

    def contains(self, other: HlaAlleleProtocol | str) -> bool:
        """Return ``True`` when ``other`` is a less specific allele prefix of this one."""
        candidate = HlaAllele.from_string(other) if isinstance(other, str) else other

        if self.has_mac_code or candidate.has_mac_code:
            return False

        return str(self).startswith(str(candidate))

    @property
    def has_mac_code(self) -> bool:
        """Whether the allele specificity ends with a MAC code suffix."""
        return not is_blank(self.mac_code)

    def with_display_field_count(self, value: int) -> HlaAllele:
        """Return a new allele configured to display the specified number of fields."""
        return HlaAllele(self.locus, self.specificity, self.field_count, self.mac_code, value)

    @property
    def is_drb345(self) -> bool:
        """Whether the allele maps to one of the DRB3/4/5 loci."""
        return self.locus.is_drb345

    @property
    def is_negative(self) -> bool:
        """Whether the specificity indicates a negative typing result."""
        return self.specificity.upper() in _NEGATIVE_SPECIFICITIES

    @property
    def is_missing(self) -> bool:
        """Whether the specificity denotes a missing or absent allele."""
        return is_blank(self.specificity) or self.specificity.upper() in _MISSING_SPECIFICITIES

    @property
    def is_class_i(self) -> bool:
        """Whether the allele belongs to an HLA Class I locus."""
        return self.locus.is_class_i

    @property
    def is_class_ii(self) -> bool:
        """Whether the allele belongs to an HLA Class II locus."""
        return self.locus.is_class_ii
