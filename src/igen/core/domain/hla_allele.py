"""Domain model and helpers for parsing and reasoning about HLA alleles."""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Optional, Protocol, runtime_checkable

from igen.shared.string_utils import is_blank

from igen.core.enum import HlaLocus

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
    suffix: Optional[str]

    @classmethod
    def from_string(cls, allele: str) -> "HlaAlleleProtocol":
        """Parse a raw allele string into an instance."""
        ...

    @classmethod
    def is_valid_allele(cls, allele: str) -> bool:
        """Return ``True`` when the provided string is a well-formed HLA allele."""
        ...

    def clone(self) -> "HlaAlleleProtocol":
        """Return a shallow copy of the allele."""
        ...

    def __str__(self) -> str:
        """Return the canonical allele string."""
        ...

    @property
    def allele(self) -> str:
        """Expose the canonical allele string as a property."""
        ...

    def display(self, force_truncate: bool = False, keep_suffix: bool = False) -> str:
        """Return the allele string truncated to ``display_field_count`` fields.

        When ``force_truncate`` is ``True`` the return value is always truncated,
        even if the original specificity already ends with a lettered suffix.
        When ``keep_suffix`` is ``True`` the specificity suffix is always kept.
        """
        ...

    def display_specificity(self, force_truncate: bool = False, keep_suffix: bool = False) -> str:
        """Return the truncated specificity with locus-specific adjustments.

        When ``force_truncate`` is ``True`` the specificity is always truncated,
        even if the original string would normally be left intact.
        When ``keep_suffix`` is ``True`` the specificity suffix is always kept.
        """
        ...

    def contains(self, other: "HlaAlleleProtocol | str") -> bool:
        """Return ``True`` when ``other`` is a less specific allele prefix of this one."""
        ...

    @property
    def has_mac_code(self) -> bool:
        """Whether the allele specificity ends with a MAC code suffix."""
        ...

    def with_display_field_count(self, value: int) -> "HlaAlleleProtocol":
        """Return a new allele configured to display the specified number of fields."""
        ...

    @property
    def has_suffix(self) -> bool:
        """Whether the allele has a suffix."""
        ...

    def without_suffix(self) -> "HlaAlleleProtocol":
        """Return a copy of the allele without any suffix."""
        ...

    def as_resolution(self, n_field: int, keep_suffix: bool = True) -> "HlaAlleleProtocol":
        """Return a copy of the allele truncated to ``n_field`` fields."""
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

    @property
    def allelic_group(self) -> str:
        """Return the allele string truncated to its first (allelic group) field."""
        ...

    @property
    def has_suffix(self) -> bool:
        """Whether the allele has a suffix."""
        ...

    @property
    def is_null(self) -> bool:
        """Whether the allele carries the null-expression ``N`` suffix."""
        ...

    @property
    def is_low(self) -> bool:
        """Whether the allele carries the low-expression ``L`` suffix."""
        ...

    @property
    def is_questionable(self) -> bool:
        """Whether the allele carries the questionable-expression ``Q`` suffix."""
        ...

    @property
    def is_low_resolution(self) -> bool:
        """Whether the allele is typed only to the low-resolution (one-field) level."""
        ...

    @property
    def is_mid_resolution(self) -> bool:
        """Whether the allele is typed to two fields with an appended MAC code."""
        ...

    @property
    def is_high_resolution(self) -> bool:
        """Whether the allele is typed to at least two fields without a MAC code."""
        ...


def _is_valid_specificity(specificity: str) -> bool:
    normalized = specificity.upper()

    if normalized in _SPECIFICITY_EXCEPTIONS:
        return True

    return bool(_SPECIFICITY_REGEX.fullmatch(specificity))


def _parse_locus(locus: str) -> HlaLocus:
    norm_locus = locus.upper()
    parsed = HlaLocus.from_value(norm_locus)
    if parsed is not None:
        return parsed

    alternative = _LOCUS_ALTERNATIVE_MAP.get(norm_locus)
    if alternative is not None:
        return alternative

    raise ValueError(f"Unknown locus: {norm_locus}")


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
    suffix: Optional[str] = None

    def __post_init__(self):
        """Validate specificity component immediately after initialization."""
        if not _is_valid_specificity(self.specificity):
            raise ValueError(f"Invalid specificity: {self.specificity}")

    @classmethod
    def from_string(cls, allele: str) -> "HlaAllele":
        """Parse a raw allele string (e.g. ``'A*01:01:01'``) into an instance."""
        locus_token, specificity = _extract_parts(allele)
        locus = _parse_locus(locus_token)
        field_count = cls.get_field_count(specificity)
        mac_code = cls.get_mac_code(specificity)
        suffix = cls.get_suffix(specificity)
        return cls(locus=locus, specificity=specificity, field_count=field_count, mac_code=mac_code, suffix=suffix)

    @staticmethod
    def get_locus_str(allele: str) -> str:
        locus, _ = _extract_parts(allele)
        return locus

    @staticmethod
    def get_locus(allele: str) -> HlaLocus:
        return _parse_locus(HlaAllele.get_locus_str(allele))

    @staticmethod
    def get_specificity(allele: str) -> str:
        _, specificity = _extract_parts(allele)
        return specificity

    @staticmethod
    def get_field_count(specificity: str) -> int:
        return specificity.count(":") + 1

    @staticmethod
    def get_mac_code(specificity: str) -> Optional[str]:
        last_field = specificity.split(":")[-1]
        return last_field.upper() if len(last_field) >= 2 and last_field.isalpha() else None

    @staticmethod
    def get_suffix(specificity: str) -> Optional[str]:
        last_field = specificity.split(":")[-1]
        return (
            last_field[-1].upper()
            if len(last_field) >= 3 and not last_field.isalpha() and last_field[-1].isalpha()
            else None
        )

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

    def clone(self) -> "HlaAllele":
        """Return a shallow copy of the current allele."""
        return HlaAllele(
            self.locus, self.specificity, self.field_count, self.mac_code, self.display_field_count, self.suffix
        )

    def __str__(self) -> str:
        """Return the canonical allele string (e.g. ``'A*01:01'``)."""
        return f"{self.locus.value.value}*{self.specificity}"

    @property
    def allele(self) -> str:
        """Expose the canonical allele string as a property."""
        return str(self)

    def display(self, force_truncate: bool = False, keep_suffix: bool = False) -> str:
        """Return the allele string truncated to ``display_field_count`` fields.

        Setting ``force_truncate`` forces truncation even when the specificity ends
        with a trailing letter that we would normally preserve.
        """
        return self._reduce_specificity(str(self), force_truncate, keep_suffix)

    def display_specificity(self, force_truncate: bool = False, keep_suffix: bool = False) -> str:
        """Return the truncated specificity, applying DRB345 aliases when needed.

        Setting ``force_truncate`` forces truncation even when the specificity ends
        with a trailing letter that would normally keep the string intact.
        """
        reduced = self._reduce_specificity(self.specificity, force_truncate, keep_suffix)

        if not self.is_drb345:
            return reduced

        alias = _LOCUS_DISPLAY_MAP.get(self.locus)
        if is_blank(alias):
            return reduced

        return f"{alias}*{reduced}"

    def _reduce_specificity(self, specificity: str, force_truncate: bool = False, keep_suffix: bool = False) -> str:
        if not force_truncate and self.has_suffix:
            return specificity

        parts = specificity.split(":")
        reduced = ":".join(parts[: self.display_field_count])

        return reduced + self.suffix if keep_suffix else reduced

    def contains(self, other: "HlaAlleleProtocol | str") -> bool:
        """Return ``True`` when ``other`` is a less specific allele prefix of this one."""
        candidate = HlaAllele.from_string(other) if isinstance(other, str) else other

        if self.has_mac_code or candidate.has_mac_code:
            return False

        return str(self).startswith(str(candidate))

    @property
    def has_mac_code(self) -> bool:
        """Whether the allele specificity ends with a MAC code suffix."""
        return not is_blank(self.mac_code)

    def with_display_field_count(self, value: int) -> "HlaAllele":
        """Return a new allele configured to display the specified number of fields."""
        return HlaAllele(self.locus, self.specificity, self.field_count, self.mac_code, value, self.suffix)

    def without_suffix(self) -> "HlaAllele":
        """Return a copy of the allele without any suffix."""
        if not self.has_suffix:
            return self.clone()

        return HlaAllele(
            self.locus,
            self.specificity[:-1],
            self.field_count,
            self.mac_code,
            self.display_field_count,
            None,
        )

    def as_resolution(self, n_field: int, keep_suffix: bool = True) -> "HlaAllele":
        """Return a copy of the allele truncated to ``n_field`` fields."""
        if n_field < 1:
            raise ValueError("n_field must be >= 1")

        if n_field > self.field_count:
            return self.clone()

        specificity = self._reduce_specificity(self.specificity, force_truncate=True)
        suffix = self.suffix if keep_suffix and self.has_suffix else None

        return HlaAllele(
            self.locus,
            specificity,
            n_field,
            HlaAllele.get_mac_code(specificity),
            min(self.display_field_count, n_field),
            suffix,
        )

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

    @property
    def allelic_group(self) -> str:
        """Return the allele string truncated to its first (allelic group) field."""
        return self.with_display_field_count(1).display(force_truncate=True)

    @property
    def has_suffix(self) -> bool:
        """Whether the allele has a suffix."""
        return not is_blank(self.suffix)

    @property
    def is_null(self) -> bool:
        """Whether the allele carries the null-expression ``N`` suffix."""
        return self.suffix == "N"

    @property
    def is_low(self) -> bool:
        """Whether the allele carries the low-expression ``L`` suffix."""
        return self.suffix == "L"

    @property
    def is_questionable(self) -> bool:
        """Whether the allele carries the questionable-expression ``Q`` suffix."""
        return self.suffix == "Q"

    @property
    def is_low_resolution(self) -> bool:
        """Whether the allele is typed only to the low-resolution (one-field) level."""
        return self.field_count == 1

    @property
    def is_mid_resolution(self) -> bool:
        """Whether the allele is typed to two fields with an appended MAC code."""
        return self.field_count == 2 and self.has_mac_code

    @property
    def is_high_resolution(self) -> bool:
        """Whether the allele is typed to at least two fields without a MAC code."""
        return self.field_count >= 2 and not self.has_mac_code
