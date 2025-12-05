"""Domain-level helper for managing collections of HLA alleles as haplotypes."""

from __future__ import annotations

from typing import Dict, Iterable, Iterator, Optional, Protocol, Sequence, Tuple, Union, runtime_checkable

from igen.core.enum import HlaLocus, HlaLocusEnum

from .hla_allele import HlaAllele, HlaAlleleProtocol

AlleleInput = Union[Iterable[HlaAlleleProtocol], "HlaHaplotypeProtocol", str]
LocusInput = Union[HlaLocus, HlaLocusEnum, str]


@runtime_checkable
class HlaHaplotypeProtocol(Protocol):
    """Structural typing contract for HLA haplotype models."""

    @classmethod
    def create(cls, alleles: AlleleInput) -> "HlaHaplotypeProtocol":
        """Construct a haplotype from iterables, strings, or another haplotype."""
        ...

    def get(self, locus: LocusInput, exact: bool = False) -> Optional[HlaAlleleProtocol]:
        """Retrieve the allele registered for the provided locus."""
        ...

    def set(self, locus: LocusInput, allele: HlaAlleleProtocol) -> "HlaHaplotypeProtocol":
        """Return a new haplotype with the allele stored for ``locus``."""
        ...

    def has(self, locus: LocusInput, exact: bool = False) -> bool:
        """Whether an allele is registered for the provided locus."""
        ...

    def swap(
        self, haplotype: "HlaHaplotypeProtocol", locus: LocusInput
    ) -> Tuple["HlaHaplotypeProtocol", "HlaHaplotypeProtocol"]:
        """Swap the allele at ``locus`` between two haplotypes."""
        ...

    def swap_all(
        self, haplotype: "HlaHaplotypeProtocol", loci: Sequence[LocusInput]
    ) -> Tuple["HlaHaplotypeProtocol", "HlaHaplotypeProtocol"]:
        """Swap all specified loci between this haplotype and another."""
        ...

    def clone(self) -> "HlaHaplotypeProtocol":
        """Return a shallow copy of the haplotype."""
        ...

    def concat(self, haplotype: "HlaHaplotypeProtocol") -> "HlaHaplotypeProtocol":
        """Return a new haplotype combining alleles from both inputs."""
        ...

    def __str__(self) -> str:
        """Render the haplotype as a '+'-joined list of allele strings."""
        ...

    def __iter__(self) -> Iterator[HlaAlleleProtocol]: ...

    @property
    def alleles(self) -> list[HlaAlleleProtocol]:
        """Return the collection of alleles in insertion order."""
        ...


def _coerce_locus(locus: LocusInput) -> HlaLocus:
    if isinstance(locus, HlaLocus):
        return locus

    coerced = HlaLocus.from_value(locus)
    if coerced is None:
        raise ValueError(f"Unknown locus: {locus}")

    return coerced


def _normalize_locus(locus: HlaLocus) -> HlaLocus:
    return HlaLocus.DRB345 if locus.is_drb345 else locus


class HlaHaplotype(HlaHaplotypeProtocol):
    """Immutable mapping of HLA loci to their corresponding alleles."""

    __slots__ = ("_allele_map",)

    def __init__(self, alleles: Iterable[HlaAlleleProtocol]):
        mapping: Dict[HlaLocus, HlaAlleleProtocol] = {}
        for allele in alleles:
            mapping[_normalize_locus(allele.locus)] = allele

        self._allele_map: Dict[HlaLocus, HlaAlleleProtocol] = mapping

    @classmethod
    def create(cls, alleles: AlleleInput) -> HlaHaplotype:
        """Factory mirroring the flexibility of the JavaScript implementation."""
        if isinstance(alleles, HlaHaplotype):
            return alleles.clone()

        if isinstance(alleles, str):
            allele_objects = [HlaAllele.from_string(token.strip()) for token in alleles.split("+") if token.strip()]
            return cls(allele_objects)

        return cls(alleles)

    def get(self, locus: LocusInput, exact: bool = False) -> Optional[HlaAlleleProtocol]:
        """Retrieve the allele registered for the given locus."""
        resolved = _coerce_locus(locus)

        if resolved.is_drb345:
            candidate = self._allele_map.get(HlaLocus.DRB345)
            if candidate is None:
                return None

            if exact and candidate.locus != resolved:
                return None

            return candidate

        return self._allele_map.get(resolved)

    def set(self, locus: LocusInput, allele: HlaAlleleProtocol) -> HlaHaplotype:
        """Return a new haplotype with the allele stored for ``locus``."""
        _locus = _coerce_locus(locus)
        norm_locus = _normalize_locus(_locus)

        if norm_locus != _normalize_locus(allele.locus):
            return self

        mapping = dict(self._allele_map)
        mapping[norm_locus] = allele
        return HlaHaplotype(mapping.values())

    def has(self, locus: LocusInput, exact: bool = False) -> bool:
        """Whether an allele is registered for the provided locus."""
        return self.get(locus, exact) is not None

    def swap(self, haplotype: "HlaHaplotypeProtocol", locus: LocusInput) -> Tuple["HlaHaplotype", "HlaHaplotype"]:
        """Swap the allele at ``locus`` between two haplotypes, returning clones."""
        first, second = self.clone(), haplotype.clone()

        allele1 = first.get(locus)
        allele2 = second.get(locus)

        if allele1 is None or allele2 is None:
            return first, second

        first = first.set(locus, allele2)
        second = second.set(locus, allele1)

        return first, second

    def swap_all(
        self, haplotype: "HlaHaplotypeProtocol", loci: Sequence[LocusInput]
    ) -> Tuple["HlaHaplotype", "HlaHaplotype"]:
        """Swap all loci listed in ``loci`` between this haplotype and another."""
        pair: Tuple[HlaHaplotype, HlaHaplotype] = (self, haplotype)

        for locus in loci:
            pair = pair[0].swap(pair[1], locus)

        return pair

    def clone(self) -> HlaHaplotype:
        """Return a shallow copy of the haplotype."""
        return HlaHaplotype(self.alleles)

    def concat(self, haplotype: "HlaHaplotypeProtocol") -> HlaHaplotype:
        """Return a new haplotype combining all alleles from both inputs."""
        return HlaHaplotype([*self.alleles, *haplotype.alleles])

    def __str__(self) -> str:
        """Render the haplotype as a '+'-joined list of allele strings."""
        return "+".join(allele.allele for allele in self.alleles)

    def __iter__(self) -> Iterator[HlaAlleleProtocol]:
        return iter(self.alleles)

    @property
    def alleles(self) -> list[HlaAlleleProtocol]:
        """Return the collection of alleles in insertion order."""
        return list(self._allele_map.values())
