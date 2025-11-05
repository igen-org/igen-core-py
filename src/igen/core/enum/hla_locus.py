from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import ClassVar, Iterable, Optional

from .hla_locus_chain import HlaLocusChain
from .hla_locus_group import HlaLocusGroup


class HlaLocusEnum(Enum):
    """Canonical HLA loci recognized by the core domain."""

    A = "A"
    B = "B"
    C = "C"
    DRB1 = "DRB1"
    DRB3 = "DRB3"
    DRB4 = "DRB4"
    DRB5 = "DRB5"
    DRB345 = "DRB345"
    DQB1 = "DQB1"
    DQA1 = "DQA1"
    DPB1 = "DPB1"
    DPA1 = "DPA1"


@dataclass(frozen=True)
class HlaLocus:
    value: HlaLocusEnum
    group: HlaLocusGroup
    chain_type: HlaLocusChain

    _registry: ClassVar[list["HlaLocus"]] = []

    def __post_init__(self):
        self.__class__._registry.append(self)

    def __str__(self) -> str:
        return self.value.value

    @property
    def is_drb345(self) -> bool:
        return self.value in {HlaLocusEnum.DRB345, HlaLocusEnum.DRB3, HlaLocusEnum.DRB4, HlaLocusEnum.DRB5}

    @property
    def is_alpha(self) -> bool:
        return self.chain_type == HlaLocusChain.ALPHA

    @property
    def is_beta(self) -> bool:
        return self.chain_type == HlaLocusChain.BETA

    @property
    def is_class_i(self) -> bool:
        return self.group.is_class_i

    @property
    def is_class_ii(self) -> bool:
        return self.group.is_class_ii

    @classmethod
    def _loci(cls) -> Iterable["HlaLocus"]:
        return cls._registry

    @classmethod
    def from_value(cls, value: Optional[str | HlaLocusEnum]) -> Optional["HlaLocus"]:
        if value is None:
            return None

        value_str = value.value if isinstance(value, HlaLocusEnum) else str(value)
        for locus in cls._registry:
            if locus.value.value == value_str:
                return locus
        return None


HlaLocus.A = HlaLocus(HlaLocusEnum.A, HlaLocusGroup.A, HlaLocusChain.ALPHA)
HlaLocus.B = HlaLocus(HlaLocusEnum.B, HlaLocusGroup.B, HlaLocusChain.ALPHA)
HlaLocus.C = HlaLocus(HlaLocusEnum.C, HlaLocusGroup.C, HlaLocusChain.ALPHA)
HlaLocus.DRB1 = HlaLocus(HlaLocusEnum.DRB1, HlaLocusGroup.DR, HlaLocusChain.BETA)
HlaLocus.DRB3 = HlaLocus(HlaLocusEnum.DRB3, HlaLocusGroup.DR, HlaLocusChain.BETA)
HlaLocus.DRB4 = HlaLocus(HlaLocusEnum.DRB4, HlaLocusGroup.DR, HlaLocusChain.BETA)
HlaLocus.DRB5 = HlaLocus(HlaLocusEnum.DRB5, HlaLocusGroup.DR, HlaLocusChain.BETA)
HlaLocus.DRB345 = HlaLocus(HlaLocusEnum.DRB345, HlaLocusGroup.DR, HlaLocusChain.BETA)
HlaLocus.DQB1 = HlaLocus(HlaLocusEnum.DQB1, HlaLocusGroup.DQ, HlaLocusChain.BETA)
HlaLocus.DQA1 = HlaLocus(HlaLocusEnum.DQA1, HlaLocusGroup.DQ, HlaLocusChain.ALPHA)
HlaLocus.DPB1 = HlaLocus(HlaLocusEnum.DPB1, HlaLocusGroup.DP, HlaLocusChain.BETA)
HlaLocus.DPA1 = HlaLocus(HlaLocusEnum.DPA1, HlaLocusGroup.DP, HlaLocusChain.ALPHA)

HlaLocus.LOCI = tuple(HlaLocus._loci())
