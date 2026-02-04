from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import ClassVar

from .hla_locus_class import HlaLocusClass


class HlaLocusGroupEnum(Enum):
    """Groups of HLA loci organized by shared biological characteristics."""

    A = "A"
    B = "B"
    C = "C"
    DR = "DR"
    DQ = "DQ"
    DP = "DP"


@dataclass(frozen=True)
class HlaLocusGroup:
    value: HlaLocusGroupEnum
    has_alpha: bool
    has_beta: bool
    classification: HlaLocusClass

    A: ClassVar["HlaLocusGroup"]
    B: ClassVar["HlaLocusGroup"]
    C: ClassVar["HlaLocusGroup"]
    DR: ClassVar["HlaLocusGroup"]
    DQ: ClassVar["HlaLocusGroup"]
    DP: ClassVar["HlaLocusGroup"]

    GROUPS: ClassVar[tuple["HlaLocusGroup", ...]]

    @property
    def is_class_i(self) -> bool:
        return self.classification == HlaLocusClass.I

    @property
    def is_class_ii(self) -> bool:
        return self.classification == HlaLocusClass.II


HlaLocusGroup.A = HlaLocusGroup(HlaLocusGroupEnum.A, True, False, HlaLocusClass.I)
HlaLocusGroup.B = HlaLocusGroup(HlaLocusGroupEnum.B, True, False, HlaLocusClass.I)
HlaLocusGroup.C = HlaLocusGroup(HlaLocusGroupEnum.C, True, False, HlaLocusClass.I)
HlaLocusGroup.DR = HlaLocusGroup(HlaLocusGroupEnum.DR, False, True, HlaLocusClass.II)
HlaLocusGroup.DQ = HlaLocusGroup(HlaLocusGroupEnum.DQ, True, True, HlaLocusClass.II)
HlaLocusGroup.DP = HlaLocusGroup(HlaLocusGroupEnum.DP, True, True, HlaLocusClass.II)

HlaLocusGroup.GROUPS = (
    HlaLocusGroup.A,
    HlaLocusGroup.B,
    HlaLocusGroup.C,
    HlaLocusGroup.DR,
    HlaLocusGroup.DQ,
    HlaLocusGroup.DP,
)
