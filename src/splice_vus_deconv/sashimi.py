"""1Mb Regulatory Window Sashimi Splice-Junction & Donor/Acceptor Perturbation."""

from __future__ import annotations
from dataclasses import dataclass
import math
from typing import Sequence


@dataclass(frozen=True)
class SpliceJunction:
    chrom: str
    donor_pos: int
    acceptor_pos: int
    ref_reads: float
    alt_reads: float
    skip_ref_reads: float
    skip_alt_reads: float

    @property
    def span_bp(self) -> int:
        return abs(self.acceptor_pos - self.donor_pos)


# Position-weight log-odds proxy for canonical U2 5' donor (MAG|GTRAGT)
_DONOR_CONSENSUS = (
    {"A": 0.35, "C": 0.30, "G": 0.20, "T": 0.15},
    {"A": 0.62, "C": 0.12, "G": 0.14, "T": 0.12},
    {"A": 0.09, "C": 0.04, "G": 0.80, "T": 0.07},
    {"A": 0.005, "C": 0.005, "G": 0.985, "T": 0.005},
    {"A": 0.005, "C": 0.005, "G": 0.005, "T": 0.985},
    {"A": 0.54, "C": 0.03, "G": 0.39, "T": 0.04},
    {"A": 0.71, "C": 0.08, "G": 0.12, "T": 0.09},
    {"A": 0.07, "C": 0.06, "G": 0.81, "T": 0.06},
    {"A": 0.16, "C": 0.15, "G": 0.22, "T": 0.47},
)


def compute_delta_psi(junction: SpliceJunction, pseudocount: float = 1.0) -> float:
    """Compute exon-inclusion Delta-PSI from Sashimi inclusion and skipping junction arcs."""
    if pseudocount <= 0:
        raise ValueError("pseudocount must be strictly positive")
    psi_ref = (junction.ref_reads + pseudocount) / (
        junction.ref_reads + 2.0 * junction.skip_ref_reads + 2.0 * pseudocount
    )
    psi_alt = (junction.alt_reads + pseudocount) / (
        junction.alt_reads + 2.0 * junction.skip_alt_reads + 2.0 * pseudocount
    )
    return round(psi_alt - psi_ref, 6)


def score_donor_acceptor_shift(ref_9mer: str, alt_9mer: str) -> float:
    """Return log2 information-content disruption across a 9-nt 5' splice donor motif."""
    r = ref_9mer.upper().replace("U", "T")
    a = alt_9mer.upper().replace("U", "T")
    if len(r) != 9 or len(a) != 9:
        raise ValueError("5' donor motif sequences must be exactly 9 nucleotides")
    bg = 0.25
    ll_ref = 0.0
    ll_alt = 0.0
    for idx, pwm in enumerate(_DONOR_CONSENSUS):
        if r[idx] not in pwm or a[idx] not in pwm:
            raise ValueError(f"Invalid nucleotide at position {idx}")
        ll_ref += math.log2(pwm[r[idx]] / bg)
        ll_alt += math.log2(pwm[a[idx]] / bg)
    return round(ll_alt - ll_ref, 6)


def aggregate_1mb_sashimi_window(
    junctions: Sequence[SpliceJunction], variant_pos: int, decay_scale_bp: float = 250_000.0
) -> float:
    """Aggregate distance-attenuated |Delta-PSI| across a 1Mb cis-regulatory interval."""
    if not junctions:
        return 0.0
    best = 0.0
    for j in junctions:
        nearest = min(abs(variant_pos - j.donor_pos), abs(variant_pos - j.acceptor_pos))
        if nearest > 1_000_000:
            continue
        atten = math.exp(-nearest / decay_scale_bp)
        score = abs(compute_delta_psi(j)) * atten
        if score > best:
            best = score
    return round(best, 6)
