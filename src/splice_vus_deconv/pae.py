"""AlphaFold Isoform Pairwise Aligned Error (PAE) & Integrated VUS Triage."""

from __future__ import annotations
import math
from typing import Sequence


def compute_isoform_pae_penalty(
    wt_pae_matrix: Sequence[Sequence[float]],
    alt_pae_matrix: Sequence[Sequence[float]],
    retained_indices: Sequence[int],
    max_pae_angstrom: float = 31.75,
) -> float:
    """Quantify inter-residue packing destabilization induced by exon skipping/retention."""
    if not retained_indices:
        return 1.0
    deltas: list[float] = []
    for r_idx, wt_r in enumerate(retained_indices):
        for c_idx, wt_c in enumerate(retained_indices):
            if r_idx == c_idx:
                continue
            wt_val = wt_pae_matrix[wt_r][wt_c]
            alt_val = alt_pae_matrix[r_idx][c_idx]
            deltas.append(max(0.0, alt_val - wt_val) / max_pae_angstrom)
    if not deltas:
        return 0.0
    mean_shift = sum(deltas) / len(deltas)
    return round(min(1.0, mean_shift), 6)


def prioritize_variant_panel(records: Sequence[dict]) -> list[dict]:
    """Rank non-coding VUS records by joint Sashimi, Deconvolution, and PAE impact."""
    ranked = []
    for rec in records:
        sashimi = abs(float(rec["sashimi_delta_psi"]))
        motif = abs(float(rec.get("motif_bits_shift", 0.0))) / 10.0
        deconv = abs(float(rec["deconv_weighted_psi"]))
        pae = float(rec["isoform_pae_penalty"])
        logit = -2.4 + 3.2 * sashimi + 1.8 * deconv + 2.1 * pae + 0.9 * min(motif, 1.5)
        posterior = 1.0 / (1.0 + math.exp(-logit))
        enriched = dict(rec)
        enriched["pathogenicity_posterior"] = round(posterior, 4)
        ranked.append(enriched)
    ranked.sort(key=lambda r: r["pathogenicity_posterior"], reverse=True)
    return ranked
