"""
splice-vus-deconv
~~~~~~~~~~~~~~~~~
Cell-type-resolved prioritization of non-coding splice and regulatory variants
combining 1Mb Sashimi junction delta-PSI tracks, projected NNLS cellular
deconvolution, and AlphaFold isoform Pairwise Aligned Error (PAE) domain shift.
"""

from .sashimi import SpliceJunction, compute_delta_psi, score_donor_acceptor_shift
from .deconv import deconvolve_bulk_proportions, weight_celltype_psi
from .pae import compute_isoform_pae_penalty, prioritize_variant_panel

__version__ = "0.4.2"
__all__ = [
    "SpliceJunction",
    "compute_delta_psi",
    "score_donor_acceptor_shift",
    "deconvolve_bulk_proportions",
    "weight_celltype_psi",
    "compute_isoform_pae_penalty",
    "prioritize_variant_panel",
]
