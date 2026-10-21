"""Projected Gradient Non-Negative Least Squares (NNLS) Cellular Deconvolution.

Implements lightweight, dependency-free simplex-projected NNLS inspired by our
benchmark synthesis in Nucleic Acids Research (2024; doi:10.1093/nar/gkae267)
and perturbation consensus weighting from PINSPlus (Bioinformatics 2019).
"""

from __future__ import annotations
from typing import Sequence


def _project_simplex(vec: Sequence[float]) -> list[float]:
    """Euclidean projection onto the probability simplex {w >= 0, sum(w) = 1}."""
    n = len(vec)
    if n == 0:
        raise ValueError("Cannot project an empty vector")
    u = sorted(vec, reverse=True)
    cssv = 0.0
    rho = -1
    theta = 0.0
    for i in range(n):
        cssv += u[i]
        t = (cssv - 1.0) / (i + 1)
        if u[i] - t > 0:
            rho = i
            theta = t
    return [max(v - theta, 0.0) for v in vec]


def deconvolve_bulk_proportions(
    signature_matrix: Sequence[Sequence[float]],
    bulk_mixture: Sequence[float],
    max_iter: int = 400,
    tol: float = 1e-8,
) -> list[float]:
    """Estimate cell-type fractions w from bulk RNA-seq mixture y ~ S @ w on the simplex."""
    n_genes = len(signature_matrix)
    if n_genes == 0 or len(bulk_mixture) != n_genes:
        raise ValueError("Dimension mismatch between signature_matrix and bulk_mixture")
    n_cells = len(signature_matrix[0])
    if n_cells == 0:
        raise ValueError("signature_matrix must contain at least one cell type")

    # Compute Lipschitz step size 1 / ||S||_F^2
    frob_sq = sum(val * val for row in signature_matrix for val in row)
    lr = 1.0 / max(frob_sq, 1e-6)

    w = [1.0 / n_cells] * n_cells
    for _ in range(max_iter):
        # residual = S @ w - y
        grad = [0.0] * n_cells
        for g in range(n_genes):
            pred = sum(signature_matrix[g][c] * w[c] for c in range(n_cells))
            diff = pred - bulk_mixture[g]
            for c in range(n_cells):
                grad[c] += signature_matrix[g][c] * diff
        w_next = _project_simplex([w[c] - lr * grad[c] for c in range(n_cells)])
        delta = max(abs(w_next[c] - w[c]) for c in range(n_cells))
        w = w_next
        if delta < tol:
            break
    return [round(x, 6) for x in w]


def weight_celltype_psi(
    cell_proportions: Sequence[float], celltype_delta_psi: Sequence[float]
) -> float:
    """Compute lineage-weighted effective splice disruption across deconvolved cell states."""
    if len(cell_proportions) != len(celltype_delta_psi) or not cell_proportions:
        raise ValueError("Mismatched cell_proportions and celltype_delta_psi lengths")
    total = sum(abs(p) * abs(d) for p, d in zip(cell_proportions, celltype_delta_psi))
    return round(total, 6)
