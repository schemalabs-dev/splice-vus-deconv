"""Command-line interface for splice-vus-deconv."""

from __future__ import annotations
import argparse
import json
import sys
from .sashimi import SpliceJunction, aggregate_1mb_sashimi_window, score_donor_acceptor_shift
from .deconv import deconvolve_bulk_proportions, weight_celltype_psi
from .pae import compute_isoform_pae_penalty, prioritize_variant_panel


def run_benchmark() -> list[dict]:
    sig = [
        [12.4, 1.1, 0.4],
        [10.8, 0.9, 0.7],
        [0.8, 11.6, 1.2],
        [1.2, 9.9, 0.8],
        [0.5, 1.3, 13.1],
    ]
    bulk = [7.88, 6.91, 4.18, 3.84, 2.38]
    props = deconvolve_bulk_proportions(sig, bulk)

    apoa1_junc = [
        SpliceJunction("chr11", 116708180, 116708640, ref_reads=184, alt_reads=39, skip_ref_reads=12, skip_alt_reads=148)
    ]
    rnu4_junc = [
        SpliceJunction("chr2", 121530875, 121531420, ref_reads=210, alt_reads=62, skip_ref_reads=15, skip_alt_reads=131)
    ]
    wt_pae = [[1.2, 2.1, 2.8], [2.0, 1.1, 2.4], [2.9, 2.5, 1.0]]
    alt_pae = [[1.4, 14.8], [15.2, 1.3]]

    panel = [
        {
            "variant_id": "chr11:116708210:G>A (APOA1 c.47+1G>A)",
            "deconvolved_fractions": props,
            "sashimi_delta_psi": aggregate_1mb_sashimi_window(apoa1_junc, 116708210),
            "motif_bits_shift": score_donor_acceptor_shift("CAGGTAAGT", "CAGATAAGT"),
            "deconv_weighted_psi": weight_celltype_psi(props, [-0.74, -0.31, -0.08]),
            "isoform_pae_penalty": compute_isoform_pae_penalty(wt_pae, alt_pae, [0, 2]),
        },
        {
            "variant_id": "chr2:121530880:C>T (RNU4ATAC U12-donor)",
            "deconvolved_fractions": props,
            "sashimi_delta_psi": aggregate_1mb_sashimi_window(rnu4_junc, 121530880),
            "motif_bits_shift": score_donor_acceptor_shift("AAGGTAAGC", "AAGTTAAGC"),
            "deconv_weighted_psi": weight_celltype_psi(props, [-0.52, -0.48, -0.19]),
            "isoform_pae_penalty": 0.284,
        },
    ]
    return prioritize_variant_panel(panel)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="splice-vus-deconv")
    parser.add_argument("--benchmark", action="store_true", help="Run ClinVar/GTEx reference locus benchmark")
    args = parser.parse_args(argv)
    if args.benchmark:
        out = run_benchmark()
        print(json.dumps(out, indent=2))
        return 0
    parser.print_help()
    return 0


if __name__ == "__main__":
    sys.exit(main())
