import unittest
from splice_vus_deconv.sashimi import (
    SpliceJunction,
    compute_delta_psi,
    score_donor_acceptor_shift,
    aggregate_1mb_sashimi_window,
)
from splice_vus_deconv.deconv import deconvolve_bulk_proportions, weight_celltype_psi
from splice_vus_deconv.pae import compute_isoform_pae_penalty, prioritize_variant_panel
from splice_vus_deconv.cli import run_benchmark


class TestSpliceVusDeconv(unittest.TestCase):
    def test_canonical_donor_disruption_negative_bits(self):
        shift = score_donor_acceptor_shift("CAGGTAAGT", "CAGATAAGT")
        self.assertLess(shift, -6.0)

    def test_delta_psi_exon_skipping_direction(self):
        j = SpliceJunction("chr11", 1000, 1800, ref_reads=200, alt_reads=25, skip_ref_reads=10, skip_alt_reads=160)
        dpsi = compute_delta_psi(j)
        self.assertLess(dpsi, -0.75)

    def test_1mb_window_distance_cutoff(self):
        j = SpliceJunction("chr11", 100_000, 102_000, ref_reads=150, alt_reads=20, skip_ref_reads=10, skip_alt_reads=120)
        near = aggregate_1mb_sashimi_window([j], 100_050)
        outside = aggregate_1mb_sashimi_window([j], 1_500_000)
        self.assertGreater(near, 0.5)
        self.assertEqual(outside, 0.0)

    def test_nnls_simplex_deconvolution_recovers_true_weights(self):
        sig = [
            [10.0, 0.5, 0.2],
            [0.3, 12.0, 0.4],
            [0.2, 0.6, 9.5],
            [8.0, 2.0, 1.0],
        ]
        true_w = [0.6, 0.25, 0.15]
        bulk = [sum(sig[g][c] * true_w[c] for c in range(3)) for g in range(4)]
        est = deconvolve_bulk_proportions(sig, bulk)
        self.assertAlmostEqual(sum(est), 1.0, places=5)
        for e, t in zip(est, true_w):
            self.assertAlmostEqual(e, t, places=2)

    def test_pae_destabilization_and_ranking(self):
        wt = [[1.0, 2.0, 2.2], [2.0, 1.0, 1.9], [2.2, 1.9, 1.0]]
        alt = [[1.2, 18.5], [19.0, 1.1]]
        pen = compute_isoform_pae_penalty(wt, alt, [0, 2])
        self.assertGreater(pen, 0.45)
        bench = run_benchmark()
        self.assertEqual(len(bench), 2)
        self.assertGreater(bench[0]["pathogenicity_posterior"], bench[1]["pathogenicity_posterior"])


if __name__ == "__main__":
    unittest.main()
