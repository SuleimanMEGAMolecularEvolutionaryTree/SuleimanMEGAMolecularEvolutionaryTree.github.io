import sys
import os
import unittest
import numpy as np

# Add scripts directory to path
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"))

from phylo_pipeline import (
    jc69_distance,
    k2p_distance,
    compute_alignment_counts,
    parse_meg_file,
    canonical_split,
    get_splits_from_joins,
    format_newick,
    neighbor_joining,
    bootstrap_neighbor_joining
)


class TestPhyloModels(unittest.TestCase):

    def test_jc69_zero_diff(self):
        """JC69 distance should be 0 when there are no differences."""
        self.assertEqual(jc69_distance(0.0), 0.0)

    def test_jc69_standard_val(self):
        """Test JC69 distance with a typical proportion of differences."""
        # p = 0.1 -> d = -0.75 * ln(1 - 4/3 * 0.1) approx 0.1073256
        d = jc69_distance(0.1)
        self.assertAlmostEqual(d, 0.1073256, places=6)

    def test_jc69_infinite_limit(self):
        """JC69 distance should return infinity when base differences are >= 0.75."""
        self.assertEqual(jc69_distance(0.75), np.inf)
        self.assertEqual(jc69_distance(0.80), np.inf)

    def test_k2p_zero_diff(self):
        """K2P distance should be 0 when differences are 0."""
        self.assertEqual(k2p_distance(0.0, 0.0), 0.0)

    def test_k2p_standard_val(self):
        """Test K2P distance with typical transition/transversion counts."""
        # p = 0.05, q = 0.02 -> d = -0.5 * ln(1 - 2*0.05 - 0.02) - 0.25 * ln(1 - 2*0.02)
        # = -0.5 * ln(0.88) - 0.25 * ln(0.96) approx 0.07412
        d = k2p_distance(0.05, 0.02)
        self.assertAlmostEqual(d, 0.0741219, places=6)

    def test_k2p_infinite_limit(self):
        """K2P distance should return infinity for mathematically invalid values."""
        self.assertEqual(k2p_distance(0.4, 0.3), np.inf)  # 1 - 2*0.4 - 0.3 = -0.1 (<= 0)
        self.assertEqual(k2p_distance(0.1, 0.6), np.inf)  # 1 - 2*0.6 = -0.2 (<= 0)

    def test_alignment_counts(self):
        """Test alignment classification for matches, transitions, transversions, and gaps."""
        # Purines: A, G; Pyrimidines: C, T
        # Col 1: A / A -> Match
        # Col 2: A / G -> Transition (Purines)
        # Col 3: C / T -> Transition (Pyrimidines)
        # Col 4: A / C -> Transversion (Purine/Pyrim)
        # Col 5: G / G -> Match
        # Col 6: - / T -> Gap (Ignored in comparison sites)
        aligned_a = "AACAG-"
        aligned_b = "AGTCGT"

        counts = compute_alignment_counts(aligned_a, aligned_b)

        # Total sites without gaps should be 5
        self.assertEqual(counts["n_sites"], 5)
        # Matches: A-A (col 1), G-G (col 5) -> 2 matches
        # Transitions: A-G (col 2), C-T (col 4) -> 2 transitions
        # Transversions: A-C (col 3) -> 1 transversion
        self.assertEqual(counts["n_P"], 2)
        self.assertEqual(counts["n_Q"], 1)
        self.assertEqual(counts["p"], 2/5)
        self.assertEqual(counts["q"], 1/5)


class TestExtendedPhyloPipeline(unittest.TestCase):

    def test_canonical_split(self):
        """Canonical split should normalize regardless of which side is chosen."""
        taxa = ["A", "B", "C", "D"]
        # Split dividing {A, B} from {C, D}
        s1 = canonical_split({"A", "B"}, taxa)
        s2 = canonical_split({"C", "D"}, taxa)
        self.assertEqual(s1, s2)
        self.assertEqual(s1, frozenset({"C", "D"}))

    def test_format_newick_with_bootstrap(self):
        """Newick formatter should embed bootstrap support numbers properly."""
        taxa = ["A", "B", "C", "D"]
        joins = [
            ("A", "B", 0.05, 0.05),
            ("(A,B)", "C", 0.10, 0.10),
            ("((A,B),C)", "D", 0.20, 0.20)
        ]
        can_split = canonical_split({"A", "B"}, taxa)
        support = {can_split: 99.0}
        nwk = format_newick(joins, all_taxa=taxa, split_support=support)
        self.assertIn(")99:0.10000", nwk)
        self.assertTrue(nwk.endswith(";"))

    def test_parse_meg_file(self):
        """MEGA file parser should correctly extract valid non-empty sequences."""
        meg_path = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "data", "raw", "animal_alignments", "animals.meg"
        )
        if os.path.exists(meg_path):
            seqs = parse_meg_file(meg_path)
            self.assertGreaterEqual(len(seqs), 8)
            self.assertIn("Homo_sapiens_9QB3_b", seqs)
            self.assertEqual(len(seqs["Homo_sapiens_9QB3_b"]), 839)

    def test_bootstrap_neighbor_joining(self):
        """Bootstrap resampling should calculate support values between 0 and 100%."""
        taxa = ["A", "B", "C", "D"]
        aligned_pairs = {
            (0, 1): ("ACGTACGTACGT", "ACGTACGTACGA"),
            (0, 2): ("ACGTACGTACGT", "CCGTACGTACGA"),
            (0, 3): ("ACGTACGTACGT", "GGGTACGTACGA"),
            (1, 2): ("ACGTACGTACGA", "CCGTACGTACGA"),
            (1, 3): ("ACGTACGTACGA", "GGGTACGTACGA"),
            (2, 3): ("CCGTACGTACGA", "GGGTACGTACGA")
        }
        support = bootstrap_neighbor_joining(taxa, aligned_pairs, model="k2p", n_boot=20)
        self.assertIsInstance(support, dict)
        for pct in support.values():
            self.assertGreaterEqual(pct, 0.0)
            self.assertLessEqual(pct, 100.0)


if __name__ == "__main__":
    unittest.main()
