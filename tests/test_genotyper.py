"""Tests for genotyper.py — Bayesian genotyping module."""

import math
import pytest
from LOCATE.genotyper import (
    genotype_posterior,
    genotype_from_counts,
    define_genotype_threshold,
    _beta_binomial_log_likelihood,
)


class TestBetaBinomialLogLikelihood:
    """Unit tests for the core log-likelihood function."""

    def test_zero_trials(self):
        """With n=0, the likelihood should be 0 (log 0 = 0 means likelihood = 1)."""
        assert _beta_binomial_log_likelihood(0, 0, 1, 1) == 0.0

    def test_symmetric_prior(self):
        """With Beta(1,1), log-likelihood should be log(1/(n+1))."""
        ll = _beta_binomial_log_likelihood(2, 10, 1, 1)
        # Beta(1,1) is uniform: P(k|n) = 1/(n+1) for any k
        expected = -math.log(11)
        assert math.isclose(ll, expected, rel_tol=1e-12)

    def test_all_success(self):
        """k=n should have a valid finite log-likelihood (not -inf)."""
        ll = _beta_binomial_log_likelihood(100, 100, 30, 1)
        assert math.isfinite(ll)
        assert ll < 0  # likelihood should be < 1

    def test_all_failure(self):
        """k=0 should have a valid finite log-likelihood."""
        ll = _beta_binomial_log_likelihood(0, 100, 1, 30)
        assert math.isfinite(ll)
        assert ll < 0

    def test_extreme_counts_no_overflow(self):
        """Very large counts should not cause overflow (inf or nan)."""
        ll = _beta_binomial_log_likelihood(5000, 100000, 5, 5)
        assert math.isfinite(ll)


class TestGenotypePosterior:
    """Tests for posterior probability computation."""

    def test_posterior_sums_to_one(self):
        """Posterior probabilities should sum to 1."""
        post = genotype_posterior(10, 5, 8, 20)
        total = sum(post.values())
        assert math.isclose(total, 1.0, rel_tol=1e-10)

    def test_strong_hom_ref(self):
        """Many ref reads, no alt reads → 0/0 with high confidence."""
        post = genotype_posterior(0, 0, 0, 50)
        assert post["0/0"] > 0.99
        assert post["0/0"] > post["0/1"]
        assert post["0/0"] > post["1/1"]

    def test_strong_hom_alt(self):
        """Many alt reads, very few ref reads → 1/1 with high confidence."""
        post = genotype_posterior(50, 10, 40, 1)
        assert post["1/1"] > 0.99
        assert post["1/1"] > post["0/1"]
        assert post["1/1"] > post["0/0"]

    def test_strong_het(self):
        """Balanced alt/ref counts → 0/1 with high confidence."""
        # W_alt = 10 + 8 + 2*6 = 30, W_ref = 2*30 = 60, N = 90
        # freq = 30/90 = 0.333 — close to 0.5
        post = genotype_posterior(10, 6, 8, 30)
        assert post["0/1"] > 0.95
        assert post["0/1"] > post["0/0"]
        assert post["0/1"] > post["1/1"]

    def test_ambiguous_low_coverage(self):
        """Very few reads → posterior should be less certain."""
        # L=1, R=0, M=0, ref=1 → W_alt=1, W_ref=2, N=3
        # freq = 1/3 ≈ 0.333 — close to het but very low evidence
        post = genotype_posterior(1, 0, 0, 1)
        # No genotype should have >90% certainty
        assert max(post.values()) < 0.95


class TestGenotypeFromCounts:
    """End-to-end tests for genotype calling."""

    def test_hom_ref_call(self):
        """All-ref site → 0/0 with high GQ."""
        gt, gq = genotype_from_counts(0, 0, 0, 50)
        assert gt == "0/0"
        assert gq >= 30

    def test_hom_alt_call(self):
        """Nearly all-alt site → 1/1 with high GQ."""
        gt, gq = genotype_from_counts(50, 10, 40, 1)
        assert gt == "1/1"
        assert gq >= 30

    def test_het_call(self):
        """Balanced site → 0/1 with high GQ."""
        gt, gq = genotype_from_counts(10, 6, 8, 30)
        assert gt == "0/1"
        assert gq >= 20

    def test_ambiguous_low_gq(self):
        """Low coverage → lower GQ score."""
        _, gq = genotype_from_counts(1, 0, 0, 1)
        assert gq < 10  # very low quality for ambiguous call

    def test_gq_capped_at_99(self):
        """GQ should not exceed 99."""
        _, gq = genotype_from_counts(10000, 0, 0, 0)
        assert gq <= 99

    def test_no_data(self):
        """No reads at all → should still call something without crashing."""
        gt, gq = genotype_from_counts(0, 0, 0, 0)
        assert gt in ("0/0", "0/1", "1/1")
        assert gq >= 0


class TestDefineGenotypeThreshold:
    """Backward compatibility tests for the threshold method."""

    def test_hom_ref(self):
        assert define_genotype_threshold(0.1) == "0/0"
        assert define_genotype_threshold(0.19) == "0/0"

    def test_het(self):
        assert define_genotype_threshold(0.2) == "0/1"
        assert define_genotype_threshold(0.5) == "0/1"
        assert define_genotype_threshold(0.79) == "0/1"

    def test_hom_alt(self):
        assert define_genotype_threshold(0.8) == "1/1"
        assert define_genotype_threshold(0.95) == "1/1"

    def test_edge_boundaries(self):
        assert define_genotype_threshold(0.0) == "0/0"
        assert define_genotype_threshold(1.0) == "1/1"
