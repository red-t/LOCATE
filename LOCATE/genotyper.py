"""Bayesian genotyping for TE insertions using a Beta-Binomial model.

Three genotype-specific Beta priors on the true allele frequency θ produce
Beta-Binomial marginal likelihoods for the observed read counts. The genotype
with the highest posterior probability is called, with a Phred-scale quality
score.

Posterior:  P(G | data) ∝ P(data | G) × P(G)
            P(G) = 1/3 (uniform prior on genotypes)
"""

import math

# Beta prior parameters for each genotype: Beta(alpha, beta)
# These define the expected distribution of allele frequency θ for each state.
#   0/0 (hom ref):  θ concentrated near 0   — rare alt reads from noise
#   0/1 (het):      θ concentrated near 0.5 — balanced allelic ratio
#   1/1 (hom alt):  θ concentrated near 1   — rare ref reads from noise
PRIOR_00 = (1, 30)   # Beta(1, 30)   → mean θ ≈ 0.032
PRIOR_01 = (5, 5)    # Beta(5, 5)    → mean θ = 0.5
PRIOR_11 = (30, 1)   # Beta(30, 1)   → mean θ ≈ 0.968

GENOTYPES = ("0/0", "0/1", "1/1")
_PRIOR_MAP = {"0/0": PRIOR_00, "0/1": PRIOR_01, "1/1": PRIOR_11}


def _beta_binomial_log_likelihood(k: int, n: int, alpha: float, beta: float) -> float:
    """Log marginal likelihood of k successes in n trials under Beta(α, β).

    P(k | n, α, β) = C(n,k) × B(α+k, β+n-k) / B(α, β)

    Computed in log-space with lgamma for numerical stability.
    """
    if n == 0:
        return 0.0
    # log(C(n,k)) = log(n!) - log(k!) - log((n-k)!)
    log_comb = math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)
    # log(B(α+k, β+n-k)) = lgamma(α+k) + lgamma(β+n-k) - lgamma(α+β+n)
    log_num = math.lgamma(alpha + k) + math.lgamma(beta + n - k) - math.lgamma(alpha + beta + n)
    # log(B(α, β)) = lgamma(α) + lgamma(β) - lgamma(α+β)
    log_den = math.lgamma(alpha) + math.lgamma(beta) - math.lgamma(alpha + beta)
    return log_comb + log_num - log_den


def genotype_posterior(
    leftclip: int, spanning: int, rightclip: int, numref: int,
) -> dict[str, float]:
    """Compute posterior probabilities for each genotype from raw read counts.

    Parameters
    ----------
    leftclip : int
        Number of left-clip supporting reads.
    spanning : int
        Number of spanning/mid-insert supporting reads.
    rightclip : int
        Number of right-clip supporting reads.
    numref : int
        Number of reference-spanning reads.

    Returns
    -------
    dict[str, float]
        Mapping of genotype -> posterior probability (summing to 1.0).
    """
    w_alt = leftclip + rightclip + 2 * spanning
    w_ref = 2 * numref
    n = w_alt + w_ref

    log_likelihoods = {}
    for gt in GENOTYPES:
        alpha, beta = _PRIOR_MAP[gt]
        log_likelihoods[gt] = _beta_binomial_log_likelihood(w_alt, n, alpha, beta)

    # log-likelihood → posterior (uniform prior P(G) = 1/3, so P(G|D) ∝ P(D|G))
    max_ll = max(log_likelihoods.values())
    unnorm = {gt: math.exp(ll - max_ll) for gt, ll in log_likelihoods.items()}
    total = sum(unnorm.values())
    return {gt: unnorm[gt] / total for gt in GENOTYPES}


def genotype_from_counts(
    leftclip: int, spanning: int, rightclip: int, numref: int,
) -> tuple[str, int]:
    """Call genotype using Bayesian inference from raw read counts.

    Parameters
    ----------
    leftclip : int
        Number of left-clip supporting reads.
    spanning : int
        Number of spanning/mid-insert supporting reads.
    rightclip : int
        Number of right-clip supporting reads.
    numref : int
        Number of reference-spanning reads.

    Returns
    -------
    tuple[str, int]
        (genotype_string, genotype_quality)
        Genotype quality is Phred-scale: -10 * log10(1 - max_posterior), capped at 99.
    """
    posteriors = genotype_posterior(leftclip, spanning, rightclip, numref)
    best_gt = max(posteriors, key=posteriors.get)
    max_prob = posteriors[best_gt]

    # Genotype quality: Phred-scale probability that the call is wrong
    gq = min(99, round(-10.0 * math.log10(max(1e-6, 1.0 - max_prob))))
    return best_gt, gq


def define_genotype_threshold(frequency: float) -> str:
    """Original hard-threshold genotyping for backward compatibility."""
    if frequency < 0.2:
        return "0/0"
    elif frequency >= 0.8:
        return "1/1"
    else:
        return "0/1"
