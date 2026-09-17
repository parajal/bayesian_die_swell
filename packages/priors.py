"""Prior distributions."""

import numpy as np
class PriorMixin:

    def log_prior(self, phi):
        """Uniform priors on physical parameters, exponential priors on sigmas."""
        theta = np.asarray(self._to_physical(phi), float)
        lo, hi = np.asarray(self._get_parameter_bounds(), float).T
        rates = np.array([self.sigma_noise_prior]
                         + [self.sigma_bias_prior] * self._infer_sigma_bias())
        params, sigmas = theta[:len(lo)], theta[len(lo):len(lo) + len(rates)]

        if (not np.isfinite(theta).all() or np.any((params < lo) | (params > hi))
                or np.any(sigmas <= 0)):
            return -np.inf
        return float(-np.log(hi - lo).sum() + np.sum(np.log(rates) - rates * sigmas))

    def _log_prior(self, rng, n_samples):
        """Draw n_samples from the prior."""
        lo, hi = np.asarray(self._get_parameter_bounds(), float).T
        rates = np.array([self.sigma_noise_prior]
                         + [self.sigma_bias_prior] * self._infer_sigma_bias())
        n, k = len(lo), len(rates)

        samples = np.zeros((n_samples, self._get_ndim()))
        samples[:, :n] = rng.uniform(lo, hi, size=(n_samples, n))
        samples[:, n:n + k] = rng.exponential(1 / rates, size=(n_samples, k))
        return samples