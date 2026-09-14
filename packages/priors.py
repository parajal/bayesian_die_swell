"""Prior distributions."""

import numpy as np


class PriorMixin:

    def log_prior(self, phi):
        """Evaluate log-prior density."""

        theta = self._to_physical(phi)

        bounds = np.asarray(self._get_parameter_bounds())
        n_params = len(bounds)

        params = theta[:n_params]
        lower = bounds[:, 0]
        upper = bounds[:, 1]

        # Uniform priors
        if np.any((params < lower) | (params > upper)):
            return -np.inf

        logp = -np.sum(np.log(upper - lower))

        # Noise standard deviation prior
        sigma_noise = theta[n_params]
        if sigma_noise <= 0:
            return -np.inf

        logp += (
            np.log(self.sigma_noise_prior)
            - self.sigma_noise_prior * sigma_noise
        )

        # Bias standard deviation prior
        if self._infer_sigma_bias():
            sigma_bias = theta[n_params + 1]

            if sigma_bias <= 0:
                return -np.inf

            logp += (
                np.log(self.sigma_bias_prior)
                - self.sigma_bias_prior * sigma_bias
            )

        return float(logp)

    def _log_prior(self, rng, n_samples):
        """Draw samples from the prior."""

        bounds = self._get_parameter_bounds()
        ndim = self._get_ndim()

        samples = np.zeros((n_samples, ndim))

        # Uniform parameters
        for j, (low, high) in enumerate(bounds):
            samples[:, j] = rng.uniform(low, high, size=n_samples)

        idx = len(bounds)

        # Exponential prior for noise
        samples[:, idx] = rng.exponential(
            scale=1.0 / self.sigma_noise_prior,
            size=n_samples,
        )

        # Optional bias parameter
        if self._infer_sigma_bias():
            samples[:, idx + 1] = rng.exponential(
                scale=1.0 / self.sigma_bias_prior,
                size=n_samples,
            )

        return samples
