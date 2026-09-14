"""Prior distributions."""
import numpy as np
class PriorMixin:

    def log_prior(self, phi):

        theta = self._to_physical(phi)
        bounds = np.asarray(self._get_parameter_bounds())

        vals = theta[: len(bounds)]
        lo, hi = bounds[:, 0], bounds[:, 1]

        if np.any((vals < lo) | (vals > hi)):
            return -np.inf

        logp = -np.sum(np.log(hi - lo))

        i = len(bounds)
        s = theta[i]

        if s <= 0:
            return -np.inf

        logp += np.log(self.sigma_noise_prior) - self.sigma_noise_prior * s

        if self._infer_sigma_bias():
            b = theta[i + 1]
            if b <= 0:
                return -np.inf
            logp += np.log(self.sigma_bias_prior) - self.sigma_bias_prior * b

        return float(logp)

    def _log_prior(self, rng, n):

        bounds = self._get_parameter_bounds()
        phi = np.zeros((n, self._get_ndim()))

        for j, (lo, hi) in enumerate(bounds):
            phi[:, j] = rng.uniform(lo, hi, n)

        i = len(bounds)

        phi[:, i] = rng.exponential(
            1 / self.sigma_noise_prior, n
        )

        if self._infer_sigma_bias():
            phi[:, i + 1] = rng.exponential(
                1 / self.sigma_bias_prior, n
            )

        return phi
