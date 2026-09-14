"""Gaussian likelihood and posterior."""

import numpy as np
from scipy.linalg import cho_factor, cho_solve


class LikelihoodMixin:

    @staticmethod
    def _log_likelihood(r, sigma=None, chol=None):
        r = np.asarray(r, float).ravel()
        n = r.size

        if chol is None:
            r2 = (r @ r) / sigma**2
            logdet = 2 * n * np.log(sigma)
        else:
            r2 = r @ cho_solve(chol, r)
            logdet = 2 * np.sum(np.log(np.diag(chol[0])))

        return -0.5 * (r2 + logdet + n * np.log(2 * np.pi))

    @staticmethod
    def _swell_height(y):
        return float(np.max(y))

    def _covariance(self, x, sigma_noise, sigma_bias):
        K = np.exp(-np.abs(np.subtract.outer(x, x)) / self.l_bias)
        return sigma_noise**2 * np.eye(len(x)) + sigma_bias**2 * K

    def _observed_curve(self):
        return np.asarray(self.y_obs_matrix[0], float)

    def _rom_curve(self, theta, n):
        kwargs = {}

        if len(theta) > 2:
            name = getattr(self, "third_parameter_name", "")
            kwargs["epsilon_val" if name == "epsilon" else "alpha_val"] = theta[2]

        y = np.asarray(
            self.rom_predict_curve(
                theta[0], theta[1],
                u_avg_val=self.u_avg_obs,
                **kwargs,
            ),
            float,
        ).ravel()

        if hasattr(self, "_obs_indices") and y.size > n:
            y = y[self._obs_indices]

        return y

    def _residual(self, theta):

        y_obs = self._observed_curve()

        if getattr(self, "model_family", None) == "tanner":
            pred = self._tanner_height(theta[0], self.u_avg_obs)
            return np.array([self._swell_height(y_obs) - pred])

        y_pred = self._rom_curve(theta, y_obs.size)

        if getattr(self, "mode", "full_curve") == "swell_height":
            return np.array([
                self._swell_height(y_obs) -
                self._swell_height(y_pred)
            ])

        return y_obs - y_pred

    def log_likelihood(self, phi):

        theta = self._to_physical(phi)
        sigma_noise, sigma_bias = self._extract_noise_bias(theta)

        r = self._residual(theta)

        if sigma_bias is None:
            ll = self._log_likelihood(r, sigma=sigma_noise)
        else:
            x = self.obs_x_coords[:len(r)]
            C = self._covariance(x, sigma_noise, sigma_bias)
            ll = self._log_likelihood(r, chol=cho_factor(C, lower=True))

        if getattr(self, "use_pressure", False):
            rp = self.pressure_obs - self.predict_pressure(theta)
            ll += self._log_likelihood(rp, sigma=sigma_noise)

        return ll

    def log_posterior(self, phi):
        lp = self.log_prior(phi)
        return lp + self.log_likelihood(phi) if np.isfinite(lp) else -np.inf