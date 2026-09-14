from pathlib import Path
import numpy as np

def _max_disp(y):
    return np.max(y) - 1.0

class DataLoaderMixin:

    def _add_noise(self, y, seed=0):
        sigma = self.sigma_noise_percent / 100 * _max_disp(y)
        noise = np.random.default_rng(seed).normal(0, sigma, y.shape)
        self.max_displacement = _max_disp(y)
        self.sigma_noise_target = sigma
        self.sigma_noise_realized = noise.std()
        return y + noise

    def _sigma_priors_from_data(self, y):
        self.beta = 1 / (0.10 * _max_disp(y))
        self.sigma_noise_prior = self.beta
        self.sigma_bias_prior = self.beta if self._infer_sigma_bias() else None

    def load_data(self, filename, u_avg_obs, pressure_filename=None):

        path = Path(self.swell_root, filename).resolve()

        y = np.atleast_2d(np.loadtxt(path))
        self._obs_indices = np.arange(0, y.shape[1], self.thin)

        self.u_avg_obs = u_avg_obs
        self.observed_uavgs = [u_avg_obs]
        self.infer_dir = path.parent

        self.y_obs_matrix_clean = y[:, self._obs_indices]
        self._sigma_priors_from_data(self.y_obs_matrix_clean)
        self.y_obs_matrix = self._add_noise(self.y_obs_matrix_clean)

        self.obs_x_coords = np.loadtxt(
            path.parent / "curve4_x.txt"
        ).ravel()[self._obs_indices]

        if getattr(self, "use_pressure", False):
            self._load_pressure_obs(pressure_filename, path.parent)

    def _load_pressure_obs(self, filename, folder):

        if filename is None:
            path = Path(folder) / "pressure.txt"
        else:
            raw = Path(filename).expanduser()
            path = raw if raw.is_absolute() else Path(self.swell_root) / raw

        p = np.loadtxt(path).ravel()

        sigma = self.sigma_noise_target
        noise = np.random.default_rng(1).normal(0, sigma, p.shape)

        self.pressure_obs_clean = p
        self.pressure_obs = p + noise