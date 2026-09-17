from pathlib import Path
import numpy as np


class DataLoaderMixin:
    """Load observed swell curves, add synthetic noise to the curve, set sigma priors."""

    def load_data(self, filename, u_avg_obs):
        path = Path(self.swell_root, filename).resolve()
        self.infer_dir = folder = path.parent
        self.u_avg_obs, self.observed_uavgs = u_avg_obs, [u_avg_obs]

        y = np.atleast_2d(np.loadtxt(path))
        self._obs_indices = idx = np.arange(0, y.shape[1], self.thin)
        self.y_obs_matrix_clean = y = y[:, idx]
        self.obs_x_coords = x = np.loadtxt(folder / "curve4_x.txt").ravel()[idx]

        disp = y.max() - 1.0
        if not (np.isfinite(disp) and disp > 0):
            raise ValueError("no swelling in observed curve (max height <= die radius)")
        self.max_displacement = disp
        self.beta = self.sigma_noise_prior = rate = 1.0 / (0.10 * disp)
        self.sigma_bias_prior = rate if self._infer_sigma_bias() else None

        self.sigma_noise_target = self.sigma_noise_percent / 100 * disp
        noise = np.random.default_rng(0).normal(0, self.sigma_noise_target, y.shape)
        self.sigma_noise_realized = noise.std()
        self.y_obs_matrix = y + noise

        info = {
            "folder": folder,
            "curve file": path.name,
            "U_avg": u_avg_obs,
            "points": f"{idx.size} (thin={self.thin})",
            "x range": f"[{x.min():.4g}, {x.max():.4g}]",
            "height range": f"[{y.min():.4g}, {y.max():.4g}]",
            "max displacement": f"{disp:.4g}",
            "sigma_noise prior": f"Exp(rate={rate:.4g}), mean {1 / rate:.4g}",
            "sigma_bias prior": "same as sigma_noise" if self.sigma_bias_prior else "not inferred",
            "curve noise": f"{self.sigma_noise_percent}% of disp, sigma {self.sigma_noise_target:.4g} "
                           f"(realized {self.sigma_noise_realized:.4g})",
        }

        if self.use_pressure:
            self.pressure_obs = self.pressure_obs_clean = np.loadtxt(folder / "pressure_drop.txt").ravel()
            self.pressure_drop_obs = float(self.pressure_obs[0])
            info["pressure"] = f"{np.array2string(self.pressure_obs, precision=4)} (no noise)"

        width = max(map(len, info))
        print("\n".join(f"{k:<{width}} : {v}" for k, v in info.items()))