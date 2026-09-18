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
        if self._infer_sigma_bias() and self.sigma_bias_pc is not None:
            # PC elicitation (Fuglstad et al. 2019): P(sigma_bias > sigma0) = alpha,
            # sigma0 a fraction of the displacement. rate = -ln(alpha) / sigma0.
            frac, alpha = self.sigma_bias_pc
            self._sigma0 = frac * disp
            self.sigma_bias_prior = -np.log(alpha) / self._sigma0
        elif self._infer_sigma_bias():
            self.sigma_bias_prior = 1.0 / (self._sigma_bias_frac * disp)
        else:
            self.sigma_bias_prior = None
        self.c_bias_prior_sd = 0.10 * disp if self._infer_mean_bias() else None
        if self._infer_l_bias() and self.l_bias_prior == "pc":
            # PC prior (Fuglstad/Simpson) on the length scale: penalizes short (rough)
            # scales toward the smooth base model, set by P(l < l0) = alpha with l0 a
            # fraction of the x-range. 1-D range: lambda = -ln(alpha) * l0**0.5.
            frac, alpha = self.l_bias_pc
            self._l0 = frac * float(x.max() - x.min())
            self.l_bias_pc_lambda = -np.log(alpha) * self._l0 ** 0.5

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
            "sigma_bias prior": (
                f"PC: P(sigma>{self._sigma0:.4g})={self.sigma_bias_pc[1]:g}, rate={self.sigma_bias_prior:.4g}"
                if self.sigma_bias_prior and self.sigma_bias_pc is not None else
                f"Exp(mean={self._sigma_bias_frac * disp:.4g})" if self.sigma_bias_prior else "not inferred"),
            "c_bias prior": f"Normal(0, {self.c_bias_prior_sd:.4g})" if self.c_bias_prior_sd else "not inferred",
            "l_bias prior": (
                f"Uniform({self.l_bias_bounds[0]:g}, {self.l_bias_bounds[1]:g})"
                if self._infer_l_bias() and self.l_bias_prior == "uniform" else
                f"PC: P(l<{self._l0:.4g})={self.l_bias_pc[1]:g}, lambda={self.l_bias_pc_lambda:.4g}"
                if self._infer_l_bias() else f"fixed {self.l_bias:g}"),
            "bias kernel": self.correlation_matrix,
            "bias anchor": f"delta({self.bias_anchor:g})=0" if self.bias_anchor is not None else "none",
            "bias flat": (f"delta'=0 on [{self.bias_flat[0]:g},{self.bias_flat[1]:g}] ({self.bias_flat_n} pts)"
                          if self.bias_flat is not None else "none"),
            "curve noise": f"{self.sigma_noise_percent}% of disp, sigma {self.sigma_noise_target:.4g} "
                           f"(realized {self.sigma_noise_realized:.4g})",
        }

        if self.use_pressure:
            self.pressure_obs = self.pressure_obs_clean = np.loadtxt(folder / "pressure_drop.txt").ravel()
            self.pressure_drop_obs = float(self.pressure_obs[0])
            info["pressure"] = f"{np.array2string(self.pressure_obs, precision=4)} (no noise)"

        width = max(map(len, info))
        print("\n".join(f"{k:<{width}} : {v}" for k, v in info.items()))