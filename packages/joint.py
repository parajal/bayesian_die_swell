"""Joint Bayesian inference over several U_avg datasets, each with its own ROM.

``JointInference`` extends the single-dataset ``ROMCurve4BayesianInference``:
one ROM is trained per dataset (per U_avg), one observed curve is loaded per
U_avg, and the joint log-likelihood is the sum of the per-curve Gaussian
log-likelihoods -- all sharing the material parameters and the noise/bias
parameters.
"""

from pathlib import Path
from typing import Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np

from .main import ROMCurve4BayesianInference
from .rom import ROM


class JointInference(ROMCurve4BayesianInference):
    """Infer shared material parameters from several curves at different U_avg."""

    def __init__(
        self,
        swell_root=None,
        train_data_rels: Sequence[str] = (),
        u_avgs: Sequence[float] = (),
        **kwargs,
    ) -> None:
        train_data_rels = list(train_data_rels)
        u_avgs = [float(u) for u in u_avgs]
        if len(train_data_rels) < 1:
            raise ValueError("Provide at least one dataset in train_data_rels.")
        if len(u_avgs) != len(train_data_rels):
            raise ValueError("u_avgs must have one value per training dataset.")

        super().__init__(swell_root=swell_root, train_data_rels=[train_data_rels[0]], **kwargs)
        if getattr(self, "use_pressure", False):
            raise NotImplementedError(
                "use_pressure is only supported for single-dataset ROMCurve4BayesianInference, "
                "not JointInference."
            )
        self.joint_data_rels = train_data_rels
        self.joint_uavgs = u_avgs
        self.roms: list[ROM] = []
        self.obs_curves: list[np.ndarray] = []
        self.obs_indices_list: list[np.ndarray] = []
        self.obs_x_list: list[np.ndarray] = []

    # ------------------------------------------------------------------ train
    def build_rom(self) -> None:
        """Train one ROM per dataset and adopt the shared material metadata."""
        self.roms = []
        for data_rel in self.joint_data_rels:
            rom = ROM(
                data_dir=(self.swell_root / data_rel).resolve(),
                filenames_train=self.filenames_train,
                filenames_test=self.filenames_test,
                method=self.method,
                scaler=self.scaler,
                eps=self.eps,
                model_family=self.model_family,
            )
            rom.lam_bounds = self.lam_bounds
            rom.beta_bounds = self.beta_bounds
            rom.third_parameter_bounds = self.third_parameter_bounds
            rom.alpha_bounds = self.third_parameter_bounds
            rom.epsilon_bounds = getattr(self, "epsilon_bounds", None)
            rom.train()
            self.roms.append(rom)

        r0 = self.roms[0]  # material grid is shared across datasets
        self.param_cols = r0.param_cols
        self.material_parameter_names = list(r0.material_parameter_names)
        self.third_parameter_name = r0.third_parameter_name
        self.lam_idx, self.beta_idx, self.alpha_idx = r0.lam_idx, r0.beta_idx, r0.alpha_idx
        self.lam_bounds, self.beta_bounds = r0.lam_bounds, r0.beta_bounds
        self.third_parameter_bounds = r0.third_parameter_bounds
        self.alpha_bounds = r0.alpha_bounds
        self.epsilon_bounds = r0.epsilon_bounds
        self.use_uavg = r0.use_uavg
        self.is_trained = True

    train = build_rom

    def validate_rom(self) -> None:
        total_n = 0
        weighted_error = 0.0
        for uavg, rom in zip(self.joint_uavgs, self.roms):
            print(f"\n=== ROM validation, U_avg={uavg:g} ===")
            rom.validate_rom()
            stats = rom.rom_test_error_summary()
            n_cases = int(stats["num_test_cases"])
            total_n += n_cases
            weighted_error += n_cases * float(stats["relative_l2_error_mean"])

        if total_n <= 0:
            raise RuntimeError("No ROM test cases were validated.")
        self._rom_val_overall_n = total_n
        self._rom_val_overall_err = weighted_error / total_n
        print("\n=== Joint ROM validation summary ===")
        print(
            f"{'Overall':>10}  {self._rom_val_overall_n:>6d}  "
            f"{self._rom_val_overall_err:>18.4e}"
        )

    # ------------------------------------------------------------------- data
    def load_data(self, filenames: Sequence[str], u_avgs: Optional[Sequence[float]] = None) -> None:
        """Load one observed curve per dataset, add noise, and set shared priors."""
        filenames = list(filenames)
        if u_avgs is not None:
            self.joint_uavgs = [float(u) for u in u_avgs]
        if len(filenames) != len(self.joint_data_rels):
            raise ValueError("Provide one observation file per dataset.")

        self.obs_curves, self.obs_indices_list, self.obs_x_list = [], [], []
        self.obs_stats = []
        clean_rows = []
        for fn in filenames:
            raw = Path(fn).expanduser()
            path = (raw if raw.is_absolute() else self.swell_root / raw).resolve()
            clean = np.atleast_2d(np.loadtxt(path, dtype=float))
            if clean.shape[0] != 1:
                raise ValueError(f"{path.name} has {clean.shape[0]} rows; provide a single-curve file.")
            idx = np.arange(0, clean.shape[1], self.thin, dtype=int)
            if idx.size < 2:
                raise ValueError("thin leaves fewer than two observation points.")
            clean = clean[:, idx]
            noisy = self._add_noise(clean, seed=0)[0]
            x_grid = self._load_x_grid(path.parent, idx, noisy.size)
            self.obs_curves.append(noisy)
            self.obs_indices_list.append(idx)
            self.obs_x_list.append(x_grid)
            clean_rows.append(clean.ravel())

            # Curve diagnostics: raw peak height, swelling rise, and end height at x = 5.
            y = clean.ravel()
            i_max = int(np.argmax(y))
            i_end = int(np.argmin(np.abs(x_grid - 5.0)))
            self.obs_stats.append(
                {
                    "name": path.name,
                    "max_height": float(y[i_max]),
                    "max_swell_disp": float(y[i_max] - 1.0),
                    "x_max": float(x_grid[i_max]),
                    "end_height": float(y[i_end]),
                }
            )

        # Shared sigma priors from the combined maximum swelling displacement.
        self._sigma_priors_from_data(np.concatenate(clean_rows)[None, :])

        # Single-curve attributes kept for compatibility with inherited helpers.
        self.y_obs_matrix = self.obs_curves[0][None, :]
        self.obs_x_coords = self.obs_x_list[0]
        self._obs_indices = self.obs_indices_list[0]
        self.u_avg_obs = self.joint_uavgs[0]
        self.observed_uavgs = list(self.joint_uavgs)
        self._obs_filename = filenames[0]
        print(
            f"Joint data: {len(self.obs_curves)} curves at U_avg={self.joint_uavgs}, "
            f"n_points={[int(c.size) for c in self.obs_curves]}, "
            f"sigma_noise prior Exp(beta={self.beta:.4g})"
        )
        for u, s in zip(self.joint_uavgs, self.obs_stats):
            print(
                f"  U_avg={u}: max h={s['max_height']:.6e} "
                f"(swell displacement={s['max_swell_disp']:.6e}, x={s['x_max']:.4f}), "
                f"h(x=5)={s['end_height']:.6e}"
            )

    @staticmethod
    def _load_x_grid(obs_dir: Path, idx: np.ndarray, n_pts: int) -> np.ndarray:
        x_path = obs_dir / "curve4_x.txt"
        if not x_path.exists():
            raise FileNotFoundError(f"Missing x-axis file: {x_path}")

        x = np.loadtxt(x_path, dtype=float).ravel()
        if x.size >= int(idx[-1]) + 1:
            x = x[idx]
        if x.size != n_pts or not np.all(np.isfinite(x)):
            raise ValueError("curve4_x.txt must be finite and match the thinned curve4_y.txt length.")
        return x

    # -------------------------------------------------------------- likelihood
    def _predict_for(self, rom: ROM, lam: float, bet: float, alp: Optional[float], uavg: float,
                     idx: np.ndarray, n_obs: int) -> np.ndarray:
        y = np.asarray(
            rom.predict(lam, bet, alpha_val=alp, u_avg_val=uavg if rom.use_uavg else None),
            dtype=float,
        ).ravel()
        if y.size != n_obs and y.size >= int(idx[-1]) + 1:
            y = y[idx]
        return y

    def log_likelihood(self, phi: np.ndarray) -> float:
        """Sum of the per-curve Gaussian log-likelihoods (shared parameters)."""
        phi = np.asarray(phi, dtype=float)
        if phi.shape != (self._get_ndim(),) or not np.all(np.isfinite(phi)):
            return -np.inf
        try:
            theta = self._to_physical(phi)
            sigma_noise, sigma_bias = self._extract_noise_bias(theta)
        except (FloatingPointError, RuntimeError, ValueError):
            return -np.inf
        if not (sigma_noise > 0.0 and np.isfinite(sigma_noise)):
            return -np.inf
        if sigma_bias is not None and not (sigma_bias > 0.0 and np.isfinite(sigma_bias)):
            return -np.inf

        n_material = len(self._get_parameter_bounds())
        lam, bet = float(theta[0]), float(theta[1])
        alp = float(theta[2]) if n_material > 2 else None
        mode = getattr(self, "mode", "full_curve")

        total = 0.0
        for rom, obs, idx, x, uavg in zip(
            self.roms, self.obs_curves, self.obs_indices_list, self.obs_x_list, self.joint_uavgs
        ):
            try:
                y_pred = self._predict_for(rom, lam, bet, alp, uavg, idx, obs.size)
            except (FloatingPointError, ValueError):
                return -np.inf
            if y_pred.shape != obs.shape or not np.all(np.isfinite(y_pred)):
                return -np.inf

            if mode == "swell_height":
                residual = np.asarray(
                    [self._swell_height(obs) - self._swell_height(y_pred)], dtype=float
                )
                x_res = np.asarray([x[int(np.argmax(obs))]], dtype=float)
            else:
                residual = obs - y_pred
                x_res = x

            if sigma_bias is None or sigma_bias == 0.0:
                ll = self._log_likelihood(residual, sigma_noise=sigma_noise)
            else:
                cov = self._covariance_matrix(x_res, sigma_noise, float(sigma_bias), self.l_bias)
                cholesky = self._cholesky(cov)
                if cholesky is None:
                    return -np.inf
                ll = self._log_likelihood(residual, cholesky=cholesky)
            if not np.isfinite(ll):
                return -np.inf
            total += ll
        return float(total)

    # ---------------------------------------------------------------- plotting
    def plot_data(self) -> None:
        fig, ax = self._curve_figure()
        for i, (obs, x, u) in enumerate(zip(self.obs_curves, self.obs_x_list, self.joint_uavgs)):
            ax.scatter(x, obs, s=16, alpha=0.8, color=f"C{i}", label=f"U_avg = {u:g}")
        self._set_curve_labels(ax)
        self._axis_style(ax, np.concatenate([np.ravel(o) for o in self.obs_curves]))
        self._curve_grid(ax)
        self._curve_legend(ax)
        self._save_current_figure("data")
        plt.show()

    def plot_rom_prediction(self) -> None:
        theta = getattr(self, "true_theta", None)
        if theta is None:
            raise ValueError("Set true_theta to use plot_rom_prediction().")
        lam, bet = float(theta[0]), float(theta[1])
        alp = float(theta[2]) if len(theta) > 2 and self.alpha_idx >= 0 else None

        fig, ax = self._curve_figure()
        for i, (rom, obs, idx, x, u) in enumerate(
            zip(self.roms, self.obs_curves, self.obs_indices_list, self.obs_x_list, self.joint_uavgs)
        ):
            y_rom = self._predict_for(rom, lam, bet, alp, u, idx, obs.size)
            ax.scatter(x, obs, s=24, alpha=0.9, color=f"C{i}", label=f"FOM+noise U={u:g}")
            ax.plot(x, y_rom, color=f"C{i}", lw=2.0, label=f"ROM U={u:g}")
        self._set_curve_labels(ax)
        self._axis_style(ax, np.concatenate([np.ravel(o) for o in self.obs_curves]))
        self._curve_grid(ax)
        self._curve_legend(ax)
        self._save_current_figure("rom_prediction")
        plt.show()

    def _predictive_band(self, rom, obs, idx, x, uavg, n_sigma, nsamples_pred, condition_discrepancy):
        """Posterior-predictive band and coverage diagnostics for one U_avg curve."""
        from scipy.stats import norm

        samples = self._posterior_samples()
        n_material = len(self._get_parameter_bounds())
        i_sn = n_material
        infer_bias = self._infer_sigma_bias()

        rng = np.random.default_rng(0)
        n_draws = min(int(nsamples_pred), samples.shape[0])
        sel = rng.choice(samples.shape[0], size=n_draws, replace=False)
        draws = samples[sel]
        sn_draws = np.abs(draws[:, i_sn])
        sb_draws = np.abs(draws[:, i_sn + 1]) if infer_bias else np.zeros(n_draws)
        l_bias = float(getattr(self, "l_bias", 1.0) or 1.0)
        corr = self._bias_correlation_matrix(x, l_bias) if infer_bias else None

        latent_mean = np.empty((n_draws, obs.size), dtype=float)
        Y_rep = np.empty((n_draws, obs.size), dtype=float)
        for k, s in enumerate(draws):
            alp = float(s[2]) if n_material > 2 else None
            g = self._predict_for(rom, float(s[0]), float(s[1]), alp, uavg, idx, obs.size)
            sn, sb = float(sn_draws[k]), float(sb_draws[k])
            if not infer_bias or sb <= 0.0:
                delta = delta_mean = np.zeros_like(g)
            elif condition_discrepancy:
                A = sb * sb * corr
                Sigma = 0.5 * (A + A.T) + sn * sn * np.eye(len(g))
                Sinv = np.linalg.pinv(Sigma, hermitian=True)
                delta_mean = A @ (Sinv @ (obs - g))
                delta = delta_mean + self._correlated_normal(rng, A - A @ (Sinv @ A), 1.0)
            else:
                delta = self._correlated_normal(rng, corr, sb)
                delta_mean = np.zeros_like(g)
            latent_mean[k] = g + delta_mean
            Y_rep[k] = g + delta + rng.standard_normal(len(g)) * sn

        tail = 100.0 * float(norm.cdf(n_sigma))
        mean_pred = latent_mean.mean(0)
        pred_lo = np.percentile(Y_rep, 100.0 - tail, axis=0)
        pred_hi = np.percentile(Y_rep, tail, axis=0)
        sigma_total = np.maximum(np.std(Y_rep, axis=0), 1e-8)
        zres = (obs - mean_pred) / sigma_total
        diag = dict(
            coverage=float(np.mean((obs >= pred_lo) & (obs <= pred_hi))),
            rms_z=float(np.sqrt(np.mean(zres**2))),
            max_abs_z=float(np.max(np.abs(zres))),
            mean_z=float(np.mean(zres)),
            nominal_coverage=float(2.0 * norm.cdf(n_sigma) - 1.0),
        )
        return mean_pred, pred_lo, pred_hi, diag

    def plot_posterior_predictive(
        self,
        n_sigma: float = 1.96,
        nsamples_pred: int = 5000,
        condition_discrepancy: bool = False,
    ) -> None:
        """Posterior predictive band per U_avg curve, in the reference style."""
        from matplotlib.lines import Line2D
        from scipy.stats import norm

        self.pp_diagnostics = {}
        pct = int(round(100 * (2.0 * norm.cdf(n_sigma) - 1.0)))
        for rom, obs, idx, x, u in zip(
            self.roms, self.obs_curves, self.obs_indices_list, self.obs_x_list, self.joint_uavgs
        ):
            mean_pred, pred_lo, pred_hi, diag = self._predictive_band(
                rom, obs, idx, x, u, n_sigma, nsamples_pred, condition_discrepancy
            )
            self.pp_diagnostics[u] = diag

            fig, ax = self._curve_figure()
            ax.fill_between(x, pred_lo, pred_hi, color="steelblue", alpha=0.25)
            ax.plot(x, mean_pred, color="steelblue", lw=1.5, zorder=4)
            ax.scatter(x, obs, color="black", s=12, zorder=5, alpha=0.8,
                       edgecolors="black", linewidths=0.5, label="_nolegend_")
            self._set_curve_labels(ax)
            self._axis_style(ax, np.concatenate([obs, pred_hi]))
            self._curve_grid(ax)
            band = Line2D([0], [0], color="steelblue", lw=6, alpha=0.25, label="Posterior predictive")
            h, lbl = ax.get_legend_handles_labels()
            self._curve_legend(ax, h + [band], lbl + ["Posterior predictive"])
            print(f"\nPosterior predictive (U={u:g}, {pct}% band)")
            print(f"  coverage={diag['coverage']:.1%}  rms_z={diag['rms_z']:.2f}  "
                  f"max|z|={diag['max_abs_z']:.2f}  mean_z={diag['mean_z']:.2f}")
            self._save_current_figure(f"posterior_predictive_U{u:g}")
            plt.show()
