"""Minimal plotting helpers for ROM Bayesian inference."""

from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 30,
    "axes.labelsize": 30,
    "xtick.labelsize": 30,
    "ytick.labelsize": 30,
    "legend.fontsize": 20,
    "lines.linewidth": 2,
    "axes.linewidth": 1,
})


class PlottingMixin:
    """Plots use physical coordinates."""

    CURVE_FIGSIZE = (8, 6)
    CURVE_GRID_ALPHA = 0.3
    CURVE_LEGEND_LOC = "lower right"

    def _curve_figure(self):
        return plt.subplots(figsize=self.CURVE_FIGSIZE)

    @staticmethod
    def _set_curve_labels(ax) -> None:
        ax.set(xlabel=r"$x$", ylabel=r"$h(x)$")

    def _curve_grid(self, ax) -> None:
        ax.grid(True, alpha=self.CURVE_GRID_ALPHA)

    def _curve_legend(self, ax, handles=None, labels=None, loc: Optional[str] = None) -> None:
        loc = self.CURVE_LEGEND_LOC if loc is None else loc
        if handles is None:
            ax.legend(loc=loc, framealpha=0.9)
            return
        ax.legend(handles, labels, loc=loc, framealpha=0.9)

    def _labels(self, latex: bool = True) -> "list[str]":
        labels = list(self._get_parameter_labels(latex=latex))
        samples = getattr(self, "samples", None)
        ndim = np.asarray(samples).shape[-1] if samples is not None else len(labels)
        labels.extend(f"param_{i}" for i in range(len(labels), ndim))
        return labels[:ndim]

    def _observation_x(self, n_pts: int, x_filename: str = "curve4_x.txt") -> np.ndarray:
        obs_x = getattr(self, "obs_x_coords", None)
        if obs_x is None:
            raise RuntimeError("obs_x_coords is unset. Call load_data() with curve4_x.txt present.")
        obs_x = np.asarray(obs_x, dtype=float).ravel()
        if obs_x.shape != (n_pts,) or not np.all(np.isfinite(obs_x)):
            raise ValueError("obs_x_coords must be a finite 1D array matching observed curves.")
        return obs_x

    def _posterior_samples(self) -> np.ndarray:
        if getattr(self, "samples", None) is None:
            raise RuntimeError("Call run_mcmc() first.")
        return np.asarray(self.samples, dtype=float)

    def _plot_output_dir(self) -> Path:
        obs_file = getattr(self, "_obs_filename", None)
        if obs_file:
            folder = f"{Path(obs_file).stem}_{self._noise_label()}"
            out = self.swell_root / "plots" / folder
        else:
            out = Path(self.infer_dir)
        out.mkdir(parents=True, exist_ok=True)
        return out

    def _noise_label(self) -> str:
        value = getattr(self, "sigma_noise_percent", None)
        if value is None:
            return "noise_unknown"
        try:
            text = f"{float(value):g}"
        except (TypeError, ValueError):
            text = str(value)
        return "noise_" + text.replace(" ", "_").replace("/", "_").replace("\\", "_")

    def _save_current_figure(self, filename: str) -> None:
        noise_label = self._noise_label()
        stem = filename if filename.endswith(f"_{noise_label}") else f"{filename}_{noise_label}"
        fig = plt.gcf()
        fig.tight_layout()
        fig.savefig(self._plot_output_dir() / f"{stem}.pdf", dpi=300)

    def _axis_style(self, ax, y_value) -> None:
        """Uniform x/y limits and 3 ticks for h(x) curve plots."""
        ymax = float(np.max(y_value))
        ytop = ymax + 0.01 * ymax
        ax.set(
            xlim=(0, 5),
            xticks=np.linspace(0, 5, 3),
            ylim=(1, ytop),
            yticks=np.linspace(1, ytop, 3),
        )

    def plot_data(self, x_filename: str = "curve4_x.txt") -> None:
        if getattr(self, "y_obs_matrix", None) is None:
            raise RuntimeError("Call load_data() first.")
        y = np.atleast_2d(self.y_obs_matrix)
        uavgs = None if self.observed_uavgs is None else list(self.observed_uavgs)
        if uavgs is not None and len(uavgs) != y.shape[0]:
            raise ValueError(
                f"Loaded {y.shape[0]} observed curves but only {len(uavgs)} U_avg values. "
                "Use a single-row observation file for one U_avg."
            )
        x = self._observation_x(y.shape[1], x_filename)
        clean = getattr(self, "y_obs_matrix_clean", None)
        clean = None if clean is None else np.atleast_2d(clean)
        fig, ax = self._curve_figure()
        for i, row in enumerate(y):
            u = uavgs[i] if uavgs is not None else self.u_avg_obs
            color = f"C{i}"
            if clean is not None and i < clean.shape[0]:
                ax.plot(x, clean[i], color=color, lw=2, alpha=0.9, label="_nolegend_")
            ax.scatter(x, row, s=16, alpha=0.8, color=color, label=f"U_avg = {u:g}")
        self._set_curve_labels(ax)
        self._axis_style(ax, y)
        self._curve_grid(ax)
        self._curve_legend(ax)
        self._save_current_figure("data")
        plt.show()

    def plot_rom_prediction(
        self,
        lambda_val: Optional[float] = None,
        beta_val: Optional[float] = None,
        tanner_ratio_val: Optional[float] = None,
        n1_val: Optional[float] = None,
        alpha_val: Optional[float] = None,
        epsilon_val: Optional[float] = None,
        u_avg_val: Optional[float] = None,
        x_filename: str = "curve4_x.txt",
    ) -> None:
        """Compare ROM prediction with the loaded noisy FOM curve4."""
        if not getattr(self, "is_trained", False):
            raise RuntimeError("Call train() before plot_rom_prediction().")
        if getattr(self, "y_obs_matrix", None) is None:
            raise RuntimeError("Call load_data() before plot_rom_prediction().")

        if getattr(self, "model_family", None) == "tanner":
            theta = getattr(self, "true_theta", None)
            u_avg = float(u_avg_val if u_avg_val is not None else self.u_avg_obs)
            tau_w = 4*u_avg
            if n1_val is None:
                if tanner_ratio_val is not None:
                    n1_val = 2.0 * tau_w * float(tanner_ratio_val)
                elif theta is not None:
                    n1_val = float(theta[0])
                elif getattr(self, "map_theta", None) is not None:
                    n1_val = float(self.map_theta[0])
                else:
                    raise RuntimeError(
                        "Tanner plotting needs n1_val, tanner_ratio_val, true_theta, "
                        "or posterior MAP samples."
                    )
            y_fom = np.atleast_2d(self.y_obs_matrix)[0]
            x = self._observation_x(len(y_fom), x_filename)
            n1 = float(n1_val)
            ratio = n1 / (2.0 * tau_w)
            pred_height = self._tanner_height(n1, u_avg_val=u_avg)
            obs_height = float(np.max(y_fom))
            i_max = int(np.argmax(y_fom))
            fig, ax = self._curve_figure()
            ax.plot(x, y_fom, color="tab:blue", lw=2, label="FOM + noise")
            ax.scatter([x[i_max]], [obs_height], s=36, color="black", zorder=5, label="observed max")
            ax.axhline(pred_height, color="tab:red", lw=2, label="Tanner")
            self._set_curve_labels(ax)
            self._axis_style(ax, np.append(y_fom, pred_height))
            self._curve_grid(ax)
            self._curve_legend(ax)
            self._save_current_figure("tanner_prediction")
            plt.show()
            print(f"N1: {n1:.6e}")
            print(f"tau_w from U_avg={u_avg:g}: {tau_w:.6e}")
            print(f"N1_over_2tau_w: {ratio:.6e}")
            print(f"Tanner predicted swell height: {pred_height:.6e}")
            print(f"Observed max curve4: {obs_height:.6e}")
            print(f"Residual: {obs_height - pred_height:.6e}")
            return

        if lambda_val is None or beta_val is None:
            theta = getattr(self, "true_theta", None)
            if theta is None:
                raise ValueError("Pass lambda_val and beta_val, or set true_theta on the model.")
            lambda_val, beta_val = float(theta[0]), float(theta[1])
            if alpha_val is None and epsilon_val is None and len(theta) > 2:
                if getattr(self, "third_parameter_name", None) == "epsilon":
                    epsilon_val = float(theta[2])
                else:
                    alpha_val = float(theta[2])

        u_avg = float(u_avg_val if u_avg_val is not None else self.u_avg_obs)
        y_fom = np.atleast_2d(self.y_obs_matrix)[0]
        y_rom = self.rom_predict_curve(
            float(lambda_val),
            float(beta_val),
            u_avg_val=u_avg,
            alpha_val=(alpha_val if self.alpha_idx >= 0 else None),
            epsilon_val=(epsilon_val if self.alpha_idx >= 0 else None),
        )
        idx = getattr(self, "_obs_indices", None)
        if idx is not None and y_rom.shape != y_fom.shape:
            idx = np.asarray(idx, dtype=int)
            if y_rom.size >= int(idx[-1]) + 1:
                y_rom = y_rom[idx]
        x = self._observation_x(len(y_fom), x_filename)
        mean_rel_err = np.mean(np.abs(y_rom - y_fom) / (np.abs(y_fom) + 1e-14))
        rel_l2_curve = float(np.linalg.norm(y_rom - y_fom) / (np.linalg.norm(y_fom) + 1e-14))
        fig, ax = self._curve_figure()
        ax.scatter(x, y_fom, s=24, alpha=0.9, color="tab:blue", label="FOM + noise")
        ax.scatter(x, y_rom, s=24, alpha=0.9, color="tab:red", label="ROM")
        self._set_curve_labels(ax)
        self._axis_style(ax, np.concatenate([y_fom, y_rom]))
        self._curve_grid(ax)
        self._curve_legend(ax)
        self._save_current_figure("rom_prediction")
        plt.show()
        print(f"ROM prediction mean relative error against noisy curve4: {mean_rel_err:.6e}")
        print(f"ROM prediction relative L2 error against noisy curve4:    {rel_l2_curve:.6e}")
        if getattr(self, "use_pressure", False) and getattr(self, "pressure_obs", None) is not None:
            self.plot_pressure_prediction(
                lambda_val=float(lambda_val),
                beta_val=float(beta_val),
                alpha_val=alpha_val,
                epsilon_val=epsilon_val,
            )

    def plot_pressure_prediction(
        self,
        lambda_val: Optional[float] = None,
        beta_val: Optional[float] = None,
        alpha_val: Optional[float] = None,
        epsilon_val: Optional[float] = None,
    ) -> float:
        """Compare the pressure GPR prediction with the observed pressure.

        Plots observed vs predicted pressure at each pressure point and prints
        (and returns) the relative L2 error ``||p_rom - p_obs|| / ||p_obs||``.
        """
        if not getattr(self, "use_pressure", False):
            raise RuntimeError("plot_pressure_prediction() requires use_pressure=True.")
        if getattr(self, "pressure_model", None) is None:
            raise RuntimeError("Call build_rom() before plot_pressure_prediction().")
        p_obs = getattr(self, "pressure_obs", None)
        if p_obs is None:
            raise RuntimeError(
                "Call load_data(..., pressure_filename=...) before plot_pressure_prediction()."
            )
        p_obs = np.asarray(p_obs, dtype=float).ravel()

        if lambda_val is None or beta_val is None:
            theta = getattr(self, "true_theta", None)
            if theta is None:
                raise ValueError("Pass lambda_val and beta_val, or set true_theta on the model.")
            lambda_val, beta_val = float(theta[0]), float(theta[1])
            if alpha_val is None and epsilon_val is None and len(theta) > 2:
                if getattr(self, "third_parameter_name", None) == "epsilon":
                    epsilon_val = float(theta[2])
                else:
                    alpha_val = float(theta[2])

        theta = [float(lambda_val), float(beta_val)]
        if self.alpha_idx >= 0:
            third = epsilon_val if epsilon_val is not None else alpha_val
            if third is None:
                raise ValueError("The third material parameter is required for this model.")
            theta.append(float(third))
        p_rom = np.asarray(self.predict_pressure(theta), dtype=float).ravel()

        rel_l2 = float(np.linalg.norm(p_rom - p_obs) / (np.linalg.norm(p_obs) + 1e-14))
        mean_rel = float(np.mean(np.abs(p_rom - p_obs) / (np.abs(p_obs) + 1e-14)))
        idx = np.arange(1, p_obs.size + 1)
        fig, ax = self._curve_figure()
        ax.plot(idx, p_obs, "o-", color="tab:blue", ms=10, label="FOM + noise")
        ax.plot(idx, p_rom, "s--", color="tab:red", ms=10, label="GPR")
        ax.set(xlabel=r"pressure point", ylabel=r"$p$")
        ax.set_xticks(idx)
        self._curve_grid(ax)
        self._curve_legend(ax, loc="best")
        self._save_current_figure("pressure_prediction")
        plt.show()
        print(f"Pressure GPR mean relative error against observed pressure: {mean_rel:.6e}")
        print(f"Pressure GPR relative L2 error against observed pressure:    {rel_l2:.6e}")
        return rel_l2

    def plot_prior(self, n: int = 5000) -> None:
        bounds = self._get_parameter_bounds()
        labels = self._get_parameter_labels()
        rates = [self.sigma_noise_prior]
        if self._infer_sigma_bias():
            rates.append(self.sigma_bias_prior)
        n_axes = len(bounds) + len(rates)
        fig, axes = plt.subplots(1, n_axes, figsize=(8 * n_axes, 6), squeeze=False)
        axes = axes.ravel()
        for ax, (lo, hi), label in zip(axes, bounds, labels):
            x = np.linspace(lo, hi, n)
            density = np.full_like(x, 1 / (hi - lo))
            ax.plot(x, density, lw=2)
            ax.set(xlabel=label, ylim=(0, 1.2 * float(np.max(density))))
            ax.grid(True, alpha=0.25)
        axes[0].set_ylabel("Prior density")
        for ax, rate, label in zip(axes[len(bounds):], rates, labels[len(bounds):]):
            x = np.linspace(0, 5 / rate, n)
            ax.plot(x, rate * np.exp(-rate * x), lw=2)
            ax.set(xlabel=label, ylim=(0, 1.1 * rate))
            ax.grid(True, alpha=0.25)
        self._save_current_figure("prior")
        plt.show()

    def plot_trace_all(self, burn_in_ratio: float = 0.6) -> None:
        if getattr(self, "chain", None) is None:
            raise RuntimeError("Call run_mcmc() first.")
        chain = self._to_physical_chain(np.asarray(self.chain, dtype=float))
        burn = int(burn_in_ratio * chain.shape[0])
        labels = self._labels()[:chain.shape[2]]
        fig, axes = plt.subplots(
            chain.shape[2], 1, sharex=True, figsize=(9, 2.4 * chain.shape[2]), squeeze=False
        )
        for i, ax in enumerate(axes.ravel()):
            ax.plot(chain[:, :, i], lw=0.7, alpha=0.55)
            ax.axvline(burn, color="black", ls="--", lw=1)
            ax.set_ylabel(labels[i])
            ax.grid(True, alpha=0.2)
        axes.ravel()[-1].set_xlabel("MCMC step")
        self._save_current_figure("trace_all")
        plt.show()

    def plot_corner(
        self,
        theta_true=None,
        log_scale: bool = False,
        cred_level: float = 0.95,
        label_size: float = 40,
        physical_only: bool = True,
        show_plot: bool = True,
        true_param: bool = True,
    ) -> None:
        """Corner plot with the prior on the diagonals, in the reference style.

        Material parameters carry a uniform prior on ``[lo, hi]``, so the
        diagonals show a flat prior density (red) over their bounds.
        ``log_scale`` is retained for backward compatibility but no longer
        log-transforms any axis (the samples are already in physical units).
        ``physical_only=False`` also appends sigma_noise/sigma_bias.
        ``true_param=True`` (default) draws the ground-truth line (from
        ``theta_true`` or ``self.true_theta``, plus the realized sigma_noise
        when hyperparameters are shown); ``true_param=False`` hides it.
        """
        try:
            import corner
        except ImportError as exc:
            raise ImportError("corner is required for plot_corner().") from exc

        samples_full = self._posterior_samples()
        bounds = self._get_parameter_bounds()
        n_phys = len(bounds)
        ndim = n_phys if physical_only else samples_full.shape[1]
        samples = samples_full[:, :ndim].astype(float).copy()
        labels = list(self._get_parameter_labels()[:ndim])
        if true_param:
            theta_true = theta_true if theta_true is not None else getattr(self, "true_theta", None)
        else:
            theta_true = None

        truths = [None] * ndim
        if theta_true is not None:
            vals = np.atleast_1d(np.asarray(theta_true, dtype=float))
            for i in range(min(n_phys, len(vals))):
                if not np.isfinite(vals[i]):
                    continue
                if not vals[i] > 0:
                    continue
                truths[i] = float(vals[i])
        if true_param and not physical_only:
            names = self._get_parameter_labels(latex=False)[:ndim]

            def _set_truth(pname, value):
                if value is not None and pname in names:
                    truths[names.index(pname)] = float(value)

            curve_sig = getattr(self, "sigma_noise_realized", None)
            _set_truth("sigma_noise", curve_sig)
        fig = plt.figure(figsize=(8 * ndim, 8 * ndim))
        corner.corner(
            samples,
            fig=fig,
            labels=labels,
            label_kwargs=(dict(fontsize=label_size) if label_size else None),
            levels=(0.68, 0.95),
            color="black",
            hist_kwargs=dict(histtype="step", linewidth=2, density=True, color="black"),
            data_kwargs=dict(ms=1.5, alpha=0.2, color="gray"),
        )
        q_lo = 50 * (1 - cred_level)
        q_hi = 50 * (1 + cred_level)
        axes = np.array(fig.axes).reshape((ndim, ndim))
        for i in range(ndim):
            ax = axes[i, i]
            ylim = ax.get_ylim()
            xlo, xhi = ax.get_xlim()
            if truths[i] is not None:
                pad = 0.05 * (xhi - xlo)
                xlo = min(xlo, truths[i] - pad)
                xhi = max(xhi, truths[i] + pad)
                for j in range(i, ndim):
                    axes[j, i].set_xlim(xlo, xhi)
            if i < n_phys:
                lo, hi = bounds[i]
                xs = np.linspace(xlo, xhi, 400)
                dens = np.where((xs >= lo) & (xs <= hi), 1.0 / (hi - lo), 0.0)
                ax.plot(xs, dens, color="red", lw=1.5)
            ci = np.percentile(samples[:, i], [q_lo, q_hi])
            ax.axvline(ci[0], color="black", ls=":", lw=1.5)
            ax.axvline(ci[1], color="black", ls=":", lw=1.5)
            if truths[i] is not None:
                ax.axvline(truths[i], color="blue", ls=":", lw=2, label="ground truth")
                ax.legend(loc="best")
            ax.set_ylim(ylim)

        base = "corner" if log_scale else "corner_physical"
        self._save_current_figure(base if physical_only else base + "_all")
        if show_plot:
            plt.show()
            return
        plt.close(fig)

    def plot_corner_all(self, **kwargs) -> None:
        """Corner plot including the noise/bias hyperparameters."""
        kwargs.setdefault("physical_only", False)
        self.plot_corner(**kwargs)

    @staticmethod
    def _correlated_normal(rng, corr, sigma) -> np.ndarray:
        """Draw a zero-mean sample with covariance ``sigma**2 * corr``."""
        n = corr.shape[0]
        if sigma <= 0:
            return np.zeros(n)
        try:
            L = np.linalg.cholesky(corr + 1e-12 * np.eye(n))
        except np.linalg.LinAlgError:
            w, V = np.linalg.eigh(0.5 * (corr + corr.T))
            L = V @ np.diag(np.sqrt(np.maximum(w, 0)))
        return sigma * (L @ rng.standard_normal(n))

    def _predictive_summary(self, u_avg, x, obs, n_sigma, nsamples_pred, condition_discrepancy):
        """Replicate the observed curve from the posterior draws and summarise it.

        For each posterior draw the ROM prediction is evaluated, a model
        discrepancy is added (from its prior, or conditioned on the residual), and
        measurement noise is added to form a replicated dataset. Returns the band
        (mean, lo, hi) and coverage diagnostics.
        """
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
        l_bias = float(getattr(self, "l_bias", 1) or 1)
        corr = self._bias_correlation_matrix(x, l_bias) if infer_bias else None
        obs_idx = getattr(self, "_obs_indices", None)
        if obs_idx is None:
            raise RuntimeError("_obs_indices is unset. Call load_data() first.")
        obs_idx = np.asarray(obs_idx, dtype=int)

        # Evaluate the surrogate for every draw in one shot, then slice to the
        # observed (thinned) grid; the per-draw loop below only draws noise.
        curves = self._predict_curves(draws[:, :n_material], u_avg_val=u_avg)
        if curves.shape[1] != obs.size and curves.shape[1] >= int(obs_idx[-1]) + 1:
            curves = curves[:, obs_idx]

        latent_mean = np.empty((n_draws, obs.size), dtype=float)
        Y_rep = np.empty((n_draws, obs.size), dtype=float)
        for k in range(n_draws):
            g = curves[k]
            sn, sb = float(sn_draws[k]), float(sb_draws[k])
            if not infer_bias or sb <= 0:
                delta = delta_mean = np.zeros_like(g)
            elif condition_discrepancy:
                A = sb * sb * corr
                Sigma = 0.5 * (A + A.T) + sn * sn * np.eye(len(g))
                Sinv = np.linalg.pinv(Sigma, hermitian=True)
                delta_mean = A @ Sinv @ (obs - g)
                delta = delta_mean + self._correlated_normal(rng, A - A @ Sinv @ A, 1)
            else:
                delta = self._correlated_normal(rng, corr, sb)
                delta_mean = np.zeros_like(g)
            latent_mean[k] = g + delta_mean
            Y_rep[k] = g + delta + rng.standard_normal(len(g)) * sn

        tail = 100 * float(norm.cdf(n_sigma))
        mean_pred = latent_mean.mean(0)
        pred_lo = np.percentile(Y_rep, 100 - tail, axis=0)
        pred_hi = np.percentile(Y_rep, tail, axis=0)
        sigma_total = np.maximum(np.std(Y_rep, axis=0), 1e-8)
        zres = (obs - mean_pred) / sigma_total
        diagnostics = dict(
            coverage=float(np.mean((obs >= pred_lo) & (obs <= pred_hi))),
            rms_z=float(np.sqrt(np.mean(zres ** 2))),
            max_abs_z=float(np.max(np.abs(zres))),
            mean_z=float(np.mean(zres)),
            nominal_coverage=float(2 * norm.cdf(n_sigma) - 1),
        )
        return mean_pred, pred_lo, pred_hi, diagnostics

    def plot_posterior_predictive(
        self,
        n_sigma: float = 1.96,
        nsamples_pred: int = 5000,
        u_avg_val: Optional[float] = None,
        x_filename: str = "curve4_x.txt",
        condition_discrepancy: bool = False,
    ) -> None:
        """Plot the posterior predictive band against the observed curve.

        Mirrors the reference plotting style: a central band (``n_sigma`` deviates,
        1.96 -> 95%), the mean prediction line, the observed points, and printed
        coverage diagnostics stored on ``self.pp_diagnostics``.
        """
        from matplotlib.lines import Line2D

        if getattr(self, "y_obs_matrix", None) is None:
            raise RuntimeError("Call load_data() first.")
        obs = np.atleast_2d(self.y_obs_matrix)[0]
        x = self._observation_x(obs.size, x_filename)
        u_avg = float(u_avg_val if u_avg_val is not None else self.u_avg_obs)

        if getattr(self, "model_family", None) == "tanner":
            import math

            samples = self._posterior_samples()
            n_draws = min(int(nsamples_pred), samples.shape[0])
            rng = np.random.default_rng(0)
            sel = rng.choice(samples.shape[0], size=n_draws, replace=False)
            draws = samples[sel]
            n_material = len(self._get_parameter_bounds())
            i_sn = n_material
            heights = np.asarray(
                [self._tanner_height(float(s[0]), u_avg_val=u_avg) for s in draws],
                dtype=float,
            )
            sn = np.abs(draws[:, i_sn])
            if self._infer_sigma_bias():
                sn = np.sqrt(sn ** 2 + np.abs(draws[:, i_sn + 1]) ** 2)
            y_rep = heights + rng.standard_normal(n_draws) * sn
            tail = 50 * (1 + math.erf(float(n_sigma) / math.sqrt(2)))
            pred_lo = float(np.percentile(y_rep, 100 - tail))
            pred_hi = float(np.percentile(y_rep, tail))
            mean_pred = float(np.mean(heights))
            sigma_total = max(float(np.std(y_rep, ddof=1)), 1e-8)
            obs_height = float(np.max(obs))
            zres = (obs_height - mean_pred) / sigma_total
            diag = dict(
                coverage=float(pred_lo <= obs_height <= pred_hi),
                rms_z=float(abs(zres)),
                max_abs_z=float(abs(zres)),
                mean_z=float(zres),
                nominal_coverage=float(2 * tail / 100 - 1),
            )
            self.pp_diagnostics = diag
            i_max = int(np.argmax(obs))
            fig, ax = self._curve_figure()
            ax.plot(x, obs, color="black", lw=1.8, alpha=0.8, label="_nolegend_")
            ax.scatter([x[i_max]], [obs_height], color="black", s=34, zorder=5, label="observed max")
            ax.axhspan(pred_lo, pred_hi, color="steelblue", alpha=0.25)
            ax.axhline(mean_pred, color="steelblue", lw=1.8)
            self._set_curve_labels(ax)
            self._axis_style(ax, np.append(obs, pred_hi))
            self._curve_grid(ax)
            band = Line2D([0], [0], color="steelblue", lw=6, alpha=0.25, label="Tanner posterior predictive")
            h, lbl = ax.get_legend_handles_labels()
            self._curve_legend(ax, h + [band], lbl + ["Tanner posterior predictive"])
            print(f"\nTanner posterior predictive (U={u_avg:g}, {int(round(100 * diag['nominal_coverage']))}% band)")
            print(
                f"  coverage={diag['coverage']:.1%}  rms_z={diag['rms_z']:.2f}  "
                f"max|z|={diag['max_abs_z']:.2f}  mean_z={diag['mean_z']:.2f}"
            )
            self._save_current_figure("posterior_predictive")
            plt.show()
            return

        mean_pred, pred_lo, pred_hi, diag = self._predictive_summary(
            u_avg, x, obs, n_sigma, nsamples_pred, condition_discrepancy
        )
        self.pp_diagnostics = diag
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
        print(f"\nPosterior predictive (U={u_avg:g}, {int(round(100 * diag['nominal_coverage']))}% band)")
        print(
            f"  coverage={diag['coverage']:.1%}  rms_z={diag['rms_z']:.2f}  "
            f"max|z|={diag['max_abs_z']:.2f}  mean_z={diag['mean_z']:.2f}"
        )
        self._save_current_figure("posterior_predictive")
        plt.show()
        if getattr(self, "use_pressure", False) and getattr(self, "pressure_obs", None) is not None:
            self.plot_pressure_posterior_predictive(n_sigma=n_sigma, nsamples_pred=nsamples_pred)

    def _pressure_predictive_summary(self, n_sigma, nsamples_pred):
        """Replicate the observed pressure vector from the posterior draws.

        The pressure channel enters the likelihood as an i.i.d. Gaussian residual
        with the shared ``sigma_noise`` (no model-discrepancy term, unlike the
        curve). So for each posterior draw the pressure GPR mean is evaluated and
        measurement noise is added to form a replicated pressure dataset. Returns
        the observed vector, the band (mean, lo, hi) and coverage diagnostics.
        """
        from scipy.stats import norm

        p_obs = np.asarray(getattr(self, "pressure_obs", None), dtype=float).ravel()
        if p_obs.size == 0 or not np.all(np.isfinite(p_obs)):
            raise RuntimeError(
                "Observed pressure is missing; call load_data(..., pressure_filename=...)."
            )
        samples = self._posterior_samples()
        n_material = len(self._get_parameter_bounds())
        i_sn = n_material
        rng = np.random.default_rng(0)
        n_draws = min(int(nsamples_pred), samples.shape[0])
        sel = rng.choice(samples.shape[0], size=n_draws, replace=False)
        draws = samples[sel]
        sn_draws = np.abs(draws[:, i_sn])

        # One GPR evaluation for all draws; the loop below only draws noise.
        preds = self._predict_pressure_batch(draws[:, :n_material])
        latent_mean = np.empty((n_draws, p_obs.size), dtype=float)
        Y_rep = np.empty((n_draws, p_obs.size), dtype=float)
        for k in range(n_draws):
            p_pred = preds[k]
            latent_mean[k] = p_pred
            Y_rep[k] = p_pred + rng.standard_normal(p_pred.size) * float(sn_draws[k])

        tail = 100 * float(norm.cdf(n_sigma))
        mean_pred = latent_mean.mean(0)
        pred_lo = np.percentile(Y_rep, 100 - tail, axis=0)
        pred_hi = np.percentile(Y_rep, tail, axis=0)
        sigma_total = np.maximum(np.std(Y_rep, axis=0), 1e-8)
        zres = (p_obs - mean_pred) / sigma_total
        diagnostics = dict(
            coverage=float(np.mean((p_obs >= pred_lo) & (p_obs <= pred_hi))),
            rms_z=float(np.sqrt(np.mean(zres ** 2))),
            max_abs_z=float(np.max(np.abs(zres))),
            mean_z=float(np.mean(zres)),
            nominal_coverage=float(2 * norm.cdf(n_sigma) - 1),
        )
        return p_obs, mean_pred, pred_lo, pred_hi, diagnostics

    def plot_pressure_posterior_predictive(
        self,
        n_sigma: float = 1.96,
        nsamples_pred: int = 5000,
    ) -> None:
        """Plot the posterior predictive band for the pressure channel.

        Mirrors :meth:`plot_posterior_predictive` but for the direct pressure GPR
        forward model: a central band (``n_sigma`` deviates, 1.96 -> 95%), the
        mean prediction line, the observed pressure points, and printed coverage
        diagnostics stored on ``self.pressure_pp_diagnostics``.
        """
        from matplotlib.lines import Line2D

        if not getattr(self, "use_pressure", False):
            raise RuntimeError("plot_pressure_posterior_predictive() requires use_pressure=True.")
        if getattr(self, "pressure_model", None) is None:
            raise RuntimeError("Call build_rom() before plot_pressure_posterior_predictive().")

        p_obs, mean_pred, pred_lo, pred_hi, diag = self._pressure_predictive_summary(
            n_sigma, nsamples_pred
        )
        self.pressure_pp_diagnostics = diag
        idx = np.arange(1, p_obs.size + 1)
        fig, ax = self._curve_figure()
        ax.fill_between(idx, pred_lo, pred_hi, color="steelblue", alpha=0.25)
        ax.plot(idx, mean_pred, color="steelblue", lw=1.5, zorder=4)
        ax.scatter(idx, p_obs, color="black", s=34, zorder=5, alpha=0.8,
                   edgecolors="black", linewidths=0.5, label="_nolegend_")
        ax.set(xlabel=r"pressure point", ylabel=r"$p$")
        ax.set_xticks(idx)
        self._curve_grid(ax)
        band = Line2D([0], [0], color="steelblue", lw=6, alpha=0.25, label="Posterior predictive")
        h, lbl = ax.get_legend_handles_labels()
        self._curve_legend(ax, h + [band], lbl + ["Posterior predictive"], loc="best")
        print(f"\nPressure posterior predictive ({int(round(100 * diag['nominal_coverage']))}% band)")
        print(
            f"  coverage={diag['coverage']:.1%}  rms_z={diag['rms_z']:.2f}  "
            f"max|z|={diag['max_abs_z']:.2f}  mean_z={diag['mean_z']:.2f}"
        )
        self._save_current_figure("pressure_posterior_predictive")
        plt.show()
