"""Bayesian inference wrapper that uses ROM as the forward model."""

from pathlib import Path
from typing import Optional, Sequence, Tuple

import numpy as np

from .data_io import DataLoaderMixin
from .likelihood import LikelihoodMixin
from .plotting import PlottingMixin
from .priors import PriorMixin
from .rom import MODEL_PARAMETER_NAMES, ROM, normalize_model_family
from .sampler import Sampler


TANNER_PARAMETER_NAMES = ("N1",)

# ROM input columns of parameters.txt, fixed per model family.
MODEL_PARAM_COLS = {
    "tanner": [0],
    "oldroyd": [0, 1],
    "giesekus": [0, 1, 2],
    "ptt": [0, 1, 2],
}


def normalize_forward_model(model: str) -> str:
    """Return a canonical forward-model name."""
    name = str(model).strip().lower()
    if name == "tanner":
        return name
    return normalize_model_family(name)


def normalize_inference_mode(mode: str) -> str:
    """Return a canonical observation mode for the likelihood."""
    name = str(mode).strip().lower()
    if name not in ("full_curve", "swell_height"):
        raise ValueError("mode must be 'full_curve' or 'swell_height'.")
    return name


class ROMCurve4BayesianInference(
    ROM,
    DataLoaderMixin,
    PriorMixin,
    LikelihoodMixin,
    Sampler,
    PlottingMixin,
):
    """Bayesian inference whose forward model is a trained curve ROM (this object)."""

    def __init__(
        self,
        swell_root: Optional[Path] = None,
        train_data_rels: Optional[Sequence[str]] = None,
        filenames_train: Tuple[str, str] = ("curve4_y.txt", "parameters.txt"),
        filenames_test: Tuple[str, str] = ("curve4_y_test.txt", "parameters_test.txt"),
        rom_method: str = "gpr",
        scaler: Optional[str] = "minmax",
        eps: float = 1e-6,
        lambda_bounds: Optional[Tuple[float, float]] = (1.0, 10.0),
        beta_bounds: Optional[Tuple[float, float]] = (0.1, 0.9),
        alpha_bounds: Optional[Tuple[float, float]] = (0.1, 0.5),
        true_theta: Optional[Sequence[float]] = None,
        sigma_noise_percent: float = 0.0,
        sigma_bias: Optional[str] = None,
        l_bias: float = 1.0,
        thin: int = 1,
        model: str = "auto",
        mode: str = "full_curve",
        n1_bounds: Optional[Tuple[float, float]] = None,
        tanner_ratio_bounds: Optional[Tuple[float, float]] = (0.001, 10.0),
        eta0: Optional[float] = None,
        radius: float = 1.0,
        epsilon_bounds: Optional[Tuple[float, float]] = None,
        eta0_bounds: Optional[Tuple[float, float]] = None,
        use_pressure: bool = False,
        pressure_filename: str = "pressure.txt",
        augment_eta0: bool = False, seed = 42,
    ) -> None:
        model_family = normalize_forward_model(model)
        inference_mode = normalize_inference_mode(mode)
        if model_family == "tanner":
            inference_mode = "swell_height"
        if bool(use_pressure) and model_family == "tanner":
            raise ValueError("use_pressure=True is not supported for model='tanner'.")

        augment_eta0 = model_family == "oldroyd" and (bool(use_pressure) or bool(augment_eta0))
        if sigma_bias not in (None, "infer"):
            raise ValueError("sigma_bias must be None or 'infer'.")
        if true_theta is not None:
            if model_family == "tanner":
                if len(true_theta) != 1:
                    raise ValueError("true_theta must be [N1] for model='tanner'.")
            elif model_family == "auto":
                if len(true_theta) not in (2, 3):
                    raise ValueError("true_theta must be [lambda, beta] or [lambda, beta, alpha].")
            elif augment_eta0:
                if len(true_theta) != 3:
                    raise ValueError(
                        "true_theta must be [lambda, beta, eta_0] for model='oldroyd' "
                        "with use_pressure=True."
                    )
            else:
                names = MODEL_PARAMETER_NAMES[model_family]
                if len(true_theta) != len(names):
                    raise ValueError(f"true_theta must be [{', '.join(names)}] for model='{model_family}'.")

        self.swell_root = (
            Path(swell_root).expanduser().resolve() if swell_root is not None else Path.cwd().resolve()
        )
        self.train_data_rels = list(train_data_rels) if train_data_rels is not None else None
        if model_family == "tanner":
            fixed_n1_bounds = None if n1_bounds is None else tuple(float(v) for v in n1_bounds)
            if tanner_ratio_bounds is None:
                tanner_ratio_bounds = (0.0, 10.0)
            ratio_lo, ratio_hi = (float(tanner_ratio_bounds[0]), float(tanner_ratio_bounds[1]))
            if fixed_n1_bounds is None:
                lo, hi = ratio_lo, ratio_hi
            else:
                lo, hi = fixed_n1_bounds
            if not (hi > lo and lo >= 0.0):
                raise ValueError("n1_bounds must satisfy 0 <= lower < upper.")
            if not (ratio_hi > ratio_lo and ratio_lo >= 0.0):
                raise ValueError("tanner_ratio_bounds must satisfy 0 <= lower < upper.")
            self.data_dir = None
            self.filenames_train = list(filenames_train)
            self.filenames_test = list(filenames_test)
            self.method = "tanner"
            self.scaler = None
            self.eps = float(eps)
            self.model_family = "tanner"
            self.material_parameter_names = list(TANNER_PARAMETER_NAMES)
            self.n_material_params = 1
            self.third_parameter_name = None
            self.param_cols = list(MODEL_PARAM_COLS["tanner"])
            self.lam_idx = -1
            self.beta_idx = -1
            self.alpha_idx = -1
            self.use_uavg = False
            self.n1_bounds = fixed_n1_bounds
            self.tanner_ratio_bounds = (ratio_lo, ratio_hi)
            self.lam_bounds = None
            self.beta_bounds = None
            self.third_parameter_bounds = None
            self.alpha_bounds = None
            self.epsilon_bounds = None
            self.is_trained = True
        else:
            data_rel = self.train_data_rels[0] if self.train_data_rels else "datas/rom_datas"
            data_path = Path(data_rel).expanduser()
            data_path = (data_path if data_path.is_absolute() else self.swell_root / data_path).resolve()
            resolved_filenames_train = filenames_train
            if data_path.is_file() or data_path.name.startswith("curve4_y."):
                resolved_filenames_train = (data_path.name, filenames_train[1])
                data_path = data_path.parent

            # This object IS the ROM forward model; train() populates it in place.
            ROM.__init__(
                self,
                data_dir=data_path,
                filenames_train=resolved_filenames_train,
                filenames_test=filenames_test,
                method=rom_method,
                scaler=scaler,
                eps=eps,
                model_family=model_family,
            )
            # User bounds; those left None are filled from training data during train().
            self.lam_bounds = lambda_bounds
            self.beta_bounds = beta_bounds
            if augment_eta0:
                # Oldroyd-B + pressure: add eta_0 as the third material parameter.
                # Leave its bounds None (unless eta0_bounds is given) so train()
                # fills them from the training-data range.
                self._set_material_parameter_names(("lambda", "beta", "eta_0"))
                self.param_cols = [0, 1, 2]
                self.third_parameter_bounds = (
                    tuple(float(v) for v in eta0_bounds) if eta0_bounds is not None else None
                )
                self.alpha_bounds = self.third_parameter_bounds
                self.epsilon_bounds = None
            else:
                self.third_parameter_bounds = (
                    epsilon_bounds
                    if model_family == "ptt" and epsilon_bounds is not None
                    else alpha_bounds
                )
                self.alpha_bounds = self.third_parameter_bounds
                self.epsilon_bounds = self.third_parameter_bounds if model_family == "ptt" else epsilon_bounds
        self._augment_eta0 = augment_eta0
        self.eta0 = None if eta0 is None else float(eta0)
        self.radius = float(radius)
        if not (self.radius > 0.0):
            raise ValueError("radius must be positive.")
        if self.eta0 is not None and not (self.eta0 > 0.0):
            raise ValueError("eta0 must be positive when provided.")

        self.infer_dir = self.swell_root
        self.true_theta = None if true_theta is None else tuple(float(v) for v in true_theta)
        self.thin = int(thin)
        if self.thin < 1:
            raise ValueError("thin must be a positive integer.")
        self.sigma_noise_percent = float(sigma_noise_percent)
        self.sigma_bias = sigma_bias
        self.l_bias = float(l_bias)
        self.mode = inference_mode

        # Optional pressure channel: a direct GPR (parameters -> pressure vector)
        # trained alongside the curve ROM; its residual joins the likelihood.
        self.use_pressure = bool(use_pressure)
        self.pressure_train_filename = str(pressure_filename)
        self.pressure_model = None
        self.pressure_train = None
        self.pressure_obs = None
        self.pressure_obs_clean = None

        self.y_obs_matrix = None
        self.y_obs_matrix_clean = None
        self.observed_uavgs = None
        self.u_avg_obs = 1.0
        self.obs_x_coords = None
        self.sigma_noise_prior = None
        self.sigma_bias_prior = None
        self.sampler = None
        self.samples = None
        # self.chain_log = None
        # self.log_prob = None
        # self.acceptance_fraction = None
        # self.map_theta = None
        # self.burn_fraction = 0.6
        # self._last_burn_in = None
        # self.nwalkers = None
        self.seed = seed
        
    def build_rom(self) -> None:
        if getattr(self, "model_family", None) == "tanner":
            self.is_trained = True
            print("Tanner analytical forward model selected; no ROM is built.")
            return
        ROM.train(self)
        if getattr(self, "use_pressure", False):
            self._train_pressure_model()

    # ------------------------------------------------------------ pressure GPR
    def _train_pressure_model(self) -> None:
        """Train a direct GPR mapping material parameters to the pressure vector.

        Pressure snapshots are low-dimensional (a handful of points per parameter
        set), so no POD/SVD basis is used here: the GPR regresses the pressure
        vector on the (scaled) material parameters directly. Pressure always uses
        GPR, independent of ``rom_method`` used for the curve ROM.
        """
        path = (Path(self.data_dir) / self.pressure_train_filename).resolve()
        if not path.exists():
            raise FileNotFoundError(f"Missing pressure training file: {path}")
        pressure = np.atleast_2d(np.loadtxt(path, dtype=float))
        params = np.asarray(self.parameters_train, dtype=float)
        if pressure.shape[0] != params.shape[0]:
            raise ValueError(
                f"pressure file has {pressure.shape[0]} rows but parameters_train has "
                f"{params.shape[0]}; they must match."
            )
        x_train = self._transform(params[:, self.param_cols])
        self.pressure_train = pressure
        self.pressure_model = self._train_gpr(x_train, pressure)
        print(
            f"Trained pressure GPR: {pressure.shape[0]} samples x "
            f"{pressure.shape[1]} pressure points (input dim {x_train.shape[1]})."
        )

    def _pressure_input(self, theta: np.ndarray) -> np.ndarray:
        """Scaled GPR input row for material parameters ``theta`` (matches predict())."""
        params = [float(theta[0]), float(theta[1])]
        if self.alpha_idx >= 0:
            params.append(float(theta[2]))
        if getattr(self, "use_uavg", False):
            params.append(float(self.u_avg_obs))
        return self._transform(np.asarray([params], dtype=float))

    def predict_pressure(self, theta: Sequence[float]) -> np.ndarray:
        """Predict the pressure vector for material parameters ``theta`` via GPR."""
        if self.pressure_model is None:
            raise RuntimeError(
                "Pressure model is untrained. Construct with use_pressure=True and call build_rom()."
            )
        x = self._pressure_input(np.asarray(theta, dtype=float))
        pred = np.asarray(self.pressure_model.predict(x), dtype=float).ravel()
        if not np.all(np.isfinite(pred)):
            raise FloatingPointError("Pressure GPR prediction produced NaN/Inf.")
        return pred

    def validate_rom(self) -> None:
        """Report curve ROM test error, and the pressure GPR error when enabled."""
        ROM.validate_rom(self)
        if getattr(self, "use_pressure", False):
            self.validate_pressure()

    def validate_pressure(self) -> dict:
        """Report the pressure GPR relative L2 error over a test set.

        Uses ``pressure_test.txt`` in the ROM ``data_dir`` if it is present and
        matches ``parameters_test``; otherwise falls back to the training pressure
        (in-sample) and says so, so an L2 number is always printed.
        """
        if self.pressure_model is None:
            raise RuntimeError("Call build_rom() with use_pressure=True first.")
        test_path = (Path(self.data_dir) / "pressure_test.txt").resolve()
        if test_path.exists():
            pressure_ref = np.atleast_2d(np.loadtxt(test_path, dtype=float))
            params = np.asarray(self.parameters_test, dtype=float)
            label = "test"
        else:
            pressure_ref = np.atleast_2d(np.asarray(self.pressure_train, dtype=float))
            params = np.asarray(self.parameters_train, dtype=float)
            label = "train (in-sample; add pressure_test.txt for held-out error)"
        if pressure_ref.shape[0] != params.shape[0]:
            raise ValueError("pressure reference rows must match the parameter rows.")

        errors = []
        for p, p_true in zip(params, pressure_ref):
            theta = [float(p[self.lam_idx]), float(p[self.beta_idx])]
            if self.alpha_idx >= 0:
                theta.append(float(p[self.alpha_idx]))
            p_pred = self.predict_pressure(theta)
            errors.append(
                float(np.linalg.norm(p_pred - p_true) / (np.linalg.norm(p_true) + 1e-14))
            )
        mean_err = float(np.mean(errors))
        self._pressure_val_n = len(errors)
        self._pressure_val_err = mean_err
        print(
            f"\nPressure GPR error ({label}): {len(errors)} cases, "
            f"mean rel L2 = {mean_err:.6e}"
        )
        return {"num_cases": len(errors), "relative_l2_error_mean": mean_err}

    def _set_parameter_columns(self) -> None:
        """Fix ROM input columns by model family (see MODEL_PARAM_COLS)."""
        if self.model_family == "auto":
            super()._set_parameter_columns()
            return
        params = np.asarray(self.parameters_train, dtype=float)
        cols = list(MODEL_PARAM_COLS[self.model_family])
        if getattr(self, "_augment_eta0", False):
            cols = [0, 1, 2]  # lambda, beta, eta_0 (eta_0 identified via pressure)
        if params.shape[1] < len(cols):
            suffix = " with use_pressure=True" if getattr(self, "_augment_eta0", False) else ""
            raise ValueError(
                f"model='{self.model_family}'{suffix} needs {len(cols)} parameter "
                f"columns ({', '.join(self.material_parameter_names)}); the training "
                f"file has {params.shape[1]}."
            )
        self.param_cols = list(cols)
        self.alpha_idx = cols[2] if len(cols) == 3 else -1
        self.use_uavg = False
        x = params[:, self.param_cols]
        n_unique = np.unique(x, axis=0).shape[0]
        if n_unique < x.shape[0]:
            raise ValueError(
                f"Duplicate training points for param_cols={self.param_cols} "
                f"({n_unique} unique of {x.shape[0]} rows): the training data varies a "
                f"column outside param_cols. Use training data generated with "
                f"model='{self.model_family}'."
            )

    def _infer_sigma_bias(self) -> bool:
        return self.sigma_bias == "infer"

    def _get_material_parameter_names(self) -> list[str]:
        return list(getattr(self, "material_parameter_names", ["lambda", "beta"]))

    def _get_parameter_bounds(self) -> list[tuple[float, float]]:
        if getattr(self, "model_family", None) == "tanner":
            if self.n1_bounds is not None:
                lo, hi = self.n1_bounds
            else:
                ratio_lo, ratio_hi = self.tanner_ratio_bounds
                tau_w = 4 * float(self.u_avg_obs)
                lo = 2.0 * tau_w * ratio_lo
                hi = 2.0 * tau_w * ratio_hi
            return [(float(lo), float(hi))]
        if not getattr(self, "is_trained", False):
            raise ValueError("ROM bounds are not set. Call train() first.")
        bounds = [self.lam_bounds, self.beta_bounds]
        if self.alpha_idx >= 0:
            bounds.append(self.third_parameter_bounds)
        return [(float(lo), float(hi)) for lo, hi in bounds]

    def _get_ndim(self) -> int:
        return len(self._get_parameter_bounds()) + 1 + int(self._infer_sigma_bias())

    def _get_parameter_labels(self, latex: bool = True) -> list[str]:
        names = self._get_material_parameter_names() + ["sigma_noise"]
        if self._infer_sigma_bias():
            names.append("sigma_bias")
        if not latex:
            return names
        labels = {
            "lambda": r"$\lambda$",
            "beta": r"$\beta$",
            "alpha": r"$\alpha$",
            "epsilon": r"$\epsilon$",
            "eta_0": r"$\eta_0$",
            "N1": r"$N_1$",
            "sigma_noise": r"$\sigma_{\mathrm{noise}}$",
            "sigma_bias": r"$\sigma_{\mathrm{bias}}$",
        }
        return [labels[name] for name in names]

    def _extract_noise_bias(self, theta: np.ndarray) -> tuple[float, float | None]:
        idx = len(self._get_parameter_bounds())
        return float(theta[idx]), float(theta[idx + 1]) if self._infer_sigma_bias() else None

    def _to_physical(self, phi: np.ndarray) -> np.ndarray:
        return np.array(phi, dtype=float, copy=True)

    def _to_physical_chain(self, raw: np.ndarray) -> np.ndarray:
        return self._to_physical(raw)

    # --------------------------------------------------------------- N1 output
    def posterior_material_means(self, use: str = "mean") -> dict[str, float]:
        """Point estimate of the inferred material parameters as a name->value dict.

        Uses the post-burn-in flat chain in ``self.samples`` (physical units).
        ``use='mean'`` returns the posterior mean; ``use='map'`` returns the MAP.
        The hyperparameters (sigma_noise, sigma_bias) are dropped.
        """
        if getattr(self, "samples", None) is None:
            raise RuntimeError("No posterior samples yet; call run_mcmc() first.")
        n_material = len(self._get_parameter_bounds())
        if str(use).lower() == "map":
            if getattr(self, "map_theta", None) is None:
                raise RuntimeError("MAP estimate unavailable; call run_mcmc() first.")
            vec = np.asarray(self.map_theta, dtype=float)[:n_material]
        else:
            vec = np.asarray(self.samples, dtype=float)[:, :n_material].mean(axis=0)
        names = self._get_material_parameter_names()
        return {name: float(v) for name, v in zip(names, vec)}

    def _load_n1_function(self, module_name: str, func_name: str):
        """Import a ``compute_n1_*`` helper from the swell repo root (lazily).

        The ``compute_n1_*_function.py`` modules and ``rheology.py`` live at the
        repository root (the same directory passed as ``swell_root``); make sure
        it is importable, then return the requested function.
        """
        import importlib
        import sys

        roots = [Path(__file__).resolve().parents[1]]  # packages/.. == repo root
        if getattr(self, "swell_root", None) is not None:
            roots.append(Path(self.swell_root))
        for root in roots:
            root_str = str(root)
            if root_str not in sys.path:
                sys.path.insert(0, root_str)
        return getattr(importlib.import_module(module_name), func_name)

    # ------------------------------------------------------------ Tanner model
    def _tanner_swell_ratio(self, n1, u_avg_val, c=None):
        """Tanner's analytical axisymmetric die-swell ratio from ``N1``.

        ``chi = c + (1 + 1/2 S_R^2)^(1/6)`` with the recoverable shear
        ``S_R = N1 / (2 tau_w)`` and the axisymmetric wall shear stress
        ``tau_w = 4 U_avg`` (die radius normalized to 1, unit viscosity),
        matching the ``S_R`` definition in the solver and the N1 prior bounds set
        in :meth:`_get_parameter_bounds`. The additive constant ``c`` (default
        ``0.13``, the Newtonian axisymmetric swell offset) is overridable via the
        ``_tanner_c`` attribute or the argument.
        """
        c = getattr(self, "_tanner_c", 0.13) if c is None else float(c)
        tau_w = 4.0 * float(u_avg_val)
        s_r = float(n1) / (2.0 * tau_w) if tau_w != 0.0 else float("inf")
        return c + (1.0 + 0.5 * s_r**2) ** (1.0 / 6.0)

    def _tanner_height(self, n1, u_avg_val):
        """Predicted swell height for the Tanner forward model.

        The observable is the maximum free-surface height (:meth:`_swell_height`),
        which for a die radius normalized to 1 equals the swell ratio, so this
        returns :meth:`_tanner_swell_ratio` directly.
        """
        return self._tanner_swell_ratio(n1, u_avg_val)

    def compute_N1(
        self,
        rates=None,
        U_avg: Optional[float] = None,
        radius: Optional[float] = None,
        diameter: Optional[float] = None,
        use: str = "mean",
        ptt_kind: str = "exponential",
        tanner_c: Optional[float] = None,
        **rheo_kwargs,
    ) -> dict:
        """First normal-stress difference ``N1 = tau_xx - tau_yy`` at the inferred
        (posterior-mean) material parameters.

        The forward model is chosen from ``self.model_family``:

        * ``oldroyd``  -> :func:`compute_n1_oldroyd_b`  (lambda, beta, eta_0)
        * ``giesekus`` -> :func:`compute_n1_giesekus`   (lambda, beta, alpha, eta_0)
        * ``ptt``      -> :func:`compute_n1_ptt`        (lambda, beta, epsilon, eta_0)
        * ``tanner``   -> Tanner's analytical axisymmetric die-swell relation; here
          ``N1`` is itself the inferred parameter, and the swell ratio is returned
          as ``chi_tanner = tanner_c + (1 + 0.5 (N1 / (2 tau_w))^2)^(1/6)``.

        The shear rate defaults to the observed axisymmetric wall rate
        ``gamma_w = 4 * U_avg / radius`` with ``U_avg = self.u_avg_obs`` and
        ``radius = self.radius``; pass ``rates``, ``U_avg`` or ``radius/diameter``
        to override. ``use='map'`` uses the MAP estimate instead of the mean.
        ``eta_0`` is taken from the inference when it is a fitted parameter
        (Oldroyd-B with ``augment_eta0``/``use_pressure``), otherwise from
        ``self.eta0`` (falling back to ``1.0``). Returns the dict produced by the
        underlying ``compute_n1_*`` helper (keys include ``"N1"`` and ``"rates"``),
        augmented with ``"model"`` and ``"theta_mean"``.
        """
        family = getattr(self, "model_family", None)
        means = self.posterior_material_means(use=use)

        if family == "tanner":
            n1 = float(means["N1"])
            u_avg = float(self.u_avg_obs if U_avg is None else U_avg)
            r = float(self.radius if radius is None else radius)
            tau_w = 4.0 * u_avg  # axisymmetric wall shear stress (eta=1, R=1 convention)
            sr = n1 / (2.0 * tau_w) if tau_w != 0.0 else float("inf")
            chi = self._tanner_swell_ratio(n1, u_avg, c=tanner_c)
            result = {
                "model": "tanner",
                "theta_mean": means,
                "N1": np.atleast_1d(n1),
                "rates": np.atleast_1d(4.0 * u_avg / r),
                "tau_w": tau_w,
                "Sr": sr,
                "chi_tanner": chi,
            }
            self._print_N1(result)
            return result

        # Constitutive models: run the model at the inferred means.
        # Default to the observed axisymmetric wall rate unless overridden.
        if rates is None and U_avg is None and diameter is None:
            U_avg = float(self.u_avg_obs)
            radius = float(self.radius) if radius is None else float(radius)

        lam = float(means["lambda"])
        beta = float(means["beta"])
        eta0 = float(means.get("eta_0", self.eta0 if self.eta0 is not None else 1.0))
        common = dict(
            rates=rates, U_avg=U_avg, radius=radius, diameter=diameter,
            lam=lam, beta=beta, eta0=eta0,
        )

        if family == "oldroyd":
            fn = self._load_n1_function("compute_n1_oldroyd_b_function", "compute_n1_oldroyd_b")
            result = fn(**common, **rheo_kwargs)
        elif family == "giesekus":
            fn = self._load_n1_function("compute_n1_giesekus_function", "compute_n1_giesekus")
            result = fn(alpha=float(means["alpha"]), **common, **rheo_kwargs)
        elif family == "ptt":
            fn = self._load_n1_function("compute_n1_ptt_function", "compute_n1_ptt")
            result = fn(eps=float(means["epsilon"]), ptt=ptt_kind, **common, **rheo_kwargs)
        else:
            raise ValueError(
                f"compute_N1() is not defined for model='{family}'. Use one of "
                "'oldroyd', 'giesekus', 'ptt', or 'tanner'."
            )

        result["model"] = family
        result["theta_mean"] = means
        self._print_N1(result)
        return result

    @staticmethod
    def _print_N1(result: dict) -> None:
        """Print the computed N1 (one line per shear rate)."""
        rates = np.atleast_1d(np.asarray(result["rates"], dtype=float))
        n1 = np.atleast_1d(np.asarray(result["N1"], dtype=float))
        model = result.get("model", "?")
        if n1.size == 1:
            print(f"N1 ({model}) = {n1[0]:.6g}  at shear rate {rates[0]:.6g}")
        else:
            print(f"N1 ({model}):")
            for rate, value in zip(rates, n1):
                print(f"  rate = {rate:.6g}   N1 = {value:.6g}")
