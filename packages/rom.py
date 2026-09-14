"""Compact non-intrusive curve ROM: train and predict only."""

from pathlib import Path
from typing import Optional, Sequence, Tuple

from scipy.interpolate import RBFInterpolator
import numpy as np

MODEL_PARAMETER_NAMES = {
    "oldroyd": ("lambda", "beta"),
    "giesekus": ("lambda", "beta", "alpha"),
    "ptt": ("lambda", "beta", "epsilon"),
}


def normalize_model_family(model_family: str) -> str:
    """Return a canonical constitutive-model name."""
    name = str(model_family).strip().lower()
    if name not in ("auto", *MODEL_PARAMETER_NAMES):
        allowed = "', '".join(("auto", *MODEL_PARAMETER_NAMES))
        raise ValueError(f"model must be one of '{allowed}'.")
    return name


class ROM:
    """Train an SVD + RBF/GPR reduced model for curve snapshots.

    The Bayesian-inference class inherits this, so the inference object *is* the
    ROM: ``train()`` populates it in place and ``predict()`` is the forward model.
    """

    def __init__(
        self,
        data_dir: "str | Path",
        filenames_train: Sequence[str] = ("curve4_y.txt", "parameters.txt"),
        filenames_test: Sequence[str] = ("curve4_y_test.txt", "parameters_test.txt"),
        method: str = "gpr",
        scaler: Optional[str] = "minmax",
        eps: float = 1e-6,
        kernel: str = "quintic",
        smoothing: float = 0.0,
        epsilon: Optional[float] = None,
        gpr_restarts: int = 20,
        random_state: int = 42,
        model_family: str = "auto",
    ) -> None:
        self.data_dir = Path(data_dir)
        self.filenames_train = list(filenames_train)
        self.filenames_test = list(filenames_test)
        self.method = method.lower()
        self.scaler = None if scaler is None else scaler.lower()
        self.eps = float(eps)
        self.kernel = kernel
        self.smoothing = float(smoothing)
        self.epsilon = None if epsilon is None else float(epsilon)
        self.gpr_restarts = int(gpr_restarts)
        self.random_state = int(random_state)
        self.param_scaler = None
        self.model_family = normalize_model_family(model_family)
        self._auto_material_parameters = self.model_family == "auto"

        initial_model = "giesekus" if self._auto_material_parameters else self.model_family
        self._set_material_parameter_names(MODEL_PARAMETER_NAMES[initial_model])

        self.param_cols = list(range(self.n_material_params))
        self.lam_idx = 0
        self.beta_idx = 1
        self.use_uavg = False
        self.lam_bounds = None
        self.beta_bounds = None
        self.third_parameter_bounds = None
        self.alpha_bounds = None
        self.epsilon_bounds = None
        self.is_trained = False

    def _set_material_parameter_names(self, names: Sequence[str]) -> None:
        names = tuple(str(name) for name in names)
        if len(names) not in (2, 3) or names[:2] != ("lambda", "beta"):
            raise ValueError("material parameters must be [lambda, beta] or [lambda, beta, third].")
        self.material_parameter_names = list(names)
        self.n_material_params = len(names)
        self.third_parameter_name = names[2] if len(names) == 3 else None
        if len(names) == 3:
            self.alpha_idx = 2
            return
        self.alpha_idx = -1

    def train(self) -> None:
        if not all((self.data_dir / f).exists() for f in self.filenames_test):
            self.filenames_test = self.filenames_train.copy()

        self.snapshots_train, self.parameters_train = self._load_pair(self.filenames_train)
        self.snapshots_test, self.parameters_test = self._load_pair(self.filenames_test)

        self._set_parameter_columns()
        self._compute_basis()

        x_train = self._fit_transform(self.parameters_train[:, self.param_cols])
        if self.method == "rbf":
            self.model = self._train_rbf(x_train, self.coeffs)
        elif self.method == "gpr":
            self.model = self._train_gpr(x_train, self.coeffs)
        else:
            raise ValueError("method must be 'rbf' or 'gpr'.")

        self._set_rom_bounds()
        self.is_trained = True

    def predict(
        self,
        lambda_val: float,
        beta_val: float,
        alpha_val: Optional[float] = None,
        u_avg_val: Optional[float] = None,
        epsilon_val: Optional[float] = None,
    ) -> np.ndarray:
        if not self.is_trained:
            raise RuntimeError("Call train() before predict().")

        params = [float(lambda_val), float(beta_val)]
        if self.alpha_idx >= 0:
            if epsilon_val is not None:
                if self.third_parameter_name != "epsilon":
                    raise ValueError("epsilon_val is only valid for model='ptt'.")
                if alpha_val is not None and not np.isclose(float(alpha_val), float(epsilon_val)):
                    raise ValueError("Pass only one value for the third material parameter.")
                third_value = epsilon_val
            else:
                third_value = alpha_val
            if third_value is None:
                arg_name = "epsilon_val" if self.third_parameter_name == "epsilon" else "alpha_val"
                raise ValueError(
                    f"{arg_name} is required because the ROM was trained with {self.third_parameter_name}."
                )
            params.append(float(third_value))

        if self.use_uavg:
            if u_avg_val is None:
                raise ValueError("u_avg_val is required because the ROM was trained with U_avg.")
            params.append(float(u_avg_val))

        x = self._transform(np.asarray([params], dtype=float))
        curve = self._reconstruct(self._predict_coefficients(x)).ravel()
        if not np.all(np.isfinite(curve)):
            raise FloatingPointError("ROM prediction produced NaN/Inf.")
        return curve

    def rom_predict_curve(
        self,
        lambda_val: float,
        beta_val: float,
        u_avg_val: Optional[float] = None,
        alpha_val: Optional[float] = None,
        epsilon_val: Optional[float] = None,
    ) -> np.ndarray:
        """predict() using the loaded observation's U_avg as the default."""
        if u_avg_val is None:
            return self.predict(
                lambda_val,
                beta_val,
                alpha_val=alpha_val,
                u_avg_val=getattr(self, "u_avg_obs", 1),
                epsilon_val=epsilon_val,
            )
        return self.predict(
            lambda_val,
            beta_val,
            alpha_val=alpha_val,
            u_avg_val=float(u_avg_val),
            epsilon_val=epsilon_val,
        )

    def validate_rom(self) -> None:
        if not self.is_trained:
            raise RuntimeError("Call train() before validate_rom().")

        params = np.asarray(self.parameters_test, dtype=float)
        truth = np.asarray(self.snapshots_test, dtype=float)
        u_col = 3 if self.use_uavg and params.shape[1] > 3 else None
        groups = np.unique(np.round(params[:, u_col], 8)) if u_col is not None else [None]

        all_errors = []
        print(f"{'U_avg':>10}  {'n_test':>6}  {'mean rel L2 error':>18}")
        print("-" * 40)
        for u in groups:
            mask = np.isclose(params[:, u_col], u) if u is not None else np.ones(len(params), dtype=bool)
            errors = []
            for p, y_true in zip(params[mask], truth[mask]):
                if self.third_parameter_name == "epsilon" and self.alpha_idx >= 0:
                    third_kwargs = {"epsilon_val": float(p[self.alpha_idx])}
                else:
                    third_kwargs = {
                        "alpha_val": float(p[self.alpha_idx])
                        if self.alpha_idx >= 0
                        else None
                    }
                y_pred = self.predict(
                    float(p[self.lam_idx]),
                    float(p[self.beta_idx]),
                    u_avg_val=(float(u) if u is not None else None),
                    **third_kwargs,
                )
                errors.append(
                    float(np.linalg.norm(y_pred - y_true) / (np.linalg.norm(y_true) + 1e-14))
                )
            all_errors.extend(errors)
            label = f"{u:.4g}" if u is not None else "N/A"
            print(f"{label:>10}  {len(errors):>6d}  {np.mean(errors):>18.4e}")

        self._rom_val_overall_n = len(all_errors)
        self._rom_val_overall_err = float(np.mean(all_errors))
        print("-" * 40)
        print(f"{'Overall':>10}  {self._rom_val_overall_n:>6d}  {self._rom_val_overall_err:>18.4e}")

    def rom_test_error_summary(self) -> dict[str, float]:
        if not hasattr(self, "_rom_val_overall_err"):
            self.validate_rom()
        return {
            "num_test_cases": int(self._rom_val_overall_n),
            "relative_l2_error_mean": float(self._rom_val_overall_err),
        }

    def print_rom_test_error(self) -> dict[str, float]:
        stats = self.rom_test_error_summary()
        print(
            f"\nROM test error: {stats['num_test_cases']} cases, "
            f"mean rel L2 = {stats['relative_l2_error_mean']:.6e}"
        )
        return stats

    def info(self) -> None:
        print(
            f"ROM(model={self.model_family}, method={self.method}, scaler={self.scaler}, "
            f"eps={self.eps:g}, material_parameters={self.material_parameter_names}, "
            f"param_cols={self.param_cols}, trained={self.is_trained})\n"
            f"  data_dir: {self.data_dir}"
        )

    def _load_pair(self, filenames: Sequence[str]) -> Tuple[np.ndarray, np.ndarray]:
        if len(filenames) != 2:
            raise ValueError("filenames must be [curve_file, parameters_file].")
        y = np.atleast_2d(np.loadtxt(self.data_dir / filenames[0], dtype=float))
        p = np.atleast_2d(np.loadtxt(self.data_dir / filenames[1], dtype=float))
        if y.shape[0] != p.shape[0]:
            raise ValueError("Curve and parameter files must have the same number of rows.")
        return y, p

    def _set_parameter_columns(self) -> None:
        params = np.asarray(self.parameters_train, dtype=float)
        if params.shape[1] < 2:
            raise ValueError("Training parameters must contain at least lambda and beta columns.")

        self.param_cols = [0, 1]
        if self._auto_material_parameters:
            if params.shape[1] > 2 and np.ptp(params[:, 2]) > 1e-12:
                self._set_material_parameter_names(MODEL_PARAMETER_NAMES["giesekus"])
                self.param_cols.append(2)
            else:
                self._set_material_parameter_names(MODEL_PARAMETER_NAMES["oldroyd"])
        elif self.n_material_params == 3:
            if params.shape[1] <= 2:
                raise ValueError(
                    f"model='{self.model_family}' requires training parameters "
                    f"[lambda, beta, {self.third_parameter_name}]."
                )
            self.param_cols.append(2)

        if self._auto_material_parameters and params.shape[1] > 3 and np.ptp(params[:, 3]) > 1e-12:
            self.param_cols.append(3)

        self.alpha_idx = 2 if 2 in self.param_cols else -1
        self.use_uavg = 3 in self.param_cols

    def _set_rom_bounds(self) -> None:
        train = np.asarray(self.parameters_train, dtype=float)
        if self.lam_bounds is None:
            self.lam_bounds = (
                float(np.min(train[:, self.lam_idx])),
                float(np.max(train[:, self.lam_idx])),
            )
        if self.beta_bounds is None:
            self.beta_bounds = (
                float(np.min(train[:, self.beta_idx])),
                float(np.max(train[:, self.beta_idx])),
            )
        if self.alpha_idx >= 0 and self.third_parameter_bounds is None and self.alpha_bounds is not None:
            self.third_parameter_bounds = self.alpha_bounds
        if self.alpha_idx >= 0 and self.third_parameter_bounds is None:
            self.third_parameter_bounds = (
                float(np.min(train[:, self.alpha_idx])),
                float(np.max(train[:, self.alpha_idx])),
            )
        if self.alpha_idx >= 0:
            self.alpha_bounds = self.third_parameter_bounds
            if self.third_parameter_name == "epsilon":
                self.epsilon_bounds = self.third_parameter_bounds

    def _compute_basis(self) -> None:
        self.snap_mean = self.snapshots_train.mean(axis=0)
        centered = self.snapshots_train - self.snap_mean
        basis_full, svals, _ = np.linalg.svd(centered.T, full_matrices=False)
        if svals.size == 0:
            raise ValueError("Cannot compute ROM basis from empty snapshots.")
        rel = svals / svals[0] if svals[0] > 0 else svals
        nmodes = max(1, min(int(np.sum(rel >= self.eps)), basis_full.shape[1]))
        self.basis = basis_full[:, :nmodes]
        self.coeffs = centered @ self.basis

    def _fit_transform(self, x: np.ndarray) -> np.ndarray:
        x = np.asarray(x, dtype=float)
        if self.scaler is None:
            self.param_scaler = None
            return x
        if self.scaler == "minmax":
            from sklearn.preprocessing import MinMaxScaler

            self.param_scaler = MinMaxScaler()
        elif self.scaler == "standard":
            from sklearn.preprocessing import StandardScaler

            self.param_scaler = StandardScaler()
        else:
            raise ValueError("scaler must be None, 'minmax', or 'standard'.")
        return self.param_scaler.fit_transform(x)

    def _transform(self, x: np.ndarray) -> np.ndarray:
        x = np.asarray(x, dtype=float)
        if self.param_scaler is None:
            return x
        return self.param_scaler.transform(x)

    def _predict_coefficients(self, x: np.ndarray) -> np.ndarray:
        if self.method == "rbf":
            return np.asarray(self.model(x), dtype=float)
        pred = np.asarray(self.model.predict(x), dtype=float)
        if pred.ndim == 1:
            return pred.reshape(-1, 1)
        return pred

    def _reconstruct(self, coeffs: np.ndarray) -> np.ndarray:
        return coeffs @ self.basis.T + self.snap_mean

    def _train_rbf(self, x: np.ndarray, coeffs: np.ndarray) -> RBFInterpolator:
        kwargs = dict(kernel=self.kernel, smoothing=self.smoothing)
        if self.epsilon is not None:
            kwargs["epsilon"] = self.epsilon
        return RBFInterpolator(x, coeffs, **kwargs)

    def _train_gpr(self, x: np.ndarray, coeffs: np.ndarray):
        from sklearn.gaussian_process import GaussianProcessRegressor
        from sklearn.gaussian_process.kernels import ConstantKernel, RBF

        kernel = ConstantKernel(1.0, (1e-3, 1e3)) * RBF(
            length_scale=np.ones(x.shape[1]), length_scale_bounds=(1e-3, 1e3)
        )
        return GaussianProcessRegressor(
            kernel=kernel,
            alpha=1e-10,
            normalize_y=False,
            n_restarts_optimizer=self.gpr_restarts,
            random_state=self.random_state,
        ).fit(x, coeffs)
