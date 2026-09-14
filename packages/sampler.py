"""MCMC sampling and Gelman-Rubin diagnostics."""

import numpy as np
import emcee
from sklearn.cluster import KMeans


class Sampler:
    """emcee sampling, walker initialisation, and posterior diagnostics."""

    def sample_starting_points(self, nwalkers, seed=42, pool_factor=20):
        """Seed walkers from k-means centres of a prior pool, snapping any
        centre with a non-finite posterior onto the nearest valid draw."""
        rng = np.random.default_rng(seed)
        phi = self._log_prior(rng, pool_factor * nwalkers)
        phi = phi[np.all(np.isfinite(phi), axis=1)]
        if len(phi) < nwalkers:
            raise RuntimeError("Too few finite prior samples to initialize walkers.")

        centers = KMeans(nwalkers, n_init=10, random_state=seed).fit(phi).cluster_centers_
        bad = [k for k, c in enumerate(centers) if not np.isfinite(self.log_posterior(c))]
        if bad:
            good = phi[[np.isfinite(self.log_posterior(p)) for p in phi]]
            if len(good) == 0:
                raise RuntimeError("No prior draw has a finite posterior; check priors and data.")
            for k in bad:
                centers[k] = good[np.argmin(np.linalg.norm(good - centers[k], axis=1))]
        return centers

    def _warmup(self, p0, moves):
        """Run a short chain and restart from its final ensemble."""
        n_warmup = max(1, int(0.1 * self.nsteps))
        sampler = emcee.EnsembleSampler(self.nwalkers, self.ndim, self.log_posterior, moves=moves)
        p1 = np.asarray(sampler.run_mcmc(p0, n_warmup, progress=True).coords, dtype=float)
        if all(np.isfinite(self.log_posterior(row)) for row in p1):
            return p1
        return p0  # emcee never accepts a move to -inf, so p0 is still valid

    def run_mcmc(self, nwalkers=10, nsteps=20000, burn_in_ratio=0.3, seed=42, warmup=True):
        """Run emcee and return posterior samples in physical space."""
        self.nwalkers, self.nsteps = int(nwalkers), int(nsteps)
        self.ndim = self._get_ndim()
        self.burn_fraction = float(burn_in_ratio)
        burn_in = min(int(self.burn_fraction * self.nsteps), self.nsteps - 2)

        moves = [(emcee.moves.StretchMove(a=2), 1.0)]
        p0 = self.sample_starting_points(self.nwalkers, seed)
        if warmup:
            p0 = self._warmup(p0, moves)

        self.sampler = emcee.EnsembleSampler(self.nwalkers, self.ndim, self.log_posterior, moves=moves)
        self.sampler.run_mcmc(p0, self.nsteps, progress=True)

        flat = self.sampler.get_chain(discard=burn_in, flat=True)
        self.chain_log = self.sampler.get_chain()
        self.log_prob = self.sampler.get_log_prob(discard=burn_in, flat=True)
        self.samples = self._to_physical(flat)
        self.acceptance_fraction = self.sampler.acceptance_fraction
        self.map_theta = self._to_physical(flat[np.argmax(self.log_prob)])
        self._last_burn_in = burn_in
        return self.samples

    @staticmethod
    def _gelman_rubin_diag(chain):
        """R-hat for a chain of shape (n_params, n_chains, n_samples)."""
        chain = np.asarray(chain, dtype=float)
        n_params, n_chains, n_samples = chain.shape
        if n_chains < 2 or n_samples < 2:
            return np.full(n_params, np.nan)
        within = np.mean(np.var(chain, axis=2, ddof=1), axis=1)
        between = n_samples * np.var(np.mean(chain, axis=2), axis=1, ddof=1)
        var_hat = (n_samples - 1) / n_samples * within + between / n_samples
        return np.sqrt((var_hat + 1e-12) / (within + 1e-12))

    def compute_diagnostics(self, cred_level=0.95):
        """Per-parameter mean, std, credible interval and R-hat (post burn-in)."""
        if getattr(self, "chain_log", None) is None:
            raise RuntimeError("Call run_mcmc() first.")
        chain = self._to_physical_chain(np.asarray(self.chain_log, dtype=float))
        chain = chain[self._last_burn_in:].transpose(2, 1, 0)  # (ndim, nwalkers, nsteps)
        flat = chain.reshape(chain.shape[0], -1)
        lo, hi = np.percentile(flat, [50 * (1 - cred_level), 50 * (1 + cred_level)], axis=1)
        rhat = self._gelman_rubin_diag(chain)
        names = self._get_parameter_labels(latex=False)
        return {
            names[i]: {
                "mean": float(flat[i].mean()), "std": float(flat[i].std(ddof=1)),
                "ci_low": float(lo[i]), "ci_high": float(hi[i]), "Rhat": float(rhat[i]),
            }
            for i in range(len(names))
        }

    def print_diagnostics(self, cred_level=0.95):
        """Print (and return) the per-parameter posterior summary."""
        diag = self.compute_diagnostics(cred_level)
        pct = int(round(100 * cred_level))
        print(f"\n{'Parameter':<14}{'Mean':>13}{'Std':>13}{'Rhat':>9}"
              f"{f'{pct}% CI low':>15}{f'{pct}% CI high':>15}")
        for name, r in diag.items():
            print(f"{name:<14}{r['mean']:>13.4e}{r['std']:>13.4e}{r['Rhat']:>9.4f}"
                  f"{r['ci_low']:>15.4e}{r['ci_high']:>15.4e}")
        return diag
