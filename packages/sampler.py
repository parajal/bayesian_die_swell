"""Simple MCMC sampling with emcee."""

import numpy as np
import emcee
from sklearn.cluster import KMeans

class Sampler:
    """Run MCMC and compute basic diagnostics."""

    @staticmethod
    def gelman_rubin(chain):
        """chain: (ndim, nwalkers, nsamples)"""
        _, _, n = chain.shape
        W = np.mean(np.var(chain, axis=2, ddof=1), axis=1)
        B = n * np.var(np.mean(chain, axis=2), axis=1, ddof=1)
        var_hat = ((n - 1) / n) * W + B / n
        return np.sqrt(var_hat / W)

    def run_mcmc(self, nwalkers=10, nsteps=5000, burn_fraction=0.3):
        np.random.seed(self.seed)
        rng = np.random.default_rng(self.seed)
        # p0 = self._log_prior(rng, nwalkers)
        # KMeans initialization
        p0 = self.sample_starting_points(nwalkers)
        
        nburn = int(burn_fraction * nsteps)

        sampler = emcee.EnsembleSampler(nwalkers, self._get_ndim(), self.log_posterior)
        sampler.run_mcmc(p0, nsteps, progress=True)

        self.chain = sampler.get_chain()
        self.log_prob = sampler.get_log_prob(discard=nburn, flat=True)

        flat = sampler.get_chain(discard=nburn, flat=True)
        self.samples = self._to_physical(flat)
        self.map_theta = self._to_physical(flat[np.argmax(self.log_prob)])

        self.print_inference_results(nburn)
        return self.samples

    def sample_starting_points(self, nwalkers, pool_factor=20):
        """
        Generate a large pool from the prior and use KMeans
        centers as initial walker locations.
        """

        rng = np.random.default_rng(self.seed)

        pool = self._log_prior(rng, pool_factor * nwalkers)

        # keep only finite posterior points
        mask = np.array(
            [np.isfinite(self.log_posterior(p)) for p in pool]
        )
        pool = pool[mask]

        if len(pool) < nwalkers:
            raise RuntimeError(
                "Not enough valid prior samples to initialize walkers."
            )

        kmeans = KMeans(
            n_clusters=nwalkers,
            n_init=10,
            random_state=self.seed,
        )

        centers = kmeans.fit(pool).cluster_centers_

        return centers

    def print_inference_results(self, nburn):
        """Print posterior mean, std, 95% CI and R-hat."""
        rhat = self.gelman_rubin(self.chain[nburn:].transpose(2, 1, 0))
        flat = self.samples
        names = self._get_parameter_labels(latex=False)

        header = ("Parameter", "Mean", "Std", "Rhat", "95% CI Low", "95% CI High")
        print(f"{header[0]:<15}{header[1]:>12}{header[2]:>12}"
              f"{header[3]:>10}{header[4]:>15}{header[5]:>15}")

        results = {}
        for i, name in enumerate(names):
            x = flat[:, i]
            r = {
                "mean": float(np.mean(x)),
                "std": float(np.std(x, ddof=1)),
                "rhat": float(rhat[i]),
                "ci_low": float(np.percentile(x, 2.5)),
                "ci_high": float(np.percentile(x, 97.5)),
            }
            print(f"{name:<15}{r['mean']:>12.4e}{r['std']:>12.4e}"
                  f"{r['rhat']:>10.3f}{r['ci_low']:>15.4e}{r['ci_high']:>15.4e}")
            results[name] = r

        return results