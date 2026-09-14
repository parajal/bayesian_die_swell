from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from rheology import RheoData, VEModel


def _finish(ax, *, title, xlabel, ylabel, xscale="linear", yscale="linear", legend=None):
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    ax.tick_params(labelsize=12)
    if legend is not None:
        ax.legend(legend)


def rheoplot(
    timetype: str,
    rheodata: RheoData,
    vemodel: VEModel,
    flowtype: int,
    plottype: str,
):
    """Python equivalent of rheoplot.m.

    Returns a list of Matplotlib figures so callers can save or show them.
    """
    if rheodata.stress is None:
        raise ValueError("rheodata.stress must be set")

    sxx = rheodata.stress[0, :]
    sxy = rheodata.stress[1, :]
    syy = rheodata.stress[3, :]
    szz = rheodata.stress[5, :]
    figs = []

    if flowtype == 1:
        if timetype == "startup":
            if rheodata.time is None or rheodata.rate_for_startup is None:
                raise ValueError("startup plots need rheodata.time and rate_for_startup")
            time = rheodata.time
            rate = rheodata.rate_for_startup

            if plottype == "visc":
                fig, ax = plt.subplots()
                ax.plot(time, sxy / rate, linewidth=2)
                _finish(
                    ax,
                    title=r"Transient shear viscosity $\eta(t)$",
                    xlabel=r"$t$",
                    ylabel=r"$\eta$",
                )
                figs.append(fig)

                fig, ax = plt.subplots()
                ax.plot(time, (sxx - syy) / rate**2, linewidth=2)
                _finish(
                    ax,
                    title=r"Transient first normal stress coefficient $\Psi_1(t)$",
                    xlabel=r"$t$",
                    ylabel=r"$\Psi_1$",
                )
                figs.append(fig)

                fig, ax = plt.subplots()
                ax.plot(time, (sxx - syy) / sxy, linewidth=2)
                _finish(
                    ax,
                    title=r"Transient stress ratio $S(t)=N_1(t)/\tau_{xy}(t)$",
                    xlabel=r"$t$",
                    ylabel=r"$S$",
                )
                figs.append(fig)

                fig, ax = plt.subplots()
                ax.plot(time, (syy - szz) / rate**2, linewidth=2)
                _finish(
                    ax,
                    title=r"Transient second normal stress coefficient $\Psi_2(t)$",
                    xlabel=r"$t$",
                    ylabel=r"$\Psi_2$",
                )
                figs.append(fig)

            elif plottype == "stress":
                fig, ax = plt.subplots()
                ax.plot(time, sxy, linewidth=2)
                ax.plot(time, sxx - syy, linewidth=2)
                ax.plot(time, syy - szz, linewidth=2)
                _finish(
                    ax,
                    title="Transient stress",
                    xlabel=r"$t$",
                    ylabel=r"$\tau$",
                    legend=[
                        r"$\tau_{xy}$",
                        r"$N_1 = \tau_{xx}-\tau_{yy}$",
                        r"$N_2 = \tau_{yy}-\tau_{zz}$",
                    ],
                )
                figs.append(fig)

        elif timetype == "steady":
            if rheodata.rates is None:
                raise ValueError("steady plots need rheodata.rates")
            rates = rheodata.rates

            if plottype == "visc":
                fig, ax = plt.subplots()
                ax.plot(rates, sxy / rates, linewidth=2)
                _finish(
                    ax,
                    title=r"Steady shear viscosity $\eta(\dot{\gamma})$",
                    xlabel=r"$\dot{\gamma}$",
                    ylabel=r"$\eta$",
                    xscale="log",
                    yscale="log",
                )
                figs.append(fig)

                fig, ax = plt.subplots()
                ax.plot(rates, (sxx - syy) / rates**2, linewidth=2)
                _finish(
                    ax,
                    title=r"Steady first normal stress coefficient $\Psi_1(\dot{\gamma})$",
                    xlabel=r"$\dot{\gamma}$",
                    ylabel=r"$\Psi_1$",
                    xscale="log",
                    yscale="log",
                )
                figs.append(fig)

                fig, ax = plt.subplots()
                ax.plot(rates, (sxx - syy) / sxy, linewidth=2)
                _finish(
                    ax,
                    title=r"Steady stress ratio $S(\dot{\gamma})=N_1(\dot{\gamma})/\tau_{xy}(\dot{\gamma})$",
                    xlabel=r"$\dot{\gamma}$",
                    ylabel=r"$S$",
                    xscale="log",
                    yscale="log",
                )
                figs.append(fig)

                fig, ax = plt.subplots()
                ax.plot(rates, (syy - szz) / rates**2, linewidth=2)
                _finish(
                    ax,
                    title=r"Steady second normal stress coefficient $\Psi_2(\dot{\gamma})$",
                    xlabel=r"$\dot{\gamma}$",
                    ylabel=r"$\Psi_2$",
                    xscale="log",
                    yscale="log",
                )
                figs.append(fig)

            elif plottype == "stress":
                n2 = syy - szz
                fig, ax = plt.subplots()
                ax.plot(rates, sxy, linewidth=2)
                ax.plot(rates, sxx - syy, linewidth=2)
                ax.plot(rates, np.where(n2 < 0, -n2, n2), linewidth=2)
                n2_label = (
                    r"$-N_2 = \tau_{zz}-\tau_{yy}$"
                    if np.any(n2 < 0)
                    else r"$N_2 = \tau_{yy}-\tau_{zz}$"
                )
                _finish(
                    ax,
                    title="Steady-state stress",
                    xlabel=r"$\dot{\gamma}$",
                    ylabel=r"$\tau$",
                    xscale="log",
                    yscale="log",
                    legend=[r"$\tau_{xy}$", r"$N_1 = \tau_{xx}-\tau_{yy}$", n2_label],
                )
                figs.append(fig)

        elif timetype == "startup_stress":
            if rheodata.time is None:
                raise ValueError("startup_stress plots need rheodata.time")

            if plottype == "strain":
                fig, ax = plt.subplots()
                ax.plot(rheodata.time, rheodata.strain, linewidth=2)
                _finish(
                    ax,
                    title="Transient shear strain for imposed shear stress",
                    xlabel=r"$t$",
                    ylabel=r"$\gamma$",
                )
                figs.append(fig)
            elif plottype == "rate":
                fig, ax = plt.subplots()
                ax.plot(rheodata.time, rheodata.rates, linewidth=2)
                _finish(
                    ax,
                    title="Transient strain-rate for imposed shear stress",
                    xlabel=r"$t$",
                    ylabel=r"$\dot{\gamma}$",
                )
                figs.append(fig)

    elif flowtype in (2, 3):
        if timetype == "startup":
            if rheodata.time is None or rheodata.rate_for_startup is None:
                raise ValueError("startup plots need rheodata.time and rate_for_startup")
            time = rheodata.time
            rate = rheodata.rate_for_startup

            if plottype == "visc":
                fig, ax = plt.subplots()
                ax.plot(time, (sxx - syy) / rate, linewidth=2)
                _finish(
                    ax,
                    title=r"Transient extensional viscosity $\eta_E(t)$",
                    xlabel=r"$t$",
                    ylabel=r"$\eta_E$",
                )
                figs.append(fig)
            elif plottype == "stress":
                fig, ax = plt.subplots()
                ax.plot(time, sxx - syy, linewidth=2)
                _finish(
                    ax,
                    title=r"Transient normal stress $N_1(t)=\tau_{xx}-\tau_{yy}$",
                    xlabel=r"$t$",
                    ylabel=r"$N_1$",
                )
                figs.append(fig)

        elif timetype == "steady":
            if rheodata.rates is None:
                raise ValueError("steady plots need rheodata.rates")
            rates = rheodata.rates

            if plottype == "visc":
                fig, ax = plt.subplots()
                ax.plot(rates, (sxx - syy) / rates, linewidth=2)
                _finish(
                    ax,
                    title=r"Steady extensional viscosity $\eta_E(\dot{\epsilon})$",
                    xlabel=r"$\dot{\epsilon}$",
                    ylabel=r"$\eta_E$",
                    xscale="log",
                    yscale="log",
                )
                figs.append(fig)
            elif plottype == "stress":
                fig, ax = plt.subplots()
                ax.plot(rates, (sxx - syy) / rates, linewidth=2)
                _finish(
                    ax,
                    title=r"Steady normal stress $N_1(\dot{\epsilon})=\tau_{xx}-\tau_{yy}$",
                    xlabel=r"$\dot{\epsilon}$",
                    ylabel=r"$N_1$",
                    xscale="log",
                    yscale="log",
                )
                figs.append(fig)

        elif timetype == "startup_stress":
            if rheodata.time is None:
                raise ValueError("startup_stress plots need rheodata.time")

            if plottype == "strain":
                fig, ax = plt.subplots()
                ax.plot(rheodata.time, rheodata.strain, linewidth=2)
                _finish(
                    ax,
                    title="Transient elongational strain for imposed normal stress",
                    xlabel=r"$t$",
                    ylabel=r"$\epsilon$",
                )
                figs.append(fig)
            elif plottype == "rate":
                fig, ax = plt.subplots()
                ax.plot(rheodata.time, rheodata.rates, linewidth=2)
                _finish(
                    ax,
                    title="Transient elongational rate for imposed normal stress",
                    xlabel=r"$t$",
                    ylabel=r"$\dot{\epsilon}$",
                )
                figs.append(fig)

    return figs
