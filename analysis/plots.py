import arviz as av
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import yaml
from cosmogrb.universe.survey import Survey
from posterior_predictive_check import PPC
from threeML import TimeSeriesBuilder
from threeML.plugins.OGIPLike import OGIPLike


def _band(E: np.ndarray, Ec: float, piv: float, K: float, alpha: float, beta: float):
    ret = np.zeros(E.shape)
    for i, e in enumerate(E):
        if e < (alpha - beta) * Ec:
            ret[i] = K * (e / piv) ** alpha * np.exp(-e / Ec)
        else:
            ret[i] = (
                K
                * ((alpha - beta) * Ec / piv) ** (alpha - beta)
                * np.exp(beta - alpha)
                * (e / piv) ** beta
            )
    return ret


def plot_band(width: float, siunitx: bool = False):
    K = 10
    alpha = -1
    beta = -2.2
    Ec = 1
    piv = 0.1
    x1 = 1e-2
    x2 = 1e2

    fig, (ax1, ax2) = plt.subplots(
        2, 1, sharex=True, figsize=(width, 3.5), height_ratios=(3, 2)
    )

    E = np.logspace(np.log(x1), np.log(x2), 100)
    ax1.plot(E, K * (E / piv) ** alpha * np.exp(-E / Ec), "b--")
    ax1.plot(E, _band(E, Ec, piv, K, alpha, beta), "b")

    ax1.loglog()
    ax1.set_xlim(x1, x2)
    ax1.set_ylim(1.5e-6, 1e3)
    ax1.legend(title="a)")

    ax2.plot(E, K * 1.602e-6 * E**2 * (E / piv) ** alpha * np.exp(-E / Ec), "b--")
    ax2.plot(E, 1.602e-6 * E**2 * _band(E, Ec, piv, K, alpha, beta), "b")

    ax2.loglog()
    ax2.set_ylim(8e-9, 2e-6)
    ax2.legend(title="b)")

    if siunitx:
        ax1.set_ylabel(
            r"$B$ [$\si{\per\second\per\centi\meter\squared\per\mega\electronvolt}$]"
        )
        ax2.set_xlabel(r"Photon Energy $E_\gamma$ [\si{\mega\electronvolt}]")
        ax2.set_ylabel(r"$E_\gamma^2 B$ [\si{\erg\per\second\per\centi\meter\squared}]")

    fig.align_ylabels([ax1, ax2])
    fig.tight_layout()

    return fig


def plot_light_curve_basics(
    model: str,
    data_folder: str,
    grb_name: str,
    detector: str,
    width: float,
    height: float,
    siunitx: bool = False,
):
    assert "real" in model, "Only works for real data"
    fig, ax = plt.subplots(1, 1, figsize=(width / 2, 0.7 * height))

    with open(data_folder + "self_sample/" + "dload.yml", "r") as f:
        dload = yaml.load(f, Loader=yaml.SafeLoader)

    ts = TimeSeriesBuilder.from_gbm_cspec_or_ctime(
        detector,
        cspec_or_ctime_file=dload[grb_name][detector]["cspec"],
        rsp_file=dload[grb_name][detector]["rsp"],
        verbose=False,
    )
    lc_line = ts.view_lightcurve(-10, 50).gca().lines[0]
    timebins = lc_line.get_xdata()
    cr_lc = lc_line.get_ydata()
    ax.step(timebins, cr_lc)

    ax.set_xlim(-2, 40)
    if siunitx:
        ax.set_xlabel(r"Time $t$ [\si{\second}]")
        ax.set_ylabel(r"Count Rate $n$ [\si{\per\second}]")

    fig.tight_layout()

    return fig


def plot_spectrum_basics(
    model: str,
    data_folder: str,
    grb_name: str,
    detector: str,
    width: float,
    height: float,
    siunitx: bool = False,
):
    assert "real" in model, "Only works for real data"
    fig, ax = plt.subplots(1, 1, figsize=(width / 2, 0.7 * height))

    ogip = OGIPLike(
        detector,
        observation=data_folder + grb_name + "/" + detector + ".pha",
        background=data_folder + grb_name + "/" + detector + "_bak.pha",
        response=data_folder + grb_name + "/" + detector + ".rsp",
        spectrum_number=2,
    )
    sp_lines = ogip.view_count_spectrum().gca().lines
    sp_cr_plot = sp_lines[0]
    energybins = sp_cr_plot.get_xdata()
    cr_sp = sp_cr_plot.get_ydata()
    ax.step(energybins, cr_sp)

    ax.set_xlim(energybins[0], energybins[-1])
    ax.set_xscale("log")
    ax.set_yscale("log")
    if siunitx:
        ax.set_xlabel(r"Energy $E_\gamma$ [\si{\kilo\electronvolt}]")
        ax.set_ylabel(
            r"Diff. Count Rate $n_\mathrm{diff}$ [\si{\per\second\per\kilo\electronvolt}]"
        )

    fig.tight_layout()

    return fig


def plot_weak_light_curve(
    model: str, data_folder: str, survey_name: str, width: float, siunitx: bool = False
):
    assert "simulated" in model, "Only works for simulated data"

    fig, axes = plt.subplots(
        2, 2, sharex=True, sharey="row", figsize=(width, 0.67 * width)
    )
    row = 0
    col = 0

    survey = Survey.from_file(data_folder + survey_name + ".h5")

    grb_name = "SynthGRB_3"

    dets = ["b1", "n0", "n9", "na"]

    for k, v in survey[grb_name].grb.items():
        if k in dets:
            ax = axes[row][col]

            lightcurve = v["lightcurve"]

            lightcurve.display_lightcurve(
                dt=0.5, ax=ax, lw=1, label="Total signal"
            )  # ,color='#25C68C')
            lightcurve.display_source(
                dt=0.5, ax=ax, lw=1, label="Source"
            )  # ,color="#A363DE")
            lightcurve.display_background(
                dt=0.5, ax=ax, lw=1, label="Background"
            )  # , color="#2C342E")

            lc_line = ax.lines[0]
            cr = lc_line.get_ydata()
            timebins = lc_line.get_xdata()
            cr = cr[timebins >= 0]
            timebins = timebins[timebins >= 0]
            cr = cr[timebins <= survey[grb_name].grb.duration]
            timebins = timebins[timebins <= survey[grb_name].grb.duration]
            ax.fill_between(timebins, cr, step="post", alpha=0.3, label="Selection")

            ax.set_xlim(-2, 10)
            if row != 1:
                ax.set_xlabel("")
            elif siunitx:
                ax.set_xlabel(r"Time $t$ [\si{\second}]")
            if col != 0:
                ax.set_ylabel("")
            elif siunitx:
                ax.set_ylabel(r"Count Rate $n$ [\si{\per\second}]")
            ax.legend()

            try:
                number = str(int(k[-1]))
            except ValueError:
                if k[-1] == "a":
                    number = "10"
                elif k[-1] == "b":
                    number = "11"
                else:
                    raise Exception("Unknown detector type")

            if k.startswith("b"):
                title = "BGO " + number
            elif k.startswith("n"):
                title = "NaI " + number
            else:
                raise Exception("Unknown detector type")
            ax.set_title(title)

            if col < 1:
                col += 1
            else:
                row += 1
                col = 0

    return fig


def plot_light_curve(
    model: str,
    data_folder: str,
    width: float,
    siunitx: bool = False,
    **kwargs,
):
    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(width, 0.33 * width))
    grb_name_1, grb_name_2 = kwargs["grb_names"]
    det_1, det_2 = kwargs["det_names"]

    if "simulated" in model:
        survey_name = kwargs["survey_name"]
        survey = Survey.from_file(data_folder + survey_name + ".h5")

        lightcurve = survey[grb_name_1].grb[det_1]["lightcurve"]

        lightcurve.display_lightcurve(dt=0.5, ax=ax1, lw=1, label="Total signal")
        lightcurve.display_source(dt=0.5, ax=ax1, lw=1, label="Source")
        lightcurve.display_background(dt=0.5, ax=ax1, lw=1, label="Background")

        lc_line = ax1.lines[0]
        cr = lc_line.get_ydata()
        timebins = lc_line.get_xdata()
        cr = cr[timebins >= 0]
        timebins = timebins[timebins >= 0]
        cr = cr[timebins <= survey[grb_name_1].grb.duration]
        timebins = timebins[timebins <= survey[grb_name_1].grb.duration]
        ax1.fill_between(timebins, cr, step="post", alpha=0.3, label="Selection")

        ax1.set_xlim(-2, 10)
        ax1.set_ylim(0)
        if siunitx:
            ax1.set_xlabel(r"Time $t$ [\si{\second}]")
            ax1.set_ylabel(r"Count Rate $n$ [\si{\per\second}]")
        ax1.legend()
        title = f"GRB {grb_name_1.replace('SynthGRB_', '')} ("
        if det_1.startswith("b"):
            title += "BGO " + str(int(det_1[-1])) + ")"
        elif det_1.startswith("n"):
            title += "NaI " + str(int(det_1[-1])) + ")"
        else:
            raise Exception("Unknown detector type")
        ax1.set_title(title)

        lightcurve = survey[grb_name_2].grb[det_2]["lightcurve"]

        lightcurve.display_lightcurve(dt=0.5, ax=ax2, lw=1, label="Total signal")
        lightcurve.display_source(dt=0.5, ax=ax2, lw=1, label="Source")
        lightcurve.display_background(dt=0.5, ax=ax2, lw=1, label="Background")

        lc_line = ax2.lines[0]
        cr = lc_line.get_ydata()
        timebins = lc_line.get_xdata()
        cr = cr[timebins >= 0]
        timebins = timebins[timebins >= 0]
        cr = cr[timebins <= survey[grb_name_2].grb.duration]
        timebins = timebins[timebins <= survey[grb_name_2].grb.duration]
        ax2.fill_between(timebins, cr, step="post", alpha=0.3, label="Selection")

        ax2.set_xlim(-2, 10)
        if siunitx:
            ax2.set_xlabel(r"Time $t$ [\si{\second}]")
            ax2.set_ylabel(r"")
        ax2.legend()
        title = f"GRB {grb_name_2.replace('SynthGRB_', '')} ("
        if det_2.startswith("b"):
            title += "BGO " + str(int(det_2[-1])) + ")"
        elif det_2.startswith("n"):
            title += "NaI " + str(int(det_2[-1])) + ")"
        else:
            raise Exception("Unknown detector type")
        ax2.set_title(title)

    else:
        with open(data_folder + "self_sample/" + "dload.yml") as f:
            dload = yaml.load(f, Loader=yaml.SafeLoader)

        with open(data_folder + "grb_names.yml") as f:
            real_names = yaml.load(f, Loader=yaml.SafeLoader)

        ts_1 = TimeSeriesBuilder.from_gbm_cspec_or_ctime(
            det_1,
            cspec_or_ctime_file=dload[grb_name_1][det_1]["cspec"],
            rsp_file=dload[grb_name_1][det_1]["rsp"],
            verbose=False,
        )
        lc_line = ts_1.view_lightcurve(-10, 50).gca().lines[0]
        timebins = lc_line.get_xdata()
        cr = lc_line.get_ydata()
        ax1.step(timebins, cr, label="Total signal")

        source_time = dload[grb_name_1]["source_time"]
        cr = cr[timebins >= source_time[0]]
        timebins = timebins[timebins >= source_time[0]]
        cr = cr[timebins <= source_time[1]]
        timebins = timebins[timebins <= source_time[1]]
        ax1.fill_between(timebins, cr, step="pre", alpha=0.3, label="Selection")

        ax1.set_xlim(-2, 40)
        ax1.set_ylim(0)
        if siunitx:
            ax1.set_xlabel(r"Time $t$ [\si{\second}]")
            ax1.set_ylabel(r"Count Rate $n$ [\si{\per\second}]")
        title = f"{real_names[grb_name_1]} ("
        if det_1.startswith("b"):
            title += "BGO " + str(int(det_1[-1]) + 1) + ")"
        elif det_1.startswith("n"):
            title += "NaI " + str(int(det_1[-1]) + 1) + ")"
        else:
            raise Exception("Unknown detector type")
        ax1.legend()
        ax1.set_title(title)

        ts_2 = TimeSeriesBuilder.from_gbm_cspec_or_ctime(
            det_2,
            cspec_or_ctime_file=dload[grb_name_2][det_2]["cspec"],
            rsp_file=dload[grb_name_2][det_2]["rsp"],
            verbose=False,
        )
        lc_line = ts_2.view_lightcurve(-10, 50).gca().lines[0]
        timebins = lc_line.get_xdata()
        cr = lc_line.get_ydata()
        ax2.step(timebins, cr, label="Total signal")

        source_time = dload[grb_name_2]["source_time"]
        cr = cr[timebins >= source_time[0]]
        timebins = timebins[timebins >= source_time[0]]
        cr = cr[timebins <= source_time[1]]
        timebins = timebins[timebins <= source_time[1]]
        ax2.fill_between(timebins, cr, step="pre", alpha=0.3, label="Selection")

        ax2.set_xlim(-2, 40)
        if siunitx:
            ax2.set_xlabel(r"Time $t$ [\si{\second}]")
            ax2.set_ylabel(r"")
        title = f"{real_names[grb_name_2]} ("
        if det_2.startswith("b"):
            title += "BGO " + str(int(det_2[-1]) + 1) + ")"
        elif det_2.startswith("n"):
            title += "NaI " + str(int(det_2[-1]) + 1) + ")"
        else:
            raise Exception("Unknown detector type")
        ax2.legend()
        ax2.set_title(title)

    fig.align_ylabels([ax1, ax2])
    fig.tight_layout()

    return fig


def plot_ppc(
    model: str,
    data_folder: str,
    ds,
    data: dict,
    res: av.InferenceData,
    width: float,
    siunitx: bool = False,
):
    p = PPC(data, res)
    fig, ax = plt.subplots(1, 2, sharey=True, figsize=(width, 0.33 * width))

    chain = 0
    draws = 500
    detector = 1

    if model.startswith("real"):
        grb_1 = 4
        interval_1 = 0
        grb_2 = 13
        interval_2 = 0

        with open(data_folder + "grb_names.yml") as f:
            real_names = yaml.load(f, Loader=yaml.SafeLoader)

        ax[0].set_title(real_names[list(ds._grbs.keys())[grb_1]])
        ax[1].set_title(real_names[list(ds._grbs.keys())[grb_2]])
    else:
        grb_1 = 4
        interval_1 = 0
        grb_2 = 7
        interval_2 = 0

        ax[0].set_title(f"GRB {grb_1}")
        ax[1].set_title(f"GRB {grb_2}")

    cenergies_1, ppc_counts_1 = p.ppc(
        chain, draws, interval_1, detector, grb=grb_1, interval_in_grb=True
    )
    _, ppc_1s_1, ppc_2s_1 = PPC.ppc_summary(ppc_counts_1)

    obs_mask = data["mask"][
        np.where(data["grb_id"] == grb_1 + 1)[0][interval_1], detector
    ]
    ax[0].stairs(
        data["observed_counts"][
            np.where(data["grb_id"] == grb_1 + 1)[0][interval_1], detector
        ][obs_mask > 0],
        cenergies_1,
        color="#dd8452",
    )
    ax[0].stairs(
        ppc_1s_1[0],
        cenergies_1,
        baseline=ppc_1s_1[1],
        fill=True,
        color="#55a868",
        alpha=0.5,
    )
    ax[0].stairs(
        ppc_2s_1[0],
        cenergies_1,
        baseline=ppc_2s_1[1],
        fill=True,
        color="#4c72b0",
        alpha=0.2,
    )
    if siunitx:
        ax[0].set_xlabel(r"Photon Energy $E_\gamma$ [\si{\kilo\electronvolt}]")
        ax[0].set_ylabel(r"Count Rate $n$ [\si{\per\second}]")
    ax[0].set_xlim(cenergies_1[0], cenergies_1[-1])
    ax[0].loglog()

    cenergies_2, ppc_counts_2 = p.ppc(
        chain, draws, interval_2, detector, grb=grb_2, interval_in_grb=True
    )
    _, ppc_1s_2, ppc_2s_2 = PPC.ppc_summary(ppc_counts_2)

    obs_mask = data["mask"][
        np.where(data["grb_id"] == grb_2 + 1)[0][interval_2], detector
    ]
    ax[1].stairs(
        data["observed_counts"][
            np.where(data["grb_id"] == grb_2 + 1)[0][interval_2], detector
        ][obs_mask > 0],
        cenergies_2,
        color="#dd8452",
    )
    ax[1].stairs(
        ppc_1s_2[0],
        cenergies_2,
        baseline=ppc_1s_2[1],
        fill=True,
        color="#55a868",
        alpha=0.5,
    )
    ax[1].stairs(
        ppc_2s_2[0],
        cenergies_2,
        baseline=ppc_2s_2[1],
        fill=True,
        color="#4c72b0",
        alpha=0.2,
    )
    if siunitx:
        ax[1].set_xlabel(r"Photon Energy $E_\gamma$ [\si{\kilo\electronvolt}]")
    ax[1].set_xlim(cenergies_2[0], cenergies_2[-1])
    ax[1].loglog()

    fig.align_ylabels(ax)
    fig.tight_layout()

    return fig


def plot_corner(
    model: str, res: av.InferenceData, width: float, marginals: bool = False
):
    if marginals:
        fig, axes = plt.subplots(4, 4, figsize=(width, width))
    else:
        fig, axes = plt.subplots(3, 3, figsize=(0.67 * width, 0.67 * width))

    vars = ["alpha", "log_ec", "log_Nrest", "gamma"]
    if "global" in model:
        coords = {f"{i}_dim_0": 0 for i in vars if i != "log_Nrest" and i != "gamma"}
    else:
        coords = {f"{i}_dim_0": 0 for i in vars}
    av.plot_pair(
        res,
        var_names=vars,
        coords=coords,
        ax=axes,
        marginals=marginals,
        kind="kde",
        kde_kwargs={
            "contourf_kwargs": {"cmap": "rocket_r", "colors": None},
            "contour_kwargs": {"cmap": "rocket_r", "colors": None},
        },
    )

    for axs in axes:
        for ax in axs:
            ax.grid(True)

    if marginals:
        axes[-1][0].set_xlabel(r"$\alpha$")
        axes[-1][1].set_xlabel(r"$\log E_C$")
        axes[-1][2].set_xlabel(r"$\log N_\mathrm{rest}$")
        axes[-1][3].set_xlabel(r"$\gamma$")

        axes[0][0].set_ylabel(r"$\alpha$")
        axes[1][0].set_ylabel(r"$\log E_C$")
        axes[2][0].set_ylabel(r"$\log N_\mathrm{rest}$")
        axes[3][0].set_ylabel(r"$\gamma$")
    else:
        axes[-1][0].set_xlabel(r"$\alpha$")
        axes[-1][1].set_xlabel(r"$\log E_C$")
        axes[-1][2].set_xlabel(r"$\log N_\mathrm{rest}$")

        axes[0][0].set_ylabel(r"$\log E_C$")
        axes[1][0].set_ylabel(r"$\log N_\mathrm{rest}$")
        axes[2][0].set_ylabel(r"$\gamma$")

    fig.align_ylabels(np.array(axes).T[0])
    fig.tight_layout()

    return fig


def plot_violin(
    model: str,
    data_folder: str,
    ds,
    res: av.InferenceData,
    width: float,
    height: float,
    ratio: float = 0.5,
    **kwargs,
):
    assert "global" not in model, "Global model is not supported"

    fig_gamma, ax_gamma = plt.subplots(1, 1, figsize=(ratio * width, ratio * height))
    sns.violinplot(data=res.posterior.gamma[0], ax=ax_gamma)

    if model.startswith("real"):
        with open(data_folder + "grb_names.yml") as f:
            real_names = yaml.load(f, Loader=yaml.SafeLoader)
        grb_names = [real_names[i][3:] for i in ds._grbs.keys()]
        ax_gamma.set_xticklabels(grb_names)
        fig_gamma.autofmt_xdate()
    else:
        ax_gamma.plot([-1, len(res.posterior.gamma[0, 0] + 1)], [1.5, 1.5], "--")

        survey_name = kwargs["survey_name"]
        survey = Survey.from_file(data_folder + survey_name + ".h5")
        grb_names = []
        for i in range(len(survey)):
            if f"SynthGRB_{i}" not in list(ds._grbs.keys()):
                continue
            grb_names.append(f"GRB {i}")
        ax_gamma.set_xticklabels(grb_names)
        fig_gamma.autofmt_xdate()

    ax_gamma.set_xlim(-0.5, len(res.posterior.gamma[0, 0]) - 0.5)
    ax_gamma.set_ylabel(r"$\gamma$")

    if "int" in model:
        ax_gamma.set_ylabel("")

    fig_Nrest, ax_Nrest = plt.subplots(
        1, 1, sharex=True, sharey=True, figsize=(ratio * width, ratio * height)
    )
    sns.violinplot(data=res.posterior.log_Nrest[0], ax=ax_Nrest)

    if model.startswith("real"):
        with open(data_folder + "grb_names.yml") as f:
            real_names = yaml.load(f, Loader=yaml.SafeLoader)
        grb_names = [real_names[i][3:] for i in ds._grbs.keys()]
        ax_Nrest.set_xticklabels(grb_names)
        fig_Nrest.autofmt_xdate()
    else:
        ax_Nrest.plot([-1, len(res.posterior.log_Nrest[0, 0] + 1)], [52, 52], "--")

        survey_name = kwargs["survey_name"]
        survey = Survey.from_file(data_folder + survey_name + ".h5")
        grb_names = []
        for i in range(len(survey)):
            if f"SynthGRB_{i}" not in list(ds._grbs.keys()):
                continue
            grb_names.append(f"GRB {i}")
        ax_Nrest.set_xticklabels(grb_names)
        fig_Nrest.autofmt_xdate()

    ax_Nrest.set_xlim(-0.5, len(res.posterior.log_Nrest[0, 0]) - 0.5)
    ax_Nrest.set_ylabel(r"$\log N_\mathrm{rest}$")

    if "int" in model:
        ax_Nrest.set_ylabel("")

    return fig_gamma, fig_Nrest


def plot_trace(
    model: str,
    res,
    width: float,
    divergences: str = "False",
    vars: list = ["alpha", "log_ec", "gamma", "log_Nrest"],
):
    ratio = 0.5

    rows = len(vars)
    fig, axes = plt.subplots(rows, 2, figsize=(ratio * width, 0.8 * width))

    av.plot_trace(res, var_names=vars, divergences=divergences, axes=axes)

    for i in np.array(axes).T[0]:
        i.remove()

    for i, ai in enumerate(np.array(axes).T[1]):
        ai.set_subplotspec(plt.GridSpec(rows, 1)[i : i + 1, 0:1])

    axes[0, 1].set_title(r"$\alpha$")
    axes[1, 1].set_title(r"$\log E_C$")
    axes[2, 1].set_title(r"$\gamma$")
    axes[3, 1].set_title(r"$\log N_\mathrm{rest}$")

    fig.tight_layout()

    return fig


def plot_gc(model: str, data: dict, res: av.InferenceData, width: float, height: float):
    def gc(log_epeak):
        return 52 + 1.5 * (log_epeak - 2)

    ratiow = 0.5
    fig, ax = plt.subplots(1, 1, figsize=(ratiow * width, height / 2))

    cmaps_grbs = [
        "Purples",
        "Blues",
        "Greens",
        "Oranges",
        "Reds",
        "YlOrBr",
        "OrRd",
        "PuRd",
        "RdPu",
        "BuPu",
        "GnBu",
        "PuBu",
        "YlGnBu",
        "PuBuGn",
        "BuGn",
    ]
    cmaps = [cmaps_grbs[data["grb_id"][i] - 1] for i in range(data["N_intervals"])]
    log_epeak = np.array(res.posterior.log_epeak.stack(sample=["chain", "draw"]))
    log_energy_flux = np.array(
        res.posterior.log_energy_flux.stack(sample=["chain", "draw"])
    )

    log_epeak_true = np.log10(1 + np.array(data["z"])) + log_epeak.T
    log_L = np.log10(4 * np.pi * np.square(data["dl"])) + log_energy_flux.T

    # gc_data = pd.DataFrame(
    #     {
    #         "log_epeak": log_epeak_true[0:1000].T.flatten(),
    #         "log_L": log_L[0:1000].T.flatten(),
    #         "grb": np.array(
    #             [data["grb_id"][i // 1000] for i in range(1000 * data["N_intervals"])]
    #         ),
    #     }
    # )

    # sns.kdeplot(
    #     gc_data, x="log_epeak", y="log_L", ax=ax, fill=True, hue="grb", n_levels=5
    # )

    for ep, L, cmap in zip(log_epeak_true.T, log_L.T, cmaps):
        av.plot_kde(
            ep,
            L,
            contour_kwargs={"linewidths": 0, "levels": 5},
            contourf_kwargs={"cmap": cmap, "levels": 5},
        )

    ax.grid(True)

    if "real" in model:
        if "global" in model:
            ax.set_xlim(2.3, 3.4)
        else:
            ax.set_xlim(1.75, 3.8)
        ax.set_ylim(49.6, 54.5)
    else:
        ax.plot([0, 4], [gc(0), gc(4)], "b--")
        ax.set_xlim(1, 2.8)
        ax.set_ylim(50.2, 53.7)

    ax.set_xlabel(r"$\log E_\mathrm{peak}$")
    ax.set_ylabel(r"$\log L$")

    if "int" in model or "global" in model:
        ax.set_ylabel("")
        ax.set_yticklabels("")

    # ax.get_legend().remove()

    return fig


def plot_gc_kde(model: str, res: av.InferenceData, width: float):
    if "global" in model:
        vars = ["gamma", "log_Nrest"][::-1]
    else:
        vars = ["gamma_mu_meta", "log_Nrest_mu_meta"][::-1]

    fig = plt.figure(figsize=(0.5 * width, 0.5 * width))
    gs = mpl.gridspec.GridSpec(
        2, 2, width_ratios=[5, 1], height_ratios=[1, 5], figure=fig
    )
    gs.update(hspace=0, wspace=0)
    ax = np.array([[fig.add_subplot(gs[i, j]) for i in range(2)] for j in range(2)]).T

    av.plot_pair(
        res,
        var_names=vars,
        ax=ax,
        kind="hexbin",
        gridsize=30,
        marginals=True,
        marginal_kwargs=dict(kind="hist", hist_kwargs=dict(bins=40)),
        hexbin_kwargs={"cmap": "rocket_r"},
    )
    # av.plot_kde(log_epeak, log_energy_flux, ax=ax[1,0])
    # av.plot_kde(log_epeak, ax=ax[0,0])
    # av.plot_kde(log_energy_flux, ax=ax[1,1])

    ax[0, 0].set_facecolor((1, 1, 1))
    ax[1, 1].set_facecolor((1, 1, 1))
    ax[0, 1].set_visible(False)

    ax[0, 0].set_xticks([])
    ax[0, 0].set_yticks([])
    ax[1, 1].set_xticks([])
    ax[1, 1].set_yticks([])

    ax[1, 0].grid(True)

    if "global" in model:
        ax[1, 0].set_xlabel(r"$\log N_\mathrm{rest}$")
        ax[1, 0].set_ylabel(r"$\gamma$")
    else:
        ax[1, 0].set_xlabel(r"$\mu_{\log N_\mathrm{rest}}$")
        ax[1, 0].set_ylabel(r"$\mu_\gamma$")

    fig.tight_layout()

    return fig
