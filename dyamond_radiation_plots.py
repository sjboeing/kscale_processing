"""
ACTUAL PLOTTING SCRIPTS
TODO: ADD SCATTERS, ADD PDFS
MAYBE USE TYPED DICTIONARY TO WRAP DATASETS
"""

from pathlib import Path
from dataclasses import dataclass
from itertools import product
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from dyamond_setups import dyamond_setups, DyamondSetup, yaml_load

config = yaml_load("config.yaml")
process_dir = config["process_dir"]
plots_dir = config["plots_dir"]


@dataclass
class ModelObsComparison:
    """DATACLASS TO HOLD INFORMATION RELATED TO MODEL-OBS INTERCOMPARISON"""

    # pylint: disable=<too-many-instance-attributes>
    obs_variable_name: str
    model_variable_name: str
    label_variable_name: str
    label_units: str
    vmin: float
    vmax: float
    vmin_diff: float
    vmax_diff: float


olr = ModelObsComparison(
    obs_variable_name="toa_lw_all",
    model_variable_name="toa_outgoing_longwave_flux",
    label_variable_name="TOA OLWR",
    label_units="W m^-2",
    vmin=220,
    vmax=340,
    vmin_diff=-80,
    vmax_diff=80,
)


osr = ModelObsComparison(
    obs_variable_name="toa_sw_all",
    model_variable_name="toa_outgoing_shortwave_flux",
    label_variable_name="TOA OSWR",
    label_units="W m^-2",
    vmin=40,
    vmax=260,
    vmin_diff=-80,
    vmax_diff=80,
)

comparisons = [olr, osr]


def ds_open(file_root: str) -> xr.Dataset:
    """SHORTHAND FUNCTION FOR WRITING TO NETCDF"""
    ds = xr.open_dataset(process_dir + "/" + file_root + ".nc")
    return ds


def plot_save(fig: Figure, file_root: str) -> None:
    """SHORTHAND FUNCTION FOR SAVING FIGURES"""
    fig.savefig(plots_dir + "/" + file_root + ".png")


def make_ebaf_plots(
    ds_ebaf_meanmonth: xr.Dataset,
    ds_ebaf_month: xr.Dataset,
    ds_ebaf_gal9: xr.Dataset,
    ds_ebaf_ral32: xr.Dataset,
    setup: DyamondSetup,
    comparison: ModelObsComparison,
) -> None:
    """PLOT ABSOLUTE VALUES FOR COMPARISON TO EBAF"""
    fig, axes = plt.subplots(2, 2, figsize=(16, 8), sharex=True, sharey=True)
    fig.suptitle(comparison.label_variable_name)

    plot_kwargs = {"vmin": comparison.vmin, "vmax": comparison.vmax, "cmap": "viridis"}

    ds_ebaf_month[comparison.obs_variable_name + "_mon"][0].plot.pcolormesh(
        ax=axes[0, 0], **plot_kwargs
    )

    ds_ebaf_meanmonth[comparison.obs_variable_name + "_mon"].plot.pcolormesh(
        ax=axes[0, 1], **plot_kwargs
    )

    ds_ebaf_gal9[comparison.model_variable_name].plot.pcolormesh(
        ax=axes[1, 0], **plot_kwargs
    )

    ds_ebaf_ral32[comparison.model_variable_name].plot.pcolormesh(
        ax=axes[1, 1], **plot_kwargs
    )

    axes[0, 0].set_title("CERES EBAF")
    axes[0, 1].set_title("CERES EBAF climatology")
    axes[1, 0].set_title("GAL9")
    axes[1, 1].set_title("RAL3p2")
    axes[0, 0].set_ylim(-40, 25)

    fig.tight_layout()
    plot_save(fig, comparison.obs_variable_name + "_" + setup.ebaf_str)


def make_ebaf_diff_plots(
    ds_ebaf_meanmonth: xr.Dataset,
    ds_ebaf_month: xr.Dataset,
    ds_ebaf_gal9: xr.Dataset,
    ds_ebaf_ral32: xr.Dataset,
    setup: DyamondSetup,
    comparison: ModelObsComparison,
) -> None:
    """PLOT DIFFERENCES BETWEEN MODEL AND EBAF/EBAF CLIMATOLOGY"""
    gal_minus_ebaf = ds_ebaf_gal9[comparison.model_variable_name].copy(True)
    gal_minus_ebaf[:, :] = (
        ds_ebaf_gal9[comparison.model_variable_name].values
        - ds_ebaf_month[comparison.obs_variable_name + "_mon"][0].values
    )
    gal_minus_ebafmean = ds_ebaf_gal9[comparison.model_variable_name].copy(True)
    gal_minus_ebafmean[:, :] = (
        ds_ebaf_gal9[comparison.model_variable_name].values
        - ds_ebaf_meanmonth[comparison.obs_variable_name + "_mon"].values
    )
    ral_minus_ebaf = ds_ebaf_ral32[comparison.model_variable_name].copy(True)
    ral_minus_ebaf[:, :] = (
        ds_ebaf_ral32[comparison.model_variable_name].values
        - ds_ebaf_month[comparison.obs_variable_name + "_mon"][0].values
    )
    ral_minus_ebafmean = ds_ebaf_ral32[comparison.model_variable_name].copy(True)
    ral_minus_ebafmean[:, :] = (
        ds_ebaf_ral32[comparison.model_variable_name].values
        - ds_ebaf_meanmonth[comparison.obs_variable_name + "_mon"].values
    )

    fig, axes = plt.subplots(2, 2, figsize=(16, 8), sharex=True, sharey=True)
    fig.suptitle(comparison.label_variable_name)

    plot_kwargs = {
        "vmin": comparison.vmin_diff,
        "vmax": comparison.vmax_diff,
        "cmap": "bwr",
    }

    gal_minus_ebaf.plot.pcolormesh(ax=axes[0, 0], **plot_kwargs)

    gal_minus_ebafmean.plot.pcolormesh(ax=axes[0, 1], **plot_kwargs)

    ral_minus_ebaf.plot.pcolormesh(ax=axes[1, 0], **plot_kwargs)

    ral_minus_ebafmean.plot.pcolormesh(ax=axes[1, 1], **plot_kwargs)

    axes[0, 0].set_title("GAL9-EBAF")
    axes[0, 1].set_title("GAL9-EBAF clim")
    axes[1, 0].set_title("RAL3p2-EBAF")
    axes[1, 1].set_title("RAL3p2-EBAF clim")
    axes[0, 0].set_ylim(-40, 25)

    fig.tight_layout()
    plot_save(fig, comparison.obs_variable_name + "_diff_" + setup.ebaf_str)


def make_hov_plots(
    ds_hov_syn: xr.Dataset,
    ds_hov_gal9: xr.Dataset,
    ds_hov_ral32: xr.Dataset,
    setup: DyamondSetup,
    comparison: ModelObsComparison,
) -> None:
    """PLOT HOVMOELLERS"""
    fig, axes = plt.subplots(1, 3, figsize=(12, 8), sharex=True, sharey=True)
    fig.suptitle(comparison.label_variable_name + " Hovmoeller")

    plot_kwargs = {"vmin": comparison.vmin, "vmax": comparison.vmax, "cmap": "viridis"}

    ds_hov_syn[comparison.obs_variable_name + "_1h"].plot.pcolormesh(
        ax=axes[0], **plot_kwargs
    )

    ds_hov_gal9[comparison.model_variable_name].plot.pcolormesh(
        ax=axes[1], **plot_kwargs
    )

    ds_hov_ral32[comparison.model_variable_name].plot.pcolormesh(
        ax=axes[2], **plot_kwargs
    )

    axes[0].set_title("CERES SYN")
    axes[1].set_title("GAL9")
    axes[2].set_title("RAL3p2")

    fig.tight_layout()
    plot_save(fig, comparison.obs_variable_name + "_hov_" + setup.start_end_str)


if __name__ == "__main__":
    Path(plots_dir).mkdir(parents=True, exist_ok=True)
    for this_setup, this_comparison in product(dyamond_setups, comparisons):
        # OPEN RELEVANT FILES
        this_ds_hov_syn = ds_open("syn_hov_" + this_setup.start_end_str)
        this_ds_hov_gal9 = ds_open("GAL9_interp_hov_" + this_setup.start_end_str)
        this_ds_hov_ral32 = ds_open("RAL3p2_interp_hov_" + this_setup.start_end_str)
        this_ds_ebaf_meanmonth = ds_open(
            "ebaf_meanmonth_" + this_setup.ebaf_start_pd.strftime("%m")
        )
        this_ds_ebaf_month = ds_open("ebaf_month_" + this_setup.ebaf_str)
        this_ds_ebaf_gal9 = ds_open("GAL9_interp_month_" + this_setup.ebaf_str)
        this_ds_ebaf_ral32 = ds_open("RAL3p2_interp_month_" + this_setup.ebaf_str)
        # CALL PLOT COMMANDS
        make_ebaf_plots(
            this_ds_ebaf_meanmonth,
            this_ds_ebaf_month,
            this_ds_ebaf_gal9,
            this_ds_ebaf_ral32,
            this_setup,
            this_comparison,
        )
        make_ebaf_diff_plots(
            this_ds_ebaf_meanmonth,
            this_ds_ebaf_month,
            this_ds_ebaf_gal9,
            this_ds_ebaf_ral32,
            this_setup,
            this_comparison,
        )
        make_hov_plots(
            this_ds_hov_syn,
            this_ds_hov_gal9,
            this_ds_hov_ral32,
            this_setup,
            this_comparison,
        )
