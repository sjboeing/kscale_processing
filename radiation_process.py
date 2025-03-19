"""
MODULE TO PROCESS K-SCALE DATA FOR COMPARISON TO OBS
NOTE CERES USES LAT/LON INSTEAD OF LATITUDE/LONGITUDE
"""

from pathlib import Path
from glob import glob
from typing import Tuple, Dict, List
import gc
import xarray as xr
import xesmf as xe
import numpy as np
from dyamond_setups import DyamondSetup, Basin, dyamond_setups, yaml_load, basins

config = yaml_load("config.yaml")
obs_dir = config["obs_dir"]
data_dir = config["data_dir"]
ebaf_file = config["ebaf_file"]
process_dir = config["process_dir"]
models = config["models"]
ds_ebaf = xr.open_dataset(obs_dir + ebaf_file)


def interp_to_ceres_grid(ds: xr.Dataset, ds_ceres: xr.Dataset) -> xr.Dataset:
    """DEFINE FUNCTION TO PUT ON CERES GRID"""
    grid_out = xr.Dataset(
        {
            "latitude": (["latitude"], ds_ceres["lat"].values),
            "longitude": (["longitude"], ds_ceres["lon"].values),
        }
    )
    if "longitude" in ds.dims:
        ds = ds.assign_coords(longitude=(ds.longitude + 360) % 360).sortby("longitude")
    if "latitude" in ds.dims and "longitude" in ds.dims:
        regridder = xe.Regridder(ds, grid_out, "bilinear")
        ds_interp = regridder(ds)
        del regridder
    elif "longitude" in ds.dims:
        ceres_lons = xr.DataArray(
            ds_ceres["lon"].values,  # Example longitudes
            dims="longitude",
            coords={"longitude": ds_ceres["lon"].values},
        )
        ds_interp = ds.interp(longitude=ceres_lons, method="linear")
    elif "latitude" in ds.dims:
        ceres_lats = xr.DataArray(
            ds_ceres["lat"].values,  # Example longitudes
            dims="latitude",
            coords={"latitude": ds_ceres["lat"].values},
        )
        ds_interp = ds.interp(latitude=ceres_lats, method="linear")
    else:
        raise KeyError("CERES grid needs latitude and/or longitude")
    # SORT LONGITUDE
    gc.collect()
    return ds_interp


def ds_write(ds: xr.Dataset, file_root: str) -> None:
    """SHORTHAND FUNCTION FOR WRITING TO NETCDF"""
    ds.to_netcdf(
        process_dir + "/" + file_root + ".nc",
        engine="netcdf4",
        encoding={var: {"zlib": True, "complevel": 4} for var in ds.data_vars},
    )
    del ds
    gc.collect()


def process_ceres(setup: DyamondSetup) -> None:
    """GET THE HOVMOELLERS AND MEANS FROM CERES"""
    ds_syn = xr.open_dataset(obs_dir + setup.syn_file)
    # CALCULATE HOVMOELLERS FROM CERES
    ds_syn_subset = ds_syn.sel(
        lat=slice(setup.lat_bounds[0], setup.lat_bounds[1]),
        time=slice(setup.start_pd, setup.end_pd),
    )
    ds_syn_hovmoeller = ds_syn_subset.mean(dim="lat")
    ds_syn_hourly_mean = ds_syn_hovmoeller.groupby("time.hour").mean()
    ds_write(ds_syn_hovmoeller, "syn_hov_" + setup.start_end_str)
    ds_write(ds_syn_hourly_mean, "syn_hov_hourly_mean_" + setup.start_end_str)
    # BASIN DATA
    for basin in basins:
        ds_basin_subset = ds_syn.sel(
            lat=slice(basin.lat_bounds[0], basin.lat_bounds[1]),
            lon=slice(basin.lon_bounds[0], basin.lon_bounds[1]),
            time=slice(setup.start_pd, setup.end_pd),
        )
        ds_basin_series = ds_basin_subset.mean(dim=["lat", "lon"])
        ds_write(ds_basin_series, "syn_" + basin.name + "_" + setup.start_end_str)
    # GET CERES EBAS FOR SAME MONTH AND CLIMATOLOGY FOR THAT MONTHS
    ds_ebaf_months = ds_ebaf.sel(
        time=(ds_ebaf.time.dt.month == setup.ebaf_start_pd.month)
    )
    ds_ebaf_meanmonth = ds_ebaf_months.mean(dim="time")
    ds_ebaf_month = ds_ebaf_months.sel(
        time=(ds_ebaf_months.time.dt.year == setup.ebaf_start_pd.year)
    )
    ds_write(
        ds_ebaf_meanmonth,
        "ebaf_meanmonth_" + setup.ebaf_start_pd.strftime("%m"),
    )
    ds_write(ds_ebaf_month, "ebaf_month_" + setup.ebaf_str)


def fix_um_coordinates(ds: xr.Dataset) -> xr.Dataset:
    """Rounds lat/lon values to avoid floating-point inconsistencies."""
    ds = ds.assign_coords(
        {
            "longitude": np.round(ds["longitude"].values, decimals=4),
            "latitude": np.round(ds["latitude"].values, decimals=4),
            "time": ds["time"].dt.floor("h"),
        }
    )
    return ds


def drop_inconsistent_variables(ds: xr.Dataset) -> xr.Dataset:
    """DROP METADATA WHICH IS SLIGHTLY DIFFERENT BETWEEN VARIABLES"""
    ds = ds.drop_vars(
        ["latitude_bnds", "longitude_bnds", "forecast_period"], errors="ignore"
    )
    return ds


def get_processed_data(
    setup: DyamondSetup, dir_list: List[str]
) -> Tuple[xr.Dataset, xr.Dataset, xr.Dataset, Dict[str, xr.Dataset]]:
    """GET HOVMOELLERS, MEAN BY HOUR AND MEAN"""
    hov_sets: List[xr.Dataset] = []
    mean_by_hour_sets: List[xr.Dataset] = []
    mean_sets: List[xr.Dataset] = []
    basin_sets: Dict[str, List[xr.Dataset]] = {}
    for basin in basins:
        basin_sets[basin.name] = []
    for rad_dir in dir_list:
        # CALCULATE HOVMOELLERS FROM MODEL
        ds_rad = xr.open_mfdataset(rad_dir + "/*.nc")
        ds_rad = fix_um_coordinates(ds_rad)
        ds_rad = drop_inconsistent_variables(ds_rad)
        ds_crop = ds_rad.sel(
            latitude=slice(setup.lat_bounds[0], setup.lat_bounds[1]),
            time=slice(setup.start_pd, setup.end_pd),
        )
        ds_hov = ds_crop.mean(dim="latitude")
        hov_sets.append(ds_hov)
        # FIRST CALCULATE MEAN BY HOUR
        ds_month_sel = ds_rad.sel(time=slice(setup.ebaf_start_pd, setup.ebaf_end_pd))
        ds_month_mean_by_hour = ds_month_sel.groupby("time.hour").mean()
        mean_by_hour_sets.append(ds_month_mean_by_hour)
        # THEN CALCULATE MEAN FROM HOURLY MEANS
        # SO COMPARISON IS FAIR DESPITE MISSING MIDNIGHT AT START/END
        ds_mean = ds_month_mean_by_hour.mean(dim="hour")
        mean_sets.append(ds_mean)
        # BASIN DATA
        for basin in basins:
            ds_basin_subset = ds_rad.sel(
                latitude=slice(basin.lat_bounds[0], basin.lat_bounds[1]),
                longitude=slice(basin.um_lon_bounds[0], basin.um_lon_bounds[1]),
                time=slice(setup.start_pd, setup.end_pd),
            )
            ds_basin_mean = ds_basin_subset.mean(dim=["latitude", "longitude"])
            basin_sets[basin.name].append(ds_basin_mean)
        del (
            ds_mean,
            ds_month_mean_by_hour,
            ds_hov,
            ds_crop,
            ds_rad,
            ds_basin_subset,
            ds_basin_mean,
        )
        gc.collect()
    # WRITE ON NATIVE GRID
    # HOVMOELLERS ON NATIVE GRID
    ds_allhov = xr.merge(hov_sets)
    del hov_sets
    # HOURLY MEANS ON NATIVE GRID
    ds_hourly_allmean = xr.merge(mean_by_hour_sets)
    del mean_by_hour_sets
    # MEANS ON NATIVE GRID
    ds_allmean = xr.merge(mean_sets)
    del mean_sets
    ds_allbasins: Dict[str, xr.Dataset] = {}
    for basin in basins:
        ds_allbasins[basin.name] = xr.merge(basin_sets[basin.name])
    del basin_sets
    gc.collect()
    return ds_allhov, ds_hourly_allmean, ds_allmean, ds_allbasins


def process_model(setup: DyamondSetup, model: str) -> None:
    """PROCESS THE DYAMOND DATA"""
    base_path = (
        data_dir
        + "/outdir_"
        + setup.start_str
        + "T0000Z/DMn1280GAL9/channel_n2560_"
        + model
    )
    olwr_files_path = base_path + "/single_olwr"
    qtot_files_path = base_path + "/single_qtot"
    extra_rad_path = (
        data_dir
        + "/extra_rad_diagnostics_n1280GAL9dm/outdir_"
        + setup.start_str
        + "T0000Z/CTC_N2560_"
        + model
        + "/"
    )
    # THE EXTRA VARIABLES ARE IN DIFFERENT DIRECTORIES
    extra_rad_dirs = glob(extra_rad_path + "/*")
    ds_allhov, ds_hourly_allmean, ds_allmean, ds_allbasins = get_processed_data(
        setup, [olwr_files_path, qtot_files_path] + extra_rad_dirs
    )
    ds_write(ds_allhov, model + "_hov_" + setup.start_end_str)
    ds_write(ds_hourly_allmean, model + "_hourly_month_" + setup.ebaf_str)
    ds_write(ds_allmean, model + "_month_" + setup.ebaf_str)
    # INTERPOLATE THESE TO CERES GRID
    # HOVMOELLERS ON CERES GRID
    ds_allhov_interp = interp_to_ceres_grid(ds_allhov, ds_ebaf)
    ds_write(ds_allhov_interp, model + "_interp_hov_" + setup.start_end_str)
    del ds_allhov
    # HOURLY MEANS ON CERES GRID
    ds_hourly_allmean_interp = interp_to_ceres_grid(ds_hourly_allmean, ds_ebaf)
    ds_write(
        ds_hourly_allmean_interp,
        model + "_interp_hourly_month_" + setup.ebaf_str,
    )
    del ds_hourly_allmean
    # MEANS ON CERES GRID
    ds_allmean_interp = interp_to_ceres_grid(ds_allmean, ds_ebaf)
    ds_write(ds_allmean_interp, model + "_interp_month_" + setup.ebaf_str)
    del ds_allmean
    # EXTRA: MEAN BY HOUR HOVMOELLERS ON CERES GRID
    ds_allhov_hourly_mean = ds_allhov_interp.groupby("time.hour").mean()
    ds_write(
        ds_allhov_hourly_mean,
        model + "_interp_hov_hourly_mean_" + setup.start_end_str,
    )
    for basin in basins:
        ds_write(
            ds_allbasins[basin.name],
            model + "_" + basin.name + "_" + setup.start_end_str,
        )
    # CLEANUP
    del ds_allhov_interp
    del ds_hourly_allmean_interp
    del ds_allmean_interp
    del ds_allhov_hourly_mean
    del ds_allbasins
    gc.collect()


def process_lam(setup: DyamondSetup, basin: Basin) -> None:
    """PROCESS THE LAM DATA"""
    base_path = (
        data_dir
        + "/outdir_"
        + setup.start_str
        + "T0000Z/DMn1280GAL9/"
        + basin.lam
        + "_RAL3p2"
    )
    olwr_files_path = base_path + "/single_olwr"
    qtot_files_path = base_path + "/single_qtot"
    # THE EXTRA VARIABLES ARE IN DIFFERENT DIRECTORIES
    this_basin_sets: List[xr.Dataset] = []
    for rad_dir in [olwr_files_path, qtot_files_path]:
        ds_rad = xr.open_mfdataset(rad_dir + "/*.nc")
        ds_rad = fix_um_coordinates(ds_rad)
        ds_rad = drop_inconsistent_variables(ds_rad)
        ds_basin_subset = ds_rad.sel(
            latitude=slice(basin.lat_bounds[0], basin.lat_bounds[1]),
            longitude=slice(basin.um_lon_bounds[0], basin.um_lon_bounds[1]),
            time=slice(setup.start_pd, setup.end_pd),
        )
        ds_basin_mean = ds_basin_subset.mean(dim=["latitude", "longitude"])
        this_basin_sets.append(ds_basin_mean)
        del ds_basin_mean, ds_basin_subset, ds_rad
    ds_basin = xr.merge(this_basin_sets)
    ds_write(
        ds_basin,
        "LAM_" + basin.name + "_" + setup.start_end_str,
    )
    # CLEANUP
    del ds_basin
    gc.collect()


if __name__ == "__main__":
    Path(process_dir).mkdir(parents=True, exist_ok=True)
    for this_setup in dyamond_setups:
        process_ceres(this_setup)
        for this_basin in basins:
            process_lam(this_setup, this_basin)
        for this_model in models:
            process_model(this_setup, this_model)
