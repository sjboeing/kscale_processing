"""
MODULE TO PROCESS K-SCALE DATA FOR COMPARISON TO OBS
NOTE CERES USES LAT/LON INSTEAD OF LATITUDE/LONGITUDE
"""

from pathlib import Path
from glob import glob
import xarray as xr
import xesmf as xe
import numpy as np
from dyamond_setups import DyamondSetup, dyamond_setups, yaml_load

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
        ds = ds.assign_coords(longitude=((ds.longitude + 360) % 360)).sortby(
            "longitude"
        )
    if "latitude" in ds.dims and "longitude" in ds.dims:
        regridder = xe.Regridder(ds, grid_out, "bilinear")
        ds_interp = regridder(ds)
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
    # SORT LONGITUDE
    return ds_interp


def ds_write(ds: xr.Dataset, file_root: str) -> None:
    """SHORTHAND FUNCTION FOR WRITING TO NETCDF"""
    ds.to_netcdf(process_dir + "/" + file_root + ".nc")


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
    ds_write(ds_syn_hourly_mean, "syn_hov_hourlymean_" + setup.start_end_str)
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
    """UNFORTUNATELY LATITUDE AND LONGITUDE NEED TO BE FIXED MANUALLY"""
    longitude_values = ds["longitude"].values
    latitude_values = ds["latitude"].values
    longitude_values_rounded = np.round(longitude_values, decimals=4)
    latitude_values_rounded = np.round(latitude_values, decimals=4)
    new_lat = xr.DataArray(
        latitude_values_rounded,
        dims=("latitude"),
        name="latitude",
        attrs=ds["latitude"].attrs,
    )
    new_lon = xr.DataArray(
        longitude_values_rounded,
        dims=("longitude"),
        name="longitude",
        attrs=ds["longitude"].attrs,
    )
    ds = ds.assign_coords({"longitude": new_lon})  # dict(longitude=new_lon))
    ds = ds.assign_coords({"latitude": new_lat})
    ds = ds.assign_coords(time=ds["time"].dt.floor("h"))
    return ds


def drop_inconsistent_variables(ds: xr.Dataset) -> xr.Dataset:
    """DROP METADATA WHICH IS SLIGHTLY DIFFERENT BETWEEN VARIABLES"""
    ds = ds.drop_vars(
        ["latitude_bnds", "longitude_bnds", "forecast_period"], errors="ignore"
    )
    return ds


def get_processed_data(
    setup: DyamondSetup, dir_list: list
) -> (xr.Dataset, xr.Dataset, xr.Dataset):
    """GET HOVMOELLERS, MEAN BY HOUR AND MEAN"""
    hov_sets = []
    mean_by_hour_sets = []
    mean_sets = []
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
        del ds_mean, ds_month_mean_by_hour, ds_hov, ds_crop, ds_rad
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
    return ds_allhov, ds_hourly_allmean, ds_allmean


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
    ds_allhov, ds_hourly_allmean, ds_allmean = get_processed_data(
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
    # CLEANUP
    del ds_allhov_interp
    del ds_hourly_allmean_interp
    del ds_allmean_interp
    del ds_allhov_hourly_mean


if __name__ == "__main__":
    Path(process_dir).mkdir(parents=True, exist_ok=True)
    for this_setup in dyamond_setups:
        process_ceres(this_setup)
        for this_model in models:
            process_model(this_setup, this_model)
