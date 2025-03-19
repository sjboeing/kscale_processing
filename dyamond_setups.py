"""This module defines the CTC Dyamond setups and the corresponding observational data."""

from dataclasses import dataclass
import pandas as pd
import yaml


@dataclass
class Basin:
    """Store bain information"""

    name: str
    lat_bounds: tuple[float, float]
    lon_bounds: tuple[float, float]
    lam: str

    @property
    def um_lon_bounds(self) -> tuple[float, float]:
        return tuple((ii + 180.0) % 360.0 - 180.0 for ii in self.lon_bounds)


@dataclass
class DyamondSetup:
    """We store the data in a dataclass object"""

    name: str
    start_date: str
    end_date: str
    ebaf_start_date: str
    ebaf_end_date: str
    syn_file: str
    lat_bounds: tuple[float, float]

    # Define some shorthands, possibly this is overkill
    @property
    def start_pd(self) -> pd.Timestamp:
        return pd.to_datetime(self.start_date)

    @property
    def end_pd(self) -> pd.Timestamp:
        return pd.to_datetime(self.end_date)

    @property
    def ebaf_start_pd(self) -> pd.Timestamp:
        return pd.to_datetime(self.ebaf_start_date)

    @property
    def ebaf_end_pd(self) -> pd.Timestamp:
        return pd.to_datetime(self.ebaf_end_date)

    @property
    def start_str(self) -> str:
        return self.start_pd.strftime("%Y%m%d")

    @property
    def start_end_str(self) -> str:
        return self.start_str + "_" + self.end_pd.strftime("%Y%m%d")

    @property
    def ebaf_str(self) -> str:
        return self.ebaf_start_pd.strftime("%Y%m")


# Ensure model dates does not contain midnight at start and end
# Ensure EBAF dates contain all CERES time steps
dyamond_summer = DyamondSetup(
    name="DYAMOND summer",
    start_date="2016-08-01 00:30",
    end_date="2016-09-10 23:30",
    ebaf_start_date="2016-08-01 00:15",
    ebaf_end_date="2016-08-31 23:45",
    syn_file="CERES_SYN1deg-1H_Terra-Aqua-MODIS_Ed4.1_Subset_20160801-20160930.nc",
    lat_bounds=(5.0, 15.0),
)

dyamond_winter = DyamondSetup(
    name="DYAMOND winter",
    start_date="2020-01-20 00:30",
    end_date="2020-02-28 23:30",
    ebaf_start_date="2020-02-01 00:15",
    ebaf_end_date="2020-02-28 23:45",
    syn_file="CERES_SYN1deg-1H_Terra-Aqua-MODIS_Ed4.1_Subset_20200101-20200229.nc",
    lat_bounds=(-5.0, 5.0),
)

amazon = Basin(
    name="Amazon",
    lat_bounds=(-10.0, 0.0),
    lon_bounds=(295.0, 302.0),
    lam="lam_samerica_km4p4",
)

congo = Basin(
    name="Congo",
    lat_bounds=(-5.0, 5.0),
    lon_bounds=(15.0, 25.0),
    lam="lam_africa_km4p4",
)


def yaml_load(file_name: str) -> dict:
    with open(file_name, "r", encoding="utf-8") as file:
        config = yaml.safe_load(file)
    return config


dyamond_setups = [dyamond_summer, dyamond_winter]
basins = [amazon, congo]
