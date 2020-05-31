# -*- coding: utf-8 -*-
"""
Gas and Dust Python Calculator.

Unit tests for gdpyc
"""
import numpy as np
import pytest
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.units.quantity import Quantity
from astropy.utils.exceptions import AstropyWarning


from ..core import DustMap, GasMap


def test_nh_badmap():
    coords = SkyCoord(ra=136.0, dec=44.0, unit="deg")
    with pytest.raises(ValueError):
        GasMap.nh(coords, nhmap="BADMAP")


def test_nh_DL():
    coords = SkyCoord(ra=136.0, dec=44.0, unit="deg")
    nh = GasMap.nh(coords, nhmap="DL")

    assert isinstance(nh, Quantity)
    assert nh.value > 1e19
    assert nh.value < 3e22


def test_nh_LAB():
    coords = SkyCoord(ra=136.0, dec=44.0, unit="deg")
    nh = GasMap.nh(coords, nhmap="LAB")

    assert isinstance(nh, Quantity)
    assert nh.value > 1e19
    assert nh.value < 3e22


def test_nh_HI4PI():
    coords = SkyCoord(ra=136.0, dec=44.0, unit="deg")
    nh = GasMap.nh(coords, nhmap="HI4PI")

    assert isinstance(nh, Quantity)
    assert nh.value > 1e19
    assert nh.value < 3e22


def test_nh_DL_hires():
    coords = SkyCoord(ra=136.0, dec=44.0, unit="deg")
    nh = GasMap.nh(coords, nhmap="DL", hires=True)

    assert isinstance(nh, Quantity)
    assert nh.value > 1e19
    assert nh.value < 3e22


def test_nh_LAB_hires():
    coords = SkyCoord(ra=136.0, dec=44.0, unit="deg")
    nh = GasMap.nh(coords, nhmap="LAB", hires=True)

    assert isinstance(nh, Quantity)
    assert nh.value > 1e19
    assert nh.value < 3e22


def test_nh_HI4PI_hires():
    coords = SkyCoord(ra=136.0, dec=44.0, unit="deg")
    nh = GasMap.nh(coords, nhmap="HI4PI", hires=True)

    assert isinstance(nh, Quantity)
    assert nh.value > 1e19
    assert nh.value < 3e22


def test_nh_array():
    coords = SkyCoord(
        ra=np.linspace(0, 360, num=20), dec=np.linspace(-90, 90, num=20), unit="deg"
    )
    nh = GasMap.nh(coords)

    assert len(nh) == len(coords)


def test_nhtotal():
    coords = SkyCoord(ra=136.0, dec=44.0, unit="deg")
    nh = GasMap.nhtotal(coords)

    assert isinstance(nh, Quantity)
    assert nh.value > 1e19


def test_nhtotal_hires():
    coords = SkyCoord(ra=136.0, dec=44.0, unit="deg")
    nh = GasMap.nhtotal(coords, hires=True)

    assert isinstance(nh, Quantity)
    assert nh.value > 1e19


def test_nhtotal_array():
    coords = SkyCoord(
        ra=np.linspace(0, 360, num=20), dec=np.linspace(-90, 90, num=20), unit="deg"
    )
    nh = GasMap.nhtotal(coords)

    assert len(nh) == len(coords)


def test_nhf_DL():
    coords = SkyCoord(ra=136.0, dec=44.0, unit="deg")
    nh = GasMap.nhf(coords, nhmap="DL")

    assert isinstance(nh, Quantity)
    assert nh.value > 1e19
    assert nh.value < 3e22


def test_nhf_LAB():
    coords = SkyCoord(ra=136.0, dec=44.0, unit="deg")
    nh = GasMap.nhf(coords, nhmap="LAB")

    assert isinstance(nh, Quantity)
    assert nh.value > 1e19
    assert nh.value < 3e22


def test_nhf_array():
    coords = SkyCoord(
        ra=np.linspace(0, 360, num=20), dec=np.linspace(-89.5, 89.5, num=20), unit="deg"
    )
    nh = GasMap.nhf(coords)

    assert len(nh) == len(coords)


def test_ebv_badmap():
    coords = SkyCoord(ra=136.0, dec=44.0, unit="deg")
    with pytest.raises(ValueError):
        DustMap.ebv(coords, dustmap="BADMAP")


def test_ebv_SFD():
    coords = SkyCoord(ra=136.0, dec=44.0, unit="deg")
    ebv = DustMap.ebv(coords, dustmap="SFD")

    assert ebv > 0
    assert ebv < 60


def test_ebv_Planck13():
    coords = SkyCoord(ra=136.0, dec=44.0, unit="deg")
    ebv = DustMap.ebv(coords, dustmap="Planck13")

    assert ebv > 0
    assert ebv < 125


def test_ebv_array():
    coords = SkyCoord(
        ra=np.linspace(0, 360, num=20), dec=np.linspace(-90, 90, num=20), unit="deg"
    )
    ebv = DustMap.ebv(coords)

    assert len(ebv) == len(coords)


def test_extinction_SFD():
    coords = SkyCoord(ra=136.0, dec=44.0, unit="deg")
    ext = DustMap.extinction(coords, dustmap="SFD", filters="SDSS_r")

    assert isinstance(ext, Table)
    assert ext["SDSS_r"][0] > 0


def test_extinction_Planck13():
    coords = SkyCoord(ra=136.0, dec=44.0, unit="deg")
    with pytest.warns(AstropyWarning):
        ext = DustMap.extinction(coords, dustmap="Planck13", filters="SDSS_r")

    assert isinstance(ext, Table)
    assert ext["SDSS_r"][0] > 0


def test_extinction_array():
    coords = SkyCoord(
        ra=np.linspace(0, 360, num=20), dec=np.linspace(-90, 90, num=20), unit="deg"
    )
    ext = DustMap.extinction(coords, filters="SDSS_r")

    assert len(ext) == len(coords)


def test_extinction_array_nfilters():
    filters = ["SDSS_u", "SDSS_g", "SDSS_r"]
    coords = SkyCoord(
        ra=np.linspace(0, 360, num=20), dec=np.linspace(-90, 90, num=20), unit="deg"
    )
    ext = DustMap.extinction(coords, filters=filters)

    assert len(ext) == len(coords)
    assert len(ext.colnames) == len(filters)


def test_plot_gas_map():
    plot = GasMap.plot_map(map_name="LAB")
    assert isinstance(plot, np.ma.core.MaskedArray)


def test_plot_dust_map():
    plot = DustMap.plot_map(map_name="SFD")
    assert isinstance(plot, np.ma.core.MaskedArray)
