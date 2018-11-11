"""
Gas and Dust Python Calculator
==============================

A python module for estimating extinction,
reddening and Hydrogen column densities.
"""
from __future__ import print_function

import os
import warnings
import pkg_resources

from astropy import units as u
from astropy.io import fits
from astropy.coordinates import Galactic
from astropy.table import Table, join
from astropy.wcs import WCS
from astropy_healpix import HEALPix
from regions import PixCoord
import numpy as np


class Map(object):
    """
    A class with common methods for HEALpix maps.
    """

    _data_path = pkg_resources.resource_filename('gdpyc', 'data')
    _map_type = None
    _maps = None

    @classmethod
    def show_maps(cls):
        """
        Show available maps.
        """
        print('Available maps:')
        for key in cls._maps.keys():
            print('- {} ({})'.format(key, cls._maps[key]))


    @classmethod
    def plot_map(cls, map_name, plotname=None):
        """
        Full-sky plot (mollweide projection) of the map.

        **Note:** this method needs ``healpy`` and ``matplotlib``.

        Parameters
        ----------
        map_name : ``str``
            Name of the map to be plotted. Use ``show_maps`` method
            to see a list of all available maps.
        plotname : ``str`` or ``None``, optional
            Name of the file where the plot will be saved. If ``None``,
            the plot is shown but not saved. Defaults to ``None``.

        Returns
        -------
        plot : ``numpy.ndarray``
            2D numpy array with the plot.
        """
        import healpy as hp
        import matplotlib.pyplot as plt

        cls._check_map(map_name)

        hpmapfile = '{}_{}_healpix_lowres.fits'.format(cls._map_type, map_name)
        hpmapfile = os.path.join(cls._data_path, hpmapfile)
        hpmap = hp.read_map(hpmapfile)

        title = '{} ({})'.format(map_name, cls._maps[map_name])
        if cls. _map_type == 'h1_nh':
            unit_label = 'cm-2'
            minval, maxval = 1e19, 3e22
        else:
            unit_label = 'mag'
            minval, maxval = 0, hpmap.max()

        plot = hp.mollview(hpmap, title=title, norm='hist',
                           min=minval, max=maxval, unit=unit_label,
                           return_projected_map=True)

        if plotname is None:
            plt.show()
        else:
            plt.savefig(plotname)

        return plot


    @classmethod
    def _load_map(cls, map_name, hires=False):
        # Load healpix map for map_name.
        # If hires is True and the map is not locally available,
        # it is downloaded.
        cls._check_map(map_name)

        if hires:
            resolution = 'hires'
        else:
            resolution = 'lowres'

        hmapfile = '{}_{}_healpix_{}.fits'.format(cls._map_type, map_name, resolution)
        hmapfile = os.path.join(cls._data_path, hmapfile)

        if hires and not os.path.isfile(hmapfile):
            print('High resolution map is not locally available.')
            print('Downloading {} map...'.format(map_name))
            cls._download_map(map_name)

        with fits.open(hmapfile) as hdu:
            nside = hdu[1].header['NSIDE']
            order = hdu[1].header['ORDERING']
            hpmap = hdu[1].data.field(0)

        return hpmap, nside, order


    @classmethod
    def _download_map(cls, map_name):
        from .data import get_map
        get_map(map_name, cls._data_path)


    @classmethod
    def _interpolate(cls, coords, hpmap, nside, order):
        hp = HEALPix(nside=nside, order=order, frame=Galactic())

        return hp.interpolate_bilinear_skycoord(coords, hpmap)

    @classmethod
    def _check_map(cls, map_name):
        if map_name not in cls._maps.keys():
            message = 'Map {} unknown!\nAvailable maps: {}'
            raise ValueError(message.format(map_name, cls._maps.keys()))


class GasMap(Map):
    """
    Class for HI maps.
    """

    _map_type = 'h1_nh'
    _maps = {'DL':  'Dickey & Lockman 1990, Ann. Rev. A&A, 28, 215',
             'LAB': 'Kalberla et al. 2005, A&A, 440, 775',
             'HI4PI': 'HI4PI Collaboration et al. 2016, A&A, 594, A116'}

    @classmethod
    def nh(cls, coords, nhmap='LAB', hires=False):
        """
        Hydrogen column density in the line-of-sight of `coords`,
        using LAMBDA healpix maps [1]_.

        Parameters
        ----------
        coords : ``SkyCoord``
            Astropy SkyCoord object with the line-of-sight coordinates.
        nhmap : ``str``, optional
            Name of the HI survey. Use ``show_maps`` method
            to see a list of all available maps. Defaults to 'LAB'.
        hires : ``boolean``, optional
            Use high resolution map. If the map is not available locally,
            it is downloaded. Defaults to ``False``.

        Returns
        -------
        nh : ``Quantity``
            An Astropy Quantity object with shape like `coords`, in cm**-2.

        References
        ----------
        .. [1] https://lambda.gsfc.nasa.gov/product/foreground/fg_diffuse.cfm
        """
        nh_hpmap, nside, order = cls._load_map(nhmap, hires=hires)
        nh = cls._interpolate(coords, nh_hpmap, nside, order)

        return nh * u.cm**-2


    @classmethod
    def nhtotal(cls, coords, hires=False):
        """
        Total Hydrogen column density (NHI + 2NH2) in the line-of-sight
        of `coords`, using Willingale's method [2]_.

        Parameters
        ----------
        coords : ``SkyCoord``
            Astropy SkyCoord object with the line-of-sight coordinates.
        hires : ``boolean``, optional
            Use high resolution maps. If the map is not available locally,
            it is downloaded. Defaults to ``False``.

        Returns
        -------
        nhtotal : ``Quantity``
            An Astropy Quantity object with shape like `coords`, in cm**-2.

        References
        ----------
        .. [2] Willingale et al. 2013, MNRAS, 431, 1.
        """
        nHI = cls.nh(coords, nhmap='LAB', hires=hires)
        ebv = DustMap.ebv(coords, dustmap='SFD', hires=hires)

        nH2max = 7.3e20 * u.cm**-2
        Nc = 3.0e20 * u.cm**-2
        alpha = 1.1
        nH2 = nH2max*(1 - np.exp(-nHI*ebv/Nc))**alpha

        return nHI + 2*nH2


    @classmethod
    def nhf(cls, coords, nhmap='LAB', radius=1.0*u.deg):
        """
        Hydrogen column density in the line-of-sight of `coords`,
        using HEASoft fits images (resolution of 0.675 x 0.675 deg) [3]_
        and method.

        Parameters
        ----------
        coords : ``SkyCoord``
            Astropy SkyCoord object with the line-of-sight coordinates.
        nhmap : ``str``, optional
            Name of the HI survey. Use ``show_maps`` method
            to see a list of all available maps. Defaults to 'LAB'.
        radius : ``Quantity``, optional
            Radius of the circle around `coords` where the nH value is averaged.
            An Astropy angular Quantity object, consistent with degrees.

        Returns
        -------
        nh : ``Quantity``
            An Astropy Quantity object with shape like `coords`, in cm**-2..

        References
        ----------
        .. [3] Original fits files created by K. Kuntz (LAB) and
               Steve Snowden (DL).
        """
        if nhmap not in ['LAB', 'DL']:
            raise ValueError('Only LAB and DL maps are available')

        radius = radius.to(u.deg).value
        if radius > 3:
            raise ValueError('radius must be <= 3 deg!!!')

        hmapfile = '{}_{}_heasoft.fits'.format(cls._map_type, nhmap)
        hmapfile = os.path.join(cls._data_path, hmapfile)
        hmapimage = fits.getdata(hmapfile)
        wcs = WCS(hmapfile)

        # If coords is not an array of coordinates, change to a list
        try:
            len(coords)
        except TypeError:
            coords = np.array([coords])

        # TODO: vectorize this loop
        # numpy iterator for preserving the shape of coords in nh
        nh = np.empty_like(coords)
        it = np.nditer(coords, flags=['multi_index', 'refs_ok'])
        while not it.finished:
            nh[it.multi_index] = cls._nhftools(coords[it.multi_index],
                                               hmapimage, wcs, radius)
            it.iternext()

        if len(nh) == 1:
            nh = float(nh)

        return nh * u.cm**-2


    @classmethod
    def _nhftools(cls, center_sky, image, wcs, radius=1.0):
        ### Quick and dirty python
        ### implementation of ftool's nh
        # Select a 3 x 3 deg subimage (default in ftool's nh)
        box_pixcoord = cls._subimage(center_sky, wcs, size=3.0)

        # Select only pixels with non-zero values to avoid border effects
        out_mask = image[box_pixcoord.y, box_pixcoord.x] > 0

        # Estimate sky distance of each pixel (in deg) to center
        box_skycoord = box_pixcoord.to_sky(wcs)
        distance = np.empty_like(box_skycoord)
        distance[out_mask] = box_skycoord[out_mask].separation(center_sky)
        distance[~out_mask] = -99

        # Select pixels within radius
        good_mask = np.logical_and(distance <= radius, distance > 0)
        good_nh = image[box_pixcoord[good_mask].y, box_pixcoord[good_mask].x]

        # Estimate weights
        weights = (radius - distance[good_mask])/radius

        if weights:
            nh = np.sum(good_nh * weights)/weights.sum()
        else:
            message = ('No points are within {} deg from input position. '
                       'First good point is at distance {} deg')
            message = message.format(radius, distance[out_mask].min())
            warnings.warn(RuntimeWarning(message))

            nh = np.nan

        return nh


    @staticmethod
    def _subimage(center_sky, wcs, size=3.0):
        center_x, center_y = center_sky.to_pixel(wcs)
        npixtot = size/wcs.wcs.cdelt[1]

        # Box limits
        fpix_x = int(np.round(center_x - (npixtot - 1)/2.0))
        fpix_y = int(np.round(center_y - (npixtot - 1)/2.0))
        lpix_x = int(np.round(center_x + (npixtot - 1)/2.0))
        lpix_y = int(np.round(center_y + (npixtot - 1)/2.0))

        xbox = np.arange(fpix_x, lpix_x + 1)
        ybox = np.arange(fpix_y, lpix_y + 1)

        # Define box within the pixel limits of the fits image
        mask_x = np.logical_and(xbox >= 0, xbox < wcs._naxis1)
        mask_y = np.logical_and(ybox >= 0, ybox < wcs._naxis2)

        # Grid of pixel coordinates for the subimage
        grid_x, grid_y = np.meshgrid(xbox[mask_x], ybox[mask_y])

        return PixCoord(x=grid_x, y=grid_y)


class DustMap(Map):
    """
    Class for dust maps.
    """

    _map_type = 'dust_ebv'
    _maps = {'SFD':  'Schlegel, Finkbeiner & Davis 1998, ApJ, 500, 2',
             'Planck13': 'Planck Collaboration et al. 2013, A&A, 571, A11'}

    @classmethod
    def ebv(cls, coords, dustmap='SFD', hires=False):
        """
        E(B - V) in the line-of-sight of `coords`.

        Parameters
        ----------
        coords : ``SkyCoord``
            Astropy SkyCoord object with the line-of-sight coordinates.
        nhmap : ``str``, optional
            Name of the dust survey. Use ``show_maps`` method
            to see a list of all available maps. Defaults to 'SFD'.
        hires : ``boolean``, optional
            Use high resolution map. If the map is not available locally,
            it is downloaded. Defaults to `False`.

        Returns
        -------
        ebv : ``float`` or ``numpy.ndarray``
            Float, or numpy array with shape like `coords`.
        """
        ebv_hpmap, nside, order = cls._load_map(dustmap, hires=hires)
        ebv = cls._interpolate(coords, ebv_hpmap, nside, order)

        return ebv


    @classmethod
    def extinction(cls, coords, dustmap='SFD',
                   filters='default', hires=False):
        """
        Galactic extinction in the line-of-sight of `coords` for the
        bandpasses defined in `filters`. If there are N coords and
        M filters, returns an astropy table of shape N x M.

        **Note:** E(B - V) to extinction conversion coefficients were
        estimated in [4]_ for the SFD map [5]_. Using this method with a
        different dust map is possible, but obtaining accurate values is
        highly unlikely.

        Parameters
        ----------
        coords : ``SkyCoord``
            Astropy SkyCoord object with the line-of-sight coordinates.
        nhmap : ``str``, optional
            Name of the dust survey. Use ``show_maps`` method
            to see a list of all available maps. Defaults to 'SFD'.
        filters : ``str`` or ``list``, optional
            List of the filters (or, equivalently, a single string with
            comma separated filters). Use ``show_filters`` method
            to see a list of all available maps. If filters is 'default',
            returns extinction for the five SDSS bands (u, g, r, i, z).
            Defaults to 'default'.
        hires : ``boolean``, optional
            Use high resolution map. If the map is not available locally,
            it is downloaded. Defaults to ``False``.

        Returns
        -------
        ext : ``Table``
            Astropy Table with the extinction values at `coords`
            for each selected filter.

        References
        ----------
        .. [4] Schlafly & Finkbeiner 2011, ApJ, 737, 2, 103.
        .. [5] Schlegel, Finkbeiner & Davis 1998, ApJ, 500, 2.
        """
        list_filters = cls._parse_filters(filters)
        ebv = cls.ebv(coords, dustmap=dustmap, hires=hires)
        aebv = cls._ebv_to_ext(list_filters)

        if dustmap != 'SFD':
            warnings.warn('Extinction for a dust map other than SFD.\n'
                          'Results are not reliable!!!')

        if isinstance(ebv, np.ndarray):
            ext = np.matmul(ebv[:, np.newaxis], aebv[np.newaxis, :])
        else:
            ext = ebv * aebv

        return Table(ext, names=list_filters)


    @classmethod
    def show_filters(cls, lambda_eff=False):
        """
        List of the available filters in the S&F conversion table [4]_.

        Parameters
        ----------
        lambda_eff : ``boolean``, optional
            If ``True``, also returns the effective wavelength of each filter.
            Deafults to ``False``.
        """
        aebv = cls._load_sfcoeff()

        if lambda_eff:
            return aebv['filter'], aebv['lambda_eff']
        else:
            return aebv['filter']


    @classmethod
    def _ebv_to_ext(cls, filters):
        # Conversion coefficients from E(B - V) to extinction in `filters`.
        good_filters = cls.show_filters()
        for f in filters:
            if f not in good_filters:
                raise ValueError('Unknown filter: {}'.format(f))

        aebv = cls._load_sfcoeff()

        sortidx = list(range(len(filters)))
        filters = Table([filters, sortidx], names=['filter', 'sortidx'])

        sfcoeff_filters = join(filters, aebv, keys=['filter'], join_type='left')
        sfcoeff_filters.sort('sortidx')

        return np.array(sfcoeff_filters['AEBV2']) # AEBV2 assumes RV = 3.1


    @classmethod
    def _load_sfcoeff(cls):
        # Load conversion coefficients from the SFD maps of E(B - V) to
        # extinction in several bandpasses, by Schlafly and Finkbeiner (2011).
        sfcoeff = os.path.join(cls._data_path, 'sfcoeff.fits')
        return Table.read(sfcoeff, format='fits')


    @classmethod
    def _parse_filters(cls, filters):
        if isinstance(filters, str):
            if filters == 'default':
                filters = 'SDSS_u, SDSS_g, SDSS_r, SDSS_i, SDSS_z'
            list_filters = filters.replace(' ', '').split(',')
        else:
            list_filters = filters

        return list_filters
