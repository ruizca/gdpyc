gdpyc ─ Gas and Dust Python Calculator
======================================
.. inclusion-marker-main-readme

``gdpyc`` is a Python 2/3 package for calculating Hydrogen column density
and optical extinction. ``gdpyc`` offers functionalities similar
to the `nh tool`_ included in HEASoft, or on-line web services like
HEASARC's `nH`_ or IRSA's `Galactic Dust Reddening and Extinction`_.

This package uses HEALpix maps provided by NASA's LAMBDA service. Please
cite the original papers and authors of the surveys if you use this tool
for your research (see sections Surveys and References below).

|astropy| |DOI|

Dependencies
------------
``gdpyc`` depends on:

* ``astropy``
* ``astropy-healpix``
* ``numpy`` 

Certain functionalities also requiere:

* ``healpix``
* ``matplotlib``

Installation
------------

``gdpyc`` can be easily installed using ``pip``::

    pip install gdpyc

Example
-------
A simple example of using ``gdpyc``::

    >>> from gdpyc import GasMap, DustMap
    >>> from astropy.coordinates import SkyCoord
    
    >>> coords = SkyCoord(34.0, -5.0, unit='deg')
    >>> GasMap.nh(coords, nhmap='DL')

    <Quantity 1.9278499e+20 1 / cm2>

    >>> GasMap.nh(coords, nhmap='LAB')

    <Quantity 1.9802036e+20 1 / cm2>

    >>> DustMap.ebv(coords, dustmap='Planck13')

    '0.027179038520908336'

    >>> DustMap.extinction(coords, dustmap='SFD', filters='SDSS_r')

    <Table length=1>
           SDSS_r       
          float64       
    --------------------
    0.049108389441137615

    >>> GasMap.plot_map('HI4PI')

Surveys
-------
``gdpyc`` includes several HI and dust surveys with nH and E(B - V)
estimations. We created low resolution HEALPix maps for all surveys
(NSIDE=64 ~ 1 degree pixels) by degrading the original maps using
the `ud_grade` tool from `healpy`. Only low resolution maps are
included in the installation. If the user asks for high resolution
maps (`hires` parameter, see API documentation), they are downloaded
as needed and stored for future use.

HI surveys
^^^^^^^^^^

`DL`_: Composite all-sky map of neutral Hydrogen column density (NHI),
formed from the Leiden/Dwingeloo survey data [1]_ and the composite NHI
map of [2]_. The two datasets are not matched in sensitivity or resolution;
note that discontinuities exist in the constructed composite map. 

`DL high resolution data`_ (oversampled), NSIDE=512 ~ 0.11 deg.

`LAB`_: Observations of 21-cm emission from Galactic neutral Hydrogen
over the entire sky, merging the Leiden/Dwingeloo Survey [1]_ of the sky
north of -30° with the Instituto Argentino de Radioastronomia Survey
[3]_, [4]_ of the sky south of -25°. [5]_

`LAB high resolution data`_ (oversampled), NSIDE=512 ~ 0.11 deg. [6]_

`HI4PI`_: The HI 4-PI Survey (HI4PI) is a 21-cm all-sky survey of
neutral atomic Hydrogen. It is constructed from the Effelsberg-Bonn HI
Survey (EBHIS) and the Galactic All-Sky Survey (GASS). [7]_

`HI4PI high resolution data`_, NSIDE=1024 ~ 0.06 deg.

Dust surveys
^^^^^^^^^^^^
`SFD`_: All-sky map of Galactic reddening, E(B - V), from a
composite 100 micron map formed from IRAS/ISSA maps calibrated
using DIRBE observations. [8]_

`SFD high resolution data`_ (undersampled), NSIDE=512 ~ 0.11 deg.

`Planck13`_: All-sky map of Galactic reddening, E(B - V), using
Planck-HFI and IRAS data, for extra-galactic studies. [9]_

`Planck13 high resolution data`_, NSIDE=2048 ~ 0.03 deg.

Extinction values for different filters are estimated using the E(B - V)
conversion factors presented in [10]_, assuming an extinction to
reddening ratio R=3.1 Additional factors for 2MASS, `Spitzer`-IRAC
and WISE filters are from IRSA's `Galactic Dust Reddening and Extinction`_ 
service.

References
----------
.. [1] Hartmann & Burton 1997, Cambridge University Press.
.. [2] Dickey & Lockman 1990, Ann. Rev. A&A, 28, 215.
.. [3] Arnal et al. 2000, A&AS, 142.
.. [4] Bajaja et al. 2005, A&A, 440, 2
.. [5] Kalberla et al. 2005, A&A, 440, 775.
.. [6] Land & Slosar 2007, Phys. Rev. D, 76, 8.
.. [7] HI4PI Collaboration et al. 2016, A&A, 594, A116.
.. [8] Schlegel, Finkbeiner & Davis 1998, ApJ, 500, 2.
.. [9] Planck Collaboration et al. 2013, A&A, 571, A11.
.. [10] Schlafly & Finkbeiner 2011, ApJ, 737, 2, 103.


.. _nh tool: https://heasarc.gsfc.nasa.gov/lheasoft/ftools/heasarc.html
.. _nH: https://heasarc.gsfc.nasa.gov/cgi-bin/Tools/w3nh/w3nh.pl
.. _Galactic Dust Reddening and Extinction: https://irsa.ipac.caltech.edu/applications/DUST/
.. _DL: https://lambda.gsfc.nasa.gov/product/foreground/fg_combnh_map.cfm
.. _DL high resolution data: https://lambda.gsfc.nasa.gov/product/foreground/fg_HI_get.cfm
.. _LAB: https://lambda.gsfc.nasa.gov/product/foreground/fg_LAB_HI_Survey_info.cfm
.. _LAB high resolution data: https://lambda.gsfc.nasa.gov/product/foreground/fg_LAB_HI_Survey_get.cfm
.. _HI4PI: https://lambda.gsfc.nasa.gov/product/foreground/fg_hi4pi_info.cfm
.. _HI4PI high resolution data: https://lambda.gsfc.nasa.gov/product/foreground/fg_hi4pi_get.cfm
.. _SFD: https://lambda.gsfc.nasa.gov/product/foreground/fg_ebv_map.cfm
.. _SFD high resolution data: https://lambda.gsfc.nasa.gov/product/foreground/fg_sfd_get.cfm
.. _Planck13: https://wiki.cosmos.esa.int/planckpla/index.php/CMB_and_astrophysical_component_maps#Thermal_dust_emission
.. _Planck13 high resolution data: http://pla.esac.esa.int/pla/aio/product-action?MAP.MAP_ID=HFI_CompMap_ThermalDustModel_2048_R1.20.fits

.. |astropy| image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
   :target: http://www.astropy.org/

.. |DOI| image:: https://zenodo.org/badge/156710074.svg
   :target: https://zenodo.org/badge/latestdoi/156710074