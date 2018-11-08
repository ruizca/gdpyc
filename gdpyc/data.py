"""
Gas and Dust Python Calculator.

Module for downloading and processing
the HEALpix maps used in gdpyc.
"""
import os
import pkg_resources

from astropy.io import fits
from astropy.utils.data import get_readable_fileobj


def download_map(url, data_dir, filename=None):
    if filename is None:
        filename = os.path.basename(url)
    filepath = os.path.join(data_dir, filename)

    with get_readable_fileobj(url, encoding='binary') as fs, open(filepath, 'wb') as fd:
        fd.write(fs.read())

    return filepath


def parse_DL(mapfile, outfile):
    with fits.open(mapfile) as hdu:
        header = hdu[1].header
        nhmap = hdu[1].data.field(0)

    # Flag negative values as UNSEEN in healpix convention
    nhmap[nhmap < 0] = -1.6375e+30
    header['INDXSCHM'] = 'IMPLICIT'

    col1 = fits.Column(name='TEMPERATURE', format='D', array=nhmap)
    chdu = fits.BinTableHDU.from_columns(fits.ColDefs([col1]),
                                         header=header, name='LAB')
    chdu.writeto(outfile, overwrite=True)


def parse_LAB(mapfile, outfile):
    # The original file contains a temperature map,
    # it has to be transformed to nH.
    with fits.open(mapfile) as hdu:
        header = hdu[1].header
        tmap = hdu[1].data.field(0)
        nchannels = hdu[1].data.field(1)

    # Optically thin approximation from Eq.3, Dickey & Lockman 1990.
    # The integral is approximated by <T>*nchannels*dv
    # tmap and nchannels are both nrows x 1024 arrays
    # 1.823e18
    nhmap = 1.8224e18 * tmap * nchannels * 1.030571969

    nrows, ncols = nhmap.shape
    col1 = fits.Column(name='TEMPERATURE', format='D',
                       array=nhmap.reshape(nrows*ncols))
    chdu = fits.BinTableHDU.from_columns(fits.ColDefs([col1]),
                                         header=header, name='LAB')
    chdu.writeto(outfile, overwrite=True)


def parse_HI4PI(mapfile, outfile):
    # The original fits file has 5 columns (HPXINDEX, RA2000,
    # DEC2000, GLON, GLAT, NHI). Keep only NHI.
    with fits.open(mapfile) as hdu:
        header = hdu[1].header
        nhmap = hdu[1].data.field(5)

    header['INDXSCHM'] = 'IMPLICIT'
    col1 = fits.Column(name='TEMPERATURE', format='D', array=nhmap)
    chdu = fits.BinTableHDU.from_columns(fits.ColDefs([col1]),
                                         header=header, name='HI4PI')
    chdu.writeto(outfile, overwrite=True)


def parse_Planck13(mapfile, outfile):
    # The original fits file has 8 columns. Keep only EBV.
    with fits.open(mapfile) as hdu:
        header = hdu[1].header
        dustmap = hdu[1].data.field('EBV')

    col1 = fits.Column(name='TEMPERATURE', format='D', array=dustmap)
    chdu = fits.BinTableHDU.from_columns(fits.ColDefs([col1]),
                                         header=header, name='Planck13')
    chdu.writeto(outfile, overwrite=True)


def get_map(maplabel, data_dir):

    maps = {'DL': 'https://lambda.gsfc.nasa.gov/data/foregrounds/combined_nh/lambda_combined_nh.fits',
            'LAB': 'https://lambda.gsfc.nasa.gov/data/foregrounds/HI/LAB_fullvel.fits',
            'HI4PI': 'https://lambda.gsfc.nasa.gov/data/foregrounds/HI4PI/NHI_HPX.fits',
            'SFD': 'https://lambda.gsfc.nasa.gov/data/foregrounds/SFD/lambda_sfd_ebv.fits',
            'Planck13': 'http://pla.esac.esa.int/pla/aio/product-action?MAP.MAP_ID=HFI_CompMap_ThermalDustModel_2048_R1.20.fits'
           }

    filename = '{}_{}_healpix_hires.fits'

    if maplabel == 'DL':
        mapfile = download_map(maps[maplabel], data_dir)
        newmapfile = os.path.join(data_dir, filename.format('h1_nh', maplabel))
        parse_DL(mapfile, newmapfile)
        os.remove(mapfile)

    elif maplabel == 'LAB':
        mapfile = download_map(maps[maplabel], data_dir)
        newmapfile = os.path.join(data_dir, filename.format('h1_nh', maplabel))
        parse_LAB(mapfile, newmapfile)
        os.remove(mapfile)

    elif maplabel == 'HI4PI':
        mapfile = download_map(maps[maplabel], data_dir)
        newmapfile = os.path.join(data_dir, filename.format('h1_nh', maplabel))
        parse_HI4PI(mapfile, newmapfile)
        os.remove(mapfile)

    elif maplabel == 'SFD':
        mapfile = filename.format('dust_ebv', maplabel)
        download_map(maps[maplabel], data_dir, mapfile)

    elif maplabel == 'Planck13':
        mapfile = download_map(maps[maplabel], data_dir)
        newmapfile = os.path.join(data_dir, filename.format('dust_ebv', maplabel))
        parse_Planck13(mapfile, newmapfile)
        os.remove(mapfile)

    else:
        raise ValueError('Unknown map: {}'.format(maplabel))


def degrade_map(mapfile, order=6):
    import healpy as hp

    newmapfile = mapfile.replace('_hires.fits', '_lowres.fits')

    hpmap = hp.read_map(mapfile)
    hpmap_low = hp.ud_grade(hpmap, 2**order)
    hp.write_map(newmapfile, hpmap_low, fits_IDL=False, overwrite=True)


def get_lowres_maps(data_dir):
    maps = ['DL', 'LAB', 'HI4PI', 'SFD', 'Planck13']

    mapfile = os.path.join(data_dir, '{}_{}_healpix_hires.fits')
    for key in maps:
        if key == 'SFD' or key == 'Planck13':
            degrade_map(mapfile.format('dust_ebv', key))
        else:
            degrade_map(mapfile.format('h1_nh', key))


def main():
    data_dir = pkg_resources.resource_filename('gdpyc', 'data')

    #get_map('HI4PI', data_dir)
    #get_map('DL', data_dir)
    #get_map('LAB', data_dir)
    #get_map('SFD', data_dir)
    #get_map('Planck13', data_dir)
    get_lowres_maps(data_dir)

if __name__ == '__main__':
    main()
