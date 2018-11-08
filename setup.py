from setuptools import setup

with open("README.rst", "r") as fh:
    long_description = fh.read()

setup(
    name='gdpyc',
    version='1.0',
    author='Angel Ruiz',
    author_email='angel.ruizca@gmail.com',
    description='Gas and Dust Python Calculator',
    long_description=long_description,
    long_description_content_type='text/x-rst',
    url='https://github.com/ruizca/gdpyc',
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    install_requires=['astropy', 'numpy', 'astropy-healpix', 'regions'],
    extras_require={
        'plot_map': ['matplotlib', 'healpy']
    },
    packages=['gdpyc'],
    package_dir={'gdpyc': 'gdpyc'},
    package_data={'gdpyc': ['data/*_lowres.fits', 
                            'data/*_heasoft.fits',
                            'data/sfcoeff.fits']},
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.6',
        'License :: OSI Approved :: BSD License',
        'Operating System :: POSIX :: Linux',
    ],
)
