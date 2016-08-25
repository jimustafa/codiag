from __future__ import absolute_import, division, print_function

from setuptools import find_packages
from numpy.distutils.core import setup, Extension
from numpy.distutils.system_info import get_info
from numpy.f2py.capi_maps import f2cmap_all


f2cmap_all['real'].update({
    'KIND(1.0D0)'.lower(): 'double'
    })

f2cmap_all['complex'].update({
    'KIND((1.0D0,1.0D0))'.lower(): 'complex_double'
    })


def extensions():
    lapack_opt = get_info('lapack_opt')

    if lapack_opt:
        ext1 = Extension(
            name='codiag.flib._flib',
            sources=[
                'f/nrtype.f90',
                'f/codiag.f90',
                'f/givens.f90',
                ],
            libraries=get_info('lapack_opt')['libraries'],
            library_dirs=get_info('lapack_opt')['library_dirs'],
            )

        extensions = [ext1]
    else:
        extensions = []

    return extensions


setup(
    name='codiag',
    version='0.1.0',
    description='Matrix codiagonalization library',
    author='Jamal I. Mustafa',
    author_email='jimustafa@gmail.com',
    license='BSD',
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
    packages=find_packages(exclude=['docs', 'tests']),
    setup_requires=[
        'numpy',
    ],
    install_requires=[
        'numpy',
        'scipy',
    ],
    tests_require=['pytest'],
    ext_modules=extensions(),
)
