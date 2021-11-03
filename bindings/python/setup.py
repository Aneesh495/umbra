"""
setup.py for umbra Python package.
"""

from setuptools import setup, find_packages

setup(
    name="umbra",
    version="1.0.0",
    description="High-Order Discontinuous Galerkin GRMHD Kernel and Kerr Radiative Transfer",
    author="umbra development team",
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.20.0",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Fortran",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Astronomy",
    ],
)
