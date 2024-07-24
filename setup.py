from setuptools import setup, find_packages

setup(
    name="lipss_sipe",
    version="0.58",
    packages=find_packages(),
    # required packages
    install_requires=["numpy", "matplotlib", "tqdm"],
    author="Thies Johannes Albert, Dominik Kaczmarek",
    author_email="thies.albert@uni-due.de",
    description="Python package and Mathematica Notebook to calculate the efficacy factor of LIPSS formation by SIPE",
    license="MIT",
    keywords="lipss sipe efficacy factor",
    url="https://github.com/tjalb/lipss_sipe",
)
