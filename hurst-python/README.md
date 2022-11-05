# Steps to run

In this file, the systems requirements, instalation guide and instructions for reproduction of results are presented.

## Systems Requirements

This software was tested on Windows 10, but once all requirements are fullfilled, one should be able to test it on either Linux or Mac too.

Before running, we advice to have the [most recent Python version](https://www.python.org/downloads/release/python-3110/) -- tested on v. 3.11.0 and olders.

## Instalation Guide

To run the software, one must install the Python Packages listed bellow. Make sure that all the following libraries are installed and updated:
 
- `os`
- `csv`
- `numpy`
- `sklearn`
- `math`
- `matplotlib`
- `pyplot`
- `scipy`

Follow the [Python Packaging User Guide](https://packaging.python.org/en/latest/tutorials/installing-packages/) if you have any missing package. The `pip` package manager manager is recommended, and one can insatll the packages by

`pip install package_name`

See the [pip docs](https://pip.pypa.io/en/latest/) for more details.

The install time will generally depend on your internet conection and computer hardware, but should not exceed 5 to 10 minutes depending on how many of the above packages must be installed manually.

## Demo

The hurst exponent is calculated by running  `python Hurst.py` in a terminal opened at this folder.

This will result in the cration of an output file `RESULTS.csv` at this same folder with the hurst exponent for each individual sporozoite.

## Reproducing figures

From the CSV file with individual hurst exponents, the mean hurst exponent for INV or NINV sporozoites can be calculated by selecting only those labeled as INV or NINV on first column.
This will reproduce figure 3a.

The expecting result is foun in file `EXPECTED_RESULTS.csv`
