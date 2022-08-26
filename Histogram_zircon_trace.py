#!/usr/bin/env python
# -*- coding: utf-8 -*-
#############################################################################
# Plotting an histogram on a compilation of zircon U-Pb + trace data
# Author: Antoine Triantafyllou 2022-08 - Antoine.Triantafyllou@univ-lyon1.fr
# Associated paper: Triantafyllou et al. (2022) - Geology - Add title and vol.
#############################################################################


# Import useful libraries
import pandas as pd
import matplotlib.pyplot as plt


def import_data(filename, sheet_name=0):
    """import xls file into dataframe"""
    data_in = pd.read_excel(filename, sheet_name=sheet_name)
    return data_in


###########################################################################
# this is where users should modify parameters for treatment and plotting
###########################################################################
if __name__ == '__main__':
    # Add the directory of the zircon U-Pb-trace file
    filename = r'your\directory\Detrital_zircons_UPb_trace.xlsx'

    # import all data and put it in a pd df
    df_in = import_data(filename)

    # Specify the sheet name in your xlsx file or set 0 if it's the first one
    sheet_name = 0

    # specify headers, for x put the zircon age and for y (named fact here), the trace element to be plotted as a histogram (+ bins in hist)
    fact = 'Eu/Eu_zr_n'
    x_ax = 'Age'
    number_of_bins = 160

    # clean the datafram from nan data
    df_no_fact = df_in.dropna(subset=[fact])
    df_no_age = df_no_fact.dropna(subset=[x_ax])

    # plotting parameters (limits, labels, color)
    fig, ax1 = plt.subplots(figsize=(12, 8))
    ax1.set_xlim(4250, 0)
    ax1.set_xlabel('Zircon ages bins' + '\n' + '(' + str(fact) + '_dropped - ' + str(x_ax) + '_dropped)')
    ax1.set_ylabel('Relative frequency')
    ax1.hist(df_no_age[x_ax], bins=number_of_bins, color='grey', alpha=0.3)
    plt.show()

