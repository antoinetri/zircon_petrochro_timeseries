#!/usr/bin/env python
# -*- coding: utf-8 -*-
#############################################################################
# Performing bootstrap analysis on a compilation of zircon U-Pb + trace data
# Plotting boostrap results as a timeseries (ts)
# Author: Antoine Triantafyllou 2022-08 - Antoine.Triantafyllou@univ-lyon1.fr
# Associated paper: Triantafyllou et al. (2022) - Geology - Add title and vol.
#############################################################################


# Import useful libraries
import pandas as pd
import matplotlib.pyplot as plt
import bootstrapped.bootstrap as bs
import bootstrapped.stats_functions as bs_stats


def import_data(filename, sheet_name=0):
    """import xls file into dataframe"""
    data_in = pd.read_excel(filename, sheet_name=sheet_name)
    return data_in


def df_slice(start, end, col, df_in):
    """slicing in pandas df with conditions"""
    post_start = df_in[col] >= start
    pre_end = df_in[col] <= end
    df_sliced = df_in[post_start & pre_end]
    return df_sliced


def bootstrap_results(res_boot):
    """convert class object to actual values"""
    raw = res_boot.__str__().split(' ')
    median95 = float(raw[0])
    median95_high = float(raw[5].split(')')[0])
    temp2 = raw[4].split('(')[1]
    median95_low = float(temp2.split(',')[0])
    all_results = [median95, median95_low, median95_high]
    return all_results


def bootstrap_on_slice(df_slice, y_col, num_iteration=1000):
    """Bootstrapping in sliced portions of df"""
    df_no_val = df_slice.dropna(subset=[y_col])
    num_of_item = df_no_val.shape[0]
    if df_no_val.empty:
        med_distribution = [float("nan")]
        med68_bootstrap = [float("nan"), float("nan"), float("nan")]
        med95_bootstrap = [float("nan"), float("nan"), float("nan")]
        return med_distribution, med68_bootstrap, med95_bootstrap
    else:
        # run bootstrap for hist and values, gather median and percentile values (alphas) for each iteration
        med_distribution = bs.bootstrap(df_no_val[y_col].to_numpy(), stat_func=bs_stats.median,
                                        num_iterations=num_iteration,
                                        return_distribution=True)
        med95_res = bs.bootstrap(df_no_val[y_col].to_numpy(), stat_func=bs_stats.median, alpha=0.001,
                                 num_iterations=num_iteration,
                                 return_distribution=False)
        med95_bootstrap = bootstrap_results(med95_res)

        # for 1sig - 68% int confidence
        med68_res = bs.bootstrap(df_no_val[y_col].to_numpy(), stat_func=bs_stats.median, alpha=0.3173,
                                 num_iterations=num_iteration,
                                 return_distribution=False)
        med68_bootstrap = bootstrap_results(med68_res)
        return med_distribution, med68_bootstrap, med95_bootstrap


def plot_zircon(df_roll, colow, x_axis, y_axis):
    """plotting the boostrap result as a timeseries (rolling window)"""
    fig, ax1 = plt.subplots(figsize=(12, 8))
    ax1.set_xlabel('Time [Ma]')
    ax1.set_ylabel(y_axis)
    ax1.fill_between(df_roll[x_axis], y1=df_roll['ic_ic1_low'], y2=df_roll['ic_ic1_up'], facecolor=colow, alpha=0.50)
    ax1.fill_between(df_roll[x_axis], y1=df_roll['ic_ic2_low'], y2=df_roll['ic_ic2_up'], facecolor=colow, alpha=0.30)
    line_med, = ax1.plot(df_roll[x_axis], df_roll['ic_median'], '-', linewidth=2.5, color='black', alpha=0.9)
    print(df_roll['ic_median'].mean())
    print(df_roll['ic_median'].std())
    rel_std = 100 * df_roll['ic_median'].std() / df_roll['ic_median'].mean()
    print(rel_std)

    # Set limits of your timeseries plot
    ax1.set_xlim(4500, 0)
    ax1.set_ylim(0.0, 0.6)
    # ax1.set_yscale("log")

    # add supercontinents box
    supercontinents = False
    if supercontinents:
        nuna = plt.Rectangle((1750, -1000), 350, 10000, facecolor="black", alpha=0.1)
        ax1.add_patch(nuna)

        # TODO: kenor
        kenor = plt.Rectangle((2400, -1000), 350, 10000, facecolor="black", alpha=0.1)
        ax1.add_patch(kenor)

        rodinia = plt.Rectangle((950, -1000), 310, 10000, facecolor="black", alpha=0.1)
        ax1.add_patch(rodinia)

        # TODO: gondwana
        gondwana = plt.Rectangle((460, -1000), 200, 10000, facecolor="black", alpha=0.1)
        ax1.add_patch(gondwana)

        pangea = plt.Rectangle((160, -1000), 190, 10000, facecolor="black", alpha=0.1)
        ax1.add_patch(pangea)
    else:
        print('no supercontinent')
    plt.show()



###########################################################################
# this is where users should modify parameters for treatment and plotting
###########################################################################
if __name__ == '__main__':
    # Add the directory of the zircon U-Pb-trace file
    filename = r'your\directory\Detrital_zircons_UPb_trace.xlsx'

    # Specify the sheet name in your xlsx file or set 0 if it's the first one
    sheet_name = 0

    # import all data and put it in a pd df
    df_bin = import_data(filename, sheet_name)

    # specify headers, for x put the zircon age and for y, the trace element to be bootstrapped and plotted as a ts (+ color)
    x_ax = 'Age'
    y_ax = 'Eu/Eu_zr_n'
    colow = 'orange'

    # specify the size of rolling window (bin_size), the time span to roll in (total_time) and number of iteration for bootstrap analysis (it_boot)
    bin_size = 500
    total_time = 4500
    it_boot = 500
    roll_step = 20

    # some list to concatenate then
    time_pos, med_bootstrap, ic68_btsp_low, ic68_btsp_high, ic95_btsp_low, ic95_btsp_high, num_of_item = [], [], [], [], [], [], []

    # this "for" loop iterates each time window and roll and apply bootstrap analysis on each time window (20 Myr roll here)
    for roll_step in range(-250, total_time, roll_step):
        start_lim = roll_step + 1
        print(str(start_lim) + ' Ma')
        end_lim = start_lim + bin_size
        position = start_lim + (bin_size / 2)
        df_in_bin = df_slice(start_lim, end_lim, x_ax, df_bin)
        med_distribution, med_results_68, med_results_95 = bootstrap_on_slice(df_in_bin, y_ax, it_boot)
        med_bootstrap.append(med_results_68[0])
        ic68_btsp_low.append(med_results_68[1])
        ic68_btsp_high.append(med_results_68[2])
        ic95_btsp_low.append(med_results_95[1])
        ic95_btsp_high.append(med_results_95[2])
        time_pos.append(position)

    # creating a "clean" df fullfilled with lists of bootstrap result
    df_out = pd.DataFrame(columns=['Age', 'ic_median', 'ic_ic1_low', 'ic_ic1_up', 'ic_ic2_low', 'ic_ic2_up', 'nb_folks'])
    df_out['Age'] = time_pos
    df_out['ic_median'] = med_bootstrap
    df_out['ic_ic1_low'] = ic68_btsp_low
    df_out['ic_ic1_up'] = ic68_btsp_high
    df_out['ic_ic2_low'] = ic95_btsp_low
    df_out['ic_ic2_up'] = ic95_btsp_high

    # plot the timeseries
    plot_zircon(df_out, colow, x_ax, y_ax)

