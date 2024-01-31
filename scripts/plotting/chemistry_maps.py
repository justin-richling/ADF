"""
Module: regional_map_multicase

Provides a plot with regional maps of specified variables for all cases (up to 10) in a row.

Since this is a specialty plot, it looks for several custom options to be provided in the YAML configuration file. For example, one might include this block in the YAML:

region_multicase:
    region_spec: [slat, nlat, wlon, elon]
    region_time_option: <calendar | zeroanchor>  If calendar, will look for specified years. If zeroanchor will use a nyears starting from year_offset from the beginning of timeseries
    region_start_year:
    region_end_year:
    region_nyear:
    region_year_offset:
    region_month: <NULL means look for season>
    region_season: <NULL means use annual mean>
    region_variables: <list of variables to try to use; allows for a subset of the total diag variables>

"""
#
# --- imports and configuration ---
#
from pathlib import Path
import warnings  # use to warn user about missing files.

import numpy as np
import xarray as xr
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

from plotting_functions import pres_from_hybrid, prep_contour_plot
import plotting_functions as pf


def my_formatwarning(msg, *args, **kwargs):
    # ignore everything except the message
    return str(msg) + "\n"


warnings.formatwarning = my_formatwarning

#Set seasonal ranges:
seasons = {"ANN": np.arange(1,13,1),
            "DJF": [12, 1, 2],
            "JJA": [6, 7, 8],
            "MAM": [3, 4, 5],
            "SON": [9, 10, 11]}

#Set aerosol variables with constituents
aerosol_dict = {"BC":["bc_a1", "bc_a4"],
                "POM":["pom_a1", "pom_a4"],
                "SO4":["so4_a1", "so4_a2",  "so4_a3", "so4_a5"],
                "SOA":["soa_a1", "soa_a2"],
                "DUST":["dst_a1", "dst_a2", "dst_a3"],
                "SeaSalt":["ncl_a1", "ncl_a2", "ncl_a3"]}

#
# --- Main Function Shares Name with Module: regional_map_multicase ---
#
def make_chem_maps(adfobj, diag, data_dict, case_deets):
    """
    Chemistry Map main function
        * Initially start with Aerosol Zonal maps
            - This probably can be expanded to LatLon if given single pressure levels?

        data_dict: data_dict[case_name][s][f"m{var}"]
    """

    # Notify user that script has started:
    print("\n  Generating zonal aerosol plots ...")

    #for i in data_dict:
    #    for j in data_dict[i]:
    #        for k in data_dict[i][j]:
    #            print(k,"\n\n")

    #Set plot file type:
    # -- this should be set in basic_info_dict, but is not required
    # -- So check for it, and default to png
    basic_info_dict = adfobj.read_config_var("diag_basic_info")
    plot_type = basic_info_dict.get('plot_type', 'png')

    var_list = adfobj.diag_var_list

    #Aerosol Calculations
    if diag == "aerosol":

        for var in aerosol_dict:
            #if all(elem in var for elem in constits):
            #    print(f"\t - zonal mean aerosol plots for {var}")
                        
            #If found then notify user, assuming debug log is enabled:
            adfobj.debug_log(f"zonal_mean: Found variable defaults for {var}")
            aerosol_plot(adfobj, var, data_dict, case_deets)

            #else:
            #    print(f"No constituents for {var}, moving on ...")

    
# Specific Plot Functions
#------------------------
# 

# Aerosols
def aerosol_plot(adfobj, var, data_dict, case_deets):
    """
    This assumes all composite variables are a linear combination of constituents
    """

    #Set plot file type:
    # -- this should be set in basic_info_dict, but is not required
    # -- So check for it, and default to png
    basic_info_dict = adfobj.read_config_var("diag_basic_info")
    plot_type = basic_info_dict.get('plot_type', 'png')
    redo_plot = adfobj.get_basic_info("redo_plot")

    #Set category
    web_category = "Aerosols"

    case_names = case_deets["case_names"]["cases"]

    #Loop over model cases:
    for case_idx, case_name in enumerate(case_names):
        for s in seasons:
            maerosol = 0
            oaerosol = 0
            for j in aerosol_dict[var]:
                print("blah",j)
                if (j not in data_dict[case_name][s]) or (j not in data_dict[case_deets["case_names"]["baseline"]][s]):
                    print(f"missing constituent for {var}, moving on to the next aerosol variable")
                    continue
                maerosol += data_dict[case_name][s][f"{j}"]["mdata"]
                oaerosol += data_dict[case_deets["case_names"]["baseline"]][s][f"{j}"]["odata"]

            #Gather info from data and case details
            case_nickname = case_deets["nicknames"]["cases"][case_idx]
            base_nickname = case_deets["nicknames"]["baseline"]

            case_years = [case_deets["years"]["syears"][case_idx],case_deets["years"]["eyears"][case_idx]]
            baseline_years = [case_deets["years"]["syear_baseline"],case_deets["years"]["eyear_baseline"]]

            has_lev = data_dict[case_name][s][f"m{j}"]["has_lev"]

            plot_loc = data_dict[case_name][s][f"m{j}"]["plot_loc"]
            ####

            plot_name = plot_loc / f"{var}_{s}_Zonal_Mean.{case_deets['ptype']}"
            # Check redo_plot. If set to True: remove old plot, if it already exists:
            redo_plot = adfobj.get_basic_info('redo_plot')
            if (not redo_plot) and plot_name.is_file():
                #Add already-existing plot to website (if enabled):
                adfobj.add_website_data(plot_name, var, case_name, category=web_category,
                                                season=s, plot_type="Zonal")

                #Continue to next iteration:
                #continue
                pass
            elif (redo_plot) and plot_name.is_file():
                plot_name.unlink()

            pf.plot_zonal_mean_and_save(plot_name, case_nickname, base_nickname,
                                        case_years,
                                        baseline_years,
                                        maerosol, oaerosol, has_lev,log_p=False,
                                        **case_deets["vres"])

            #Add plot to website (if enabled):
            adfobj.add_website_data(plot_name, var, case_name, category=web_category,
                                                season=s, plot_type="Zonal")

            #Create new plot with log-p:
            if has_lev:
                plot_name_log = plot_loc / f"{var}_{s}_Zonal_logp_Mean.{plot_type}"

                 # Check redo_plot. If set to True: remove old plot, if it already exists:
                if (not redo_plot) and plot_name_log.is_file():
                    #Continue to next iteration:
                    continue

                elif (redo_plot) and plot_name_log.is_file():
                    plot_name_log.unlink()
                #End if

                pf.plot_zonal_mean_and_save(plot_name_log, case_nickname, base_nickname,
                                                case_years,
                                                baseline_years,
                                                maerosol, oaerosol, has_lev, log_p=True, **case_deets["vres"])

                #Add plot to website (if enabled):
                adfobj.add_website_data(plot_name_log, f"{var}_logp", case_name, season=s, plot_type="Zonal", category="Log-P")

####



##############
#END OF SCRIPT
