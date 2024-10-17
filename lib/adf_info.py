"""
Information/Parameter (Info) class for
the Atmospheric Diagnostics Framework (ADF).

This class inherits from the AdfConfig class.
Currently this class does four things:

1.  Initializes an instance of AdfConfig.

2.  Checks for the three, required upper-level
    dictionaries specified in the config file,
    and makes copies where the variables have
    been expanded.

3.  Extract values for "compare_obs", "diag_var_list",
    and "plot_location", and provide properties to
    access these values to the rest of ADF.

4.  Set "num_procs" variable, and provde num_procs
    property to the rest of ADF.

This class also provide methods for extracting
variables from the standard, expanded config
dictionaries.
"""

#++++++++++++++++++++++++++++++
#Import standard python modules
#++++++++++++++++++++++++++++++

from pathlib import Path
import copy
import os
import getpass

#+++++++++++++++++++++++++++++++++++++++++++++++++
#import non-standard python modules, including ADF
#+++++++++++++++++++++++++++++++++++++++++++++++++

# pylint: disable=unused-import
import numpy as np
import xarray as xr
# pylint: enable=unused-import

#ADF modules:
from adf_config import AdfConfig
from adf_base   import AdfError

#+++++++++++++++++++
#Define Obs class
#+++++++++++++++++++

class AdfInfo(AdfConfig):

    """
    Information/Parameter class, which initializes
    an AdfConfig object and provides additional
    variables and methods to simplify access to the
    standard, expanded config dictionaries.
    """

    def __init__(self, config_file, debug=False):

        """
        Initalize ADF Info object.
        """

        #Initialize Config attributes:
        super().__init__(config_file, debug=debug)

        #Add basic diagnostic info to object:
        self.__basic_info = self.read_config_var('diag_basic_info', required=True)

        #Expand basic info variable strings:
        self.expand_references(self.__basic_info)

        #Add CAM climatology info to object:
        self.__cam_climo_info = self.read_config_var('diag_cam_climo', required=True)

        #Expand CAM climo info variable strings:
        self.expand_references(self.__cam_climo_info)

        # Add CVDP info to object:
        self.__cvdp_info = self.read_config_var("diag_cvdp_info")

        # Expand CVDP climo info variable strings:
        if self.__cvdp_info is not None:
            self.expand_references(self.__cvdp_info)
        # End if

        # Add MDTF info to object:
        self.__mdtf_info = self.read_config_var("diag_mdtf_info")

        if self.__mdtf_info is not None:
            if self.__mdtf_info['mdtf_run']:
                self.expand_references(self.__mdtf_info)
        # End if

        # Get the current system user
        self.__user = getpass.getuser()

        # Check if inputs are of the correct type:
        # -------------------------------------------

        #Use "cam_case_name" as the variable that sets the total number of cases:
        if isinstance(self.get_cam_info("cam_case_name", required=True), list):

            #Extract total number of test cases:
            self.__num_cases = len(self.get_cam_info("cam_case_name"))

        else:
            #Set number of cases to one:
            self.__num_cases = 1
        #End if

        #Loop over all items in config dict:
        for conf_var, conf_val in self.__cam_climo_info.items():
            # Hist_str can be a list for each case, so set it as a nested list here
            if "hist_str" in conf_var:
                self.hist_str_to_list(conf_var, conf_val)
            elif isinstance(conf_val, list):
                # If a list, then make sure it is has the correct number of entries:
                if not len(conf_val) == self.__num_cases:
                    emsg = f"diag_cam_climo config variable '{conf_var}' should have"
                    emsg += f" {self.__num_cases} entries, instead it has {len(conf_val)}"
                    self.end_diag_fail(emsg)
            else:
                #If not a list, then convert it to one:
                self.__cam_climo_info[conf_var] = [conf_val]
            # End if
        # End for
        # -------------------------------------------

        #Read hist_str (component.hist_num) from the yaml file, or set to default
        hist_str = self.__cam_climo_info['hist_str']
        #If hist_str is not present, then default to 'cam.h0':
        if not hist_str:
            hist_str = [['cam.h0']]*self.__num_cases
        #End if
        self.__hist_str = hist_str

        #Initialize ADF variable list:
        #self.__diag_var_list = self.read_config_var('diag_var_list', required=True)
        if self.read_config_var('diag_var_list'):
            self.__diag_var_list = self.read_config_var('diag_var_list')
        else:
            self.__diag_var_list = []


        #Case names:
        case_names = self.get_cam_info('cam_case_name', required=True)

        #Grab test case nickname(s)
        test_nickname_list = self.get_cam_info('case_nickname')

        if test_nickname_list:
            test_nicknames = [] #set to be an empty list
            for i,nickname in enumerate(test_nickname_list):
                if nickname is None:
                    test_nicknames.append(case_names[i])
                else:
                    test_nicknames.append(test_nickname_list[i])
                #End if
            #End for
        else:
            test_nicknames = [] #Re-set to be an empty list
            for case_name in case_names:
                test_nicknames.append(case_name)
            #End for
        #End if

        self.__base_hist_str = ""
        self.__baseline_ts_done = None
        self.__calc_bl_climo = True
        self.__bl_overwrite_ts = False
        self.__calc_baseline_ts = False
        

        #Initialize "compare_obs" variable:
        self.__compare_obs = self.get_basic_info('compare_obs')

        #Check if a CAM vs AMWG obs comparison is being performed:
        if self.__compare_obs:

            #If so, then set the baseline info to None, to ensure any scripts
            #that check this variable won't crash:
            self.__cam_bl_climo_info = None

            #Also set data name for use below:
            data_name = "Obs"
            base_nickname = "Obs"

            #Set the baseline years to empty strings:
            syear_baseline = ""
            eyear_baseline = ""

            input_ts_baseline = [None]
        #elif 'diag_cam_baseline_climo' not in self:
        #    print("no baseline case nor is this against obs")
        #    pass
        else:
            #If not, then assume a CAM vs CAM run and add CAM baseline climatology info to object:
            self.__cam_bl_climo_info = self.read_config_var('diag_cam_baseline_climo',
                                                            required=True)

            #Expand CAM baseline climo info variable strings:
            self.expand_references(self.__cam_bl_climo_info)

            #Set data name to baseline case name:
            data_name = self.get_baseline_info('cam_case_name', required=True)

            #Attempt to grab baseline start_years (not currently required):
            syear_baseline = self.get_baseline_info('start_year')
            eyear_baseline = self.get_baseline_info('end_year')

            #Get climo years for verification or assignment if missing
            baseline_hist_locs = self.get_baseline_info('cam_hist_loc')
            if baseline_hist_locs is None:
                baseline_hist_locs = [None]

            # Read hist_str (component.hist_num, eg cam.h0) from the yaml file
            baseline_hist_str = self.get_baseline_info("hist_str")
            if not isinstance(conf_val, list):
                baseline_hist_str = [baseline_hist_str]
            self.__baseline_hist_loc = baseline_hist_str

            #Check if any time series files are pre-made
            baseline_ts_done   = self.get_baseline_info("cam_ts_done")
            baseline_overwrite_ts   = self.get_baseline_info("cam_overwrite_ts")

            baseline_overwrite_climo   = self.get_baseline_info("cam_overwrite_climo")

            input_ts_baseline = self.get_baseline_info("cam_ts_loc")
            self.__bl_ts_locs = input_ts_baseline

            input_climo_baseline = self.get_baseline_info("cam_climo_loc")
            self.__bl_climo_locs = input_climo_baseline


            #calc_bl_ts = {}
            # Make new variable `calc_ts` in case the user does not want time series generation but 
            # need to use history files for diagnostics, ie MDTF, Tape Recorder, budget tables, etc.
            print(baseline_ts_done)
            print(input_ts_baseline)
            print(self.get_baseline_info("calc_cam_climo"))
            print(self.get_baseline_info("cam_hist_loc"))

            if (baseline_ts_done is None) and (input_ts_baseline is None) and (self.get_baseline_info("calc_cam_climo") is None) and (self.get_baseline_info("cam_hist_loc")):
                calc_bl_ts = False
            else:
                calc_bl_ts = True

            self.__calc_baseline_ts = {}
            self.__calc_baseline_ts[data_name] = calc_bl_ts


            print("baseline_ts_done",baseline_ts_done,"\n")
            if baseline_ts_done is None:
                baseline_ts_done = True
            self.__baseline_ts_done = baseline_ts_done
            #self.__baseline_ts_done = {data_name:baseline_ts_done}
            #input_ts_baseline = self.get_baseline_info("cam_ts_loc")


            #Check if any time series files are pre-made
            #baseline_overwrite_ts   = self.get_baseline_info("cam_overwrite_ts")
            print("baseline_overwrite_ts",baseline_overwrite_ts,"\n")
            if baseline_overwrite_ts is None:
                baseline_overwrite_ts = False
            self.__bl_overwrite_ts = baseline_overwrite_ts
            #self.__bl_overwrite_ts = {data_name:baseline_overwrite_ts}


            #Check if any time series files are pre-made
            #baseline_overwrite_ts   = self.get_baseline_info("cam_overwrite_ts")
            print("baseline_overwrite_climo",baseline_overwrite_climo,"\n")
            if baseline_overwrite_climo is None:
                baseline_overwrite_climo = False
            self.__bl_overwrite_climo = baseline_overwrite_climo



            


            if (baseline_ts_done) and (not input_ts_baseline) and (self.get_baseline_info("calc_cam_climo")):
                self.__calc_bl_climo = False
            #else:
            #    self.__calc_bl_climo = True

            """calc_bl_ts = {}
            # Make new variable `calc_ts` in case the user does not want time series generation but 
            # need to use history files for diagnostics, ie MDTF, Tape Recorder, budget tables, etc.
            print(baseline_ts_done)
            print(input_ts_baseline)
            print(self.get_baseline_info("calc_cam_climo"))
            print(self.get_baseline_info("cam_hist_loc"))

            if (baseline_ts_done is None) and (input_ts_baseline is None) and (self.get_baseline_info("calc_cam_climo") is None) and (self.get_baseline_info("cam_hist_loc")):
                calc_bl_ts[data_name] = False
            else:
                calc_bl_ts[data_name] = True"""
            """
            {"test": copy.copy(self.__calc_test_ts),
                "baseline": copy.copy(self.__calc_baseline_ts)}
            """

            #Check if time series files already exist,
            #if so don't rely on climo years from history location
            if (baseline_ts_done) and (input_ts_baseline):
                baseline_hist_locs = None

                #Grab baseline time series file location
                input_ts_loc = Path(input_ts_baseline)

                #Get years from pre-made timeseries file(s)
                found_syear_baseline, found_eyear_baseline = self.get_climo_yrs_from_ts(input_ts_loc, data_name)
                found_yr_range = np.arange(found_syear_baseline,found_eyear_baseline,1)

                #History file path isn't needed if user is running ADF directly on time series.
                #So make sure start and end year are specified:
                if syear_baseline is None:
                    msg = f"No given start year for {data_name}, "
                    msg += f"using first found year: {found_syear_baseline}"
                    print(msg)
                    syear_baseline = found_syear_baseline
                if syear_baseline not in found_yr_range:
                    msg = f"Given start year '{syear_baseline}' is not in current dataset "
                    msg += f"{data_name}, using first found year: {found_syear_baseline}\n"
                    print(msg)
                    syear_baseline = found_syear_baseline

                if eyear_baseline is None:
                    msg = f"No given end year for {data_name}, "
                    msg += f"using last found year: {found_eyear_baseline}"
                    print(msg)
                    eyear_baseline = found_eyear_baseline
                if eyear_baseline not in found_yr_range:
                    msg = f"Given end year '{eyear_baseline}' is not in current dataset "
                    msg += f"{data_name}, using first found year: {found_eyear_baseline}\n"
                    print(msg)
                    eyear_baseline = found_eyear_baseline
            # End if

            # Check if history file path exists:
            if any(baseline_hist_locs):
                #Check if user provided
                if not baseline_hist_str:
                    baseline_hist_str = ['cam.h0a']
                else:
                    #Make list if not already
                    if not isinstance(baseline_hist_str, list):
                        baseline_hist_str = [baseline_hist_str]
                #Initialize baseline history string list
                self.__base_hist_str = baseline_hist_str

                #Grab first possible hist string, just looking for years of run
                base_hist_str = baseline_hist_str[0]
                
                #hist_str = baseline_hist_str[0]
                starting_location = Path(baseline_hist_locs)
                file_list = sorted(starting_location.glob("*" + base_hist_str + ".*.nc"))
                # Partition string to find exactly where h-number is
                # This cuts the string before and after the `{hist_str}.` sub-string
                # so there will always be three parts:
                # before sub-string, sub-string, and after sub-string
                #Since the last part always includes the time range, grab that with last index (2)
                #NOTE: this is based off the current CAM file name structure in the form:
                #  $CASE.cam.h#.YYYY<other date info>.nc
                base_climo_yrs = [int(str(i).partition(f"{base_hist_str}.")[2][0:4]) for i in file_list]
                base_climo_yrs = sorted(np.unique(base_climo_yrs))

                #print("base_climo_yrs",base_climo_yrs,"\n")

                base_found_syr = int(base_climo_yrs[0])
                base_found_eyr = int(base_climo_yrs[-1])

                #Check if start or end year is missing. If so then just assume it is the
                #start or end of the entire available model data.
                if syear_baseline is None:
                    msg = f"No given start year for {data_name}, "
                    msg += f"using first found year: {base_found_syr}"
                    print(msg)
                    syear_baseline = base_found_syr
                if syear_baseline not in base_climo_yrs:
                    msg = f"Given start year '{syear_baseline}' is not in current dataset "
                    msg += f"{data_name}, using first found year: {base_climo_yrs[0]}\n"
                    print(msg)
                    syear_baseline = base_found_syr

                if eyear_baseline is None:
                    msg = f"No given end year for {data_name}, "
                    msg += f"using last found year: {base_found_eyr}"
                    print(msg)
                    eyear_baseline = base_found_eyr
                if eyear_baseline not in base_climo_yrs:
                    msg = f"Given end year '{eyear_baseline}' is not in current dataset "
                    msg += f"{data_name}, using last found year: {base_climo_yrs[-1]}\n"
                    print(msg)
                    eyear_baseline = base_found_eyr

                #Grab baseline nickname
                base_nickname = self.get_baseline_info('case_nickname')
                if base_nickname is None:
                    base_nickname = data_name
            #End if

            #Grab baseline nickname
            base_nickname = self.get_baseline_info('case_nickname')
            if base_nickname is None:
                base_nickname = data_name

            #Get integer for baseline years for searching climo files
            syear_baseline = int(syear_baseline)
            eyear_baseline = int(eyear_baseline)

            #Update baseline case name:
            data_name_long = data_name+f"_{syear_baseline}_{eyear_baseline}"
        #End if (compare_obs)

        #Initialize case nicknames:
        self.__test_nicknames = test_nicknames
        self.__base_nickname = base_nickname

        #Save starting and ending years as object variables:
        self.__syear_baseline = syear_baseline
        self.__eyear_baseline = eyear_baseline




        """
        #Check if any time series files are pre-made
        #baseline_overwrite_ts   = self.get_baseline_info("cam_overwrite_ts")
        print("baseline_overwrite_ts",baseline_overwrite_ts,"\n")
        if baseline_overwrite_ts is None:
            baseline_overwrite_ts = False
        self.__bl_overwrite_ts = {data_name:baseline_overwrite_ts}

        #Check if any time series files are pre-made
        #baseline_overwrite_ts   = self.get_baseline_info("cam_overwrite_ts")
        print("baseline_overwrite_ts",baseline_overwrite_ts,"\n")
        if baseline_overwrite_ts is None:
            baseline_overwrite_ts = False
        self.__bl_overwrite_ts = {data_name:baseline_overwrite_ts}
        """




        #Create plot location variable for potential use by the website generator.
        #Please note that this is also assumed to be the output location for the analyses scripts:
        #-------------------------------------------------------------------------
        self.__plot_location = [] #Must be a list to manage multiple cases

        #Plot directory:
        plot_dir = self.get_basic_info('cam_diag_plot_loc', required=True)
        print("plot_dir",plot_dir)

        #Case names:
        case_names = self.get_cam_info('cam_case_name', required=True)

        #Start years (not currently required):
        syears = self.get_cam_info('start_year')

        #End year (not currently rquired):
        eyears = self.get_cam_info('end_year')



        #Check if premade ts files - premade ts files
        ###########################################################
        #Start years (not currently required):
        syears = self.get_cam_info('start_year')
        if syears is None:
            syears = [None]*len(case_names)
        else:
            #Check if any time series files are pre-made
            if len(syears) == len(case_names):
                for i,syear in enumerate(syears):
                    if syear is None:
                        syears[i] = True
            else:
                print()

        self.__test_syears = {}
        for i,syear in enumerate(syears):
            self.__test_syears[case_names[i]] = syear

        test_syears = copy.copy(self.__test_syears)
        if self.__syear_baseline:
            baseline_syears = self.__syear_baseline
        else:
            baseline_syears = None

        print("test_syears",test_syears)
        #syears_dict = {"test":test_syears,"baseline":baseline_syears}
        syears_dict = {"test":test_syears,"baseline":{data_name:baseline_syears}}
        self.__syears_dict = syears_dict
        ###########################################################



        #Check if premade ts files - premade ts files
        ###########################################################
        #Start years (not currently required):
        eyears = self.get_cam_info('end_year')
        if eyears is None:
            eyears = [None]*len(case_names)
        else:
            #Check if any time series files are pre-made
            if len(eyears) == len(case_names):
                for i,eyear in enumerate(eyears):
                    if eyear is None:
                        eyears[i] = True
            else:
                print()

        self.__test_eyears = {}
        for i,eyear in enumerate(eyears):
            self.__test_eyears[case_names[i]] = eyear

        test_eyears = copy.copy(self.__test_eyears)
        if self.__eyear_baseline:
            baseline_eyears = self.__eyear_baseline
        else:
            baseline_eyears = None
        #eyears_dict = {"test":test_eyears,"baseline":baseline_eyears}
        eyears_dict = {"test":test_eyears,"baseline":{data_name:baseline_eyears}}
        self.__eyears_dict = eyears_dict
        ###########################################################






        """#Make lists of None to be iterated over for case_names
        if syears is None:
            syears = [None]*len(case_names)
        #End if
        if eyears is None:
            eyears = [None]*len(case_names)
        #End if"""



        
        """
        #Extract cam history files location:
        cam_hist_locs = self.get_cam_info('cam_hist_loc')
        if cam_hist_locs is None:
            cam_hist_locs = [None]*len(case_names)
        """

        

        #Check if premade ts files - premade ts files
        ###########################################################
        #Start years (not currently required):
        eyears = self.get_cam_info('end_year')
        cam_hist_locs = self.get_cam_info('cam_hist_loc')
        if cam_hist_locs is None:
            cam_hist_locs = [None]*len(case_names)
        else:
            #Check if any time series files are pre-made
            if len(cam_hist_locs) == len(case_names):
                for i,hist_loc in enumerate(cam_hist_locs):
                    if hist_loc is None:
                        cam_hist_locs[i] = True
            else:
                print()

        self.__test_hist_locs = {}
        for i,hist_loc in enumerate(cam_hist_locs):
            self.__test_hist_locs[case_names[i]] = hist_loc

        test_hist_locs = copy.copy(self.__test_hist_locs)
        if self.__baseline_hist_loc:
            baseline_hist_loc = self.__baseline_hist_loc
        else:
            baseline_hist_loc = None
        
        #hist_locs_dict = {"test":test_hist_locs,"baseline":baseline_hist_loc}
        hist_locs_dict = {"test":test_hist_locs,"baseline":{data_name:baseline_hist_loc}}
        self.__hist_locs_dict = hist_locs_dict
        ###########################################################













        # Read hist_str (component.hist_num, eg cam.h0) from the yaml file
        cam_hist_str = self.__hist_str

        #Check if premade ts files - premade ts files
        ###########################################################
        cam_ts_done   = self.get_cam_info("cam_ts_done")
        if cam_ts_done is None:
            cam_ts_done = [True]*len(case_names)
        else:
            #Check if any time series files are pre-made
            if len(cam_ts_done) == len(case_names):
                for i,case in enumerate(cam_ts_done):
                    if case is None:
                        cam_ts_done[i] = True
            else:
                print()

        self.__test_ts_done = {}
        for i,cam_ts in enumerate(cam_ts_done):
            self.__test_ts_done[case_names[i]] = cam_ts

        test_ts_done = copy.copy(self.__test_ts_done)
        if self.__baseline_ts_done:
            bl_ts_done = self.__baseline_ts_done
        else:
            bl_ts_done = True
        #ts_done_dict = {"test":test_ts_done,"baseline":bl_ts_done}
        ts_done_dict = {"test":test_ts_done,"baseline":{data_name:bl_ts_done}}
        
        self.__ts_done_dict = ts_done_dict
        ###########################################################







        #Check if using pre-made ts files, overwrite them - overwrite ts
        ###########################################################
        cam_overwrite_ts   = self.get_cam_info("cam_overwrite_ts")
        if cam_overwrite_ts is None:
            #cam_overwrite_ts = [False]*len(case_names)
            cam_overwrite_ts = [None]*len(case_names)
        else:
            #Check if any time series files are pre-made
            if len(cam_overwrite_ts) == len(case_names):
                for i,overwrite_ts in enumerate(cam_overwrite_ts):
                    if overwrite_ts is None:
                        #cam_overwrite_ts[i] = False
                        cam_overwrite_ts[i] = None
            else:
                print()
        self.__test_overwrite_ts = {}
        for i,cam_ts in enumerate(cam_overwrite_ts):
            self.__test_overwrite_ts[case_names[i]] = cam_ts

        test_overwrite_ts = copy.copy(self.__test_overwrite_ts)
        if self.__bl_overwrite_ts:
            bl_overwrite_ts = self.__bl_overwrite_ts
        else:
            #bl_overwrite_ts = False
            bl_overwrite_ts = None
        #overwrite_ts_dict = {"test":test_overwrite_ts,"baseline":bl_overwrite_ts}
        overwrite_ts_dict = {"test":test_overwrite_ts,"baseline":{data_name:bl_overwrite_ts}}

        self.__cam_overwrite_ts_dict = overwrite_ts_dict
        ###########################################################












        #Check if using pre-made ts files, overwrite them - overwrite ts
        ###########################################################
        cam_overwrite_climo   = self.get_cam_info("cam_overwrite_climo")
        if cam_overwrite_climo is None:
            #cam_overwrite_ts = [False]*len(case_names)
            cam_overwrite_climo = [None]*len(case_names)
        else:
            #Check if any time series files are pre-made
            if len(cam_overwrite_climo) == len(case_names):
                for i,overwrite_c in enumerate(cam_overwrite_climo):
                    if overwrite_c is None:
                        #cam_overwrite_ts[i] = False
                        cam_overwrite_climo[i] = None
            else:
                print()
        self.__overwrite_climo = {}
        for i,cam_ts in enumerate(cam_overwrite_climo):
            self.__overwrite_climo[case_names[i]] = cam_ts

        overwrite_climo = copy.copy(self.__overwrite_climo)
        if self.__bl_overwrite_climo:
            bl_overwrite_climo = self.__bl_overwrite_climo
        else:
            #bl_overwrite_ts = False
            bl_overwrite_climo = None
        #overwrite_ts_dict = {"test":test_overwrite_ts,"baseline":bl_overwrite_ts}
        overwrite_climo_dict = {"test":overwrite_climo,"baseline":{data_name:bl_overwrite_climo}}

        self.__overwrite_climo_dict = overwrite_climo_dict
        ###########################################################








        #Grab case time series file location(s) - input ts locs
        ###########################################################
        input_ts_locs = self.get_cam_info("cam_ts_loc")
        if input_ts_locs is None:
            input_ts_locs = [None]*len(case_names)
        else:
            #Check if any time series files are pre-made
            if len(input_ts_locs) == len(case_names):
                for i,case in enumerate(input_ts_locs):
                    if case is None:
                        input_ts_locs[i] = None
            else:
                print()

        self.__test_ts_locs = {}
        for i,ts_loc in enumerate(input_ts_locs):
            self.__test_ts_locs[case_names[i]] = ts_loc

        test_ts_locs = copy.copy(self.__test_ts_locs)
        if self.__bl_ts_locs:
            bl_ts_locs = self.__bl_ts_locs
        else:
            #bl_overwrite_ts = False
            bl_ts_locs = None
        #ts_locs_dict = {"test":test_ts_locs,"baseline":bl_ts_locs}
        ts_locs_dict = {"test":test_ts_locs,"baseline":{data_name:bl_ts_locs}}

        self.__cam_ts_locs_dict = ts_locs_dict
        ###########################################################












        #Grab case climo file location(s) - input ts locs
        ###########################################################
        test_climo_locs = self.get_cam_info("cam_climo_loc")
        if test_climo_locs is None:
            test_climo_locs = [None]*len(case_names)
        else:
            #Check if any time series files are pre-made
            if len(test_climo_locs) == len(case_names):
                for i,case in enumerate(test_climo_locs):
                    if case is None:
                        test_climo_locs[i] = None
            else:
                print()

        self.__test_climo_locs = {}
        for i,climo_loc in enumerate(test_climo_locs):
            self.__test_climo_locs[case_names[i]] = climo_loc

        test_climo_locs = copy.copy(self.__test_climo_locs)
        if self.__bl_climo_locs:
            bl_climo_locs = self.__bl_climo_locs
        else:
            #bl_overwrite_ts = False
            bl_climo_locs = None
        #climo_locs_dict = {"test":test_climo_locs,"baseline":bl_climo_locs}
        climo_locs_dict = {"test":test_climo_locs,"baseline":{data_name:bl_climo_locs}}

        self.__cam_climo_locs_dict = climo_locs_dict
        ###########################################################















        # Check if climos need to be calculated - calc test climo
        ###########################################################
        calc_test_climo = self.get_cam_info("calc_cam_climo")
        if calc_test_climo is None:
            #calc_test_climo = [False]*len(case_names)
            calc_test_climo = [None]*len(case_names)
        else:
            #Check if any time series files are pre-made
            if len(input_ts_locs) == len(case_names):
                for i,case in enumerate(calc_test_climo):
                    if case is None:
                        #calc_test_climo[i] = True
                        calc_test_climo[i] = None
            else:
                print()
        ###########################################################

        #Add check for obs!!!

        # Check if climo files need to be calculated - calc climo
        ###########################################################
        self.__calc_test_climo = {}
        for i in range(len(calc_test_climo)):
            #if (input_ts_locs[i]) and (not input_ts_baseline[i]) and (not calc_test_climo[i]):
            if (input_ts_locs[i]) and (not calc_test_climo[i]):
                self.__calc_test_climo[case_names[i]] = False
                #self.__calc_climo[i] = False
            else:
                self.__calc_test_climo[case_names[i]] = True

        calc_test_climo = copy.copy(self.__calc_test_climo)
        #calc_bl_climo = self.__calc_bl_climo
        if self.__calc_bl_climo:
            calc_bl_climo = self.__calc_bl_climo
        else:
            calc_bl_climo = True
        #calc_climo_dict = {"test":calc_test_climo,"baseline":calc_bl_climo}
        calc_climo_dict = {"test":calc_test_climo,"baseline":{data_name:calc_bl_climo}}

        self.__calc_climo_dict = calc_climo_dict
        ###########################################################

        #Loop over cases:
        #syears_fixed = []
        #eyears_fixed = []
        syears_fixed = {}
        eyears_fixed = {}
        ts_done = {}

        calc_test_ts = {}
        for case_idx, case_name in enumerate(case_names):

            syear = syears[case_idx]
            eyear = eyears[case_idx]


            print("OKAY, here we go...")
            try:
                test_ts_done = self.get_cam_info("cam_ts_done")[case_idx]
            except:
                test_ts_done = None

            #test_ts_done = self.get_cam_info("cam_ts_done")[case_idx]
            #test_ts_loc = self.get_cam_info("cam_ts_loc")[case_idx]
            #calc_test_climo = self.get_cam_info("calc_cam_climo")[case_idx]
            #cam_hist_loc = self.get_cam_info('cam_hist_loc')[case_idx]
            try:
                test_ts_loc = self.get_cam_info("cam_ts_loc")[case_idx]
            except:
                test_ts_loc = None
            try:
                calc_test_climo = self.get_cam_info("calc_cam_climo")[case_idx]
            except:
                calc_test_climo = None
            try:
                cam_hist_loc = self.get_cam_info('cam_hist_loc')[case_idx]
            except:
                cam_hist_loc = None
            

            print(test_ts_done)
            print(test_ts_loc)
            print(calc_test_climo)
            print(cam_hist_loc)

            #if (test_ts_done is None) and (test_ts_loc is None) and (calc_test_climo is None) and (cam_hist_loc):
            if (not test_ts_done) and (not test_ts_loc) and (not calc_test_climo) and (cam_hist_loc):
                calc_test_ts[case_name] = False
            else:
                calc_test_ts[case_name] = True



            

            #Check if time series files exist, if so don't rely on climo years
            if (cam_ts_done[case_idx]) and (input_ts_locs[case_idx]):
            #if ts_done[case_name]:
                cam_hist_locs[case_idx] = None

                #Grab case time series file location
                input_ts_loc = Path(input_ts_locs[case_idx])
                print(f"Checking existing time-series files in {input_ts_loc}")

                #Get years from pre-made timeseries file(s)
                found_syear, found_eyear = self.get_climo_yrs_from_ts(input_ts_loc, case_name)
                found_yr_range = np.arange(found_syear,found_eyear,1)

                #History file path isn't needed if user is running ADF directly on time series.
                #So make sure start and end year are specified:
                if syear is None:
                    msg = f"No given start year for {case_name}, "
                    msg += f"using first found year: {found_syear}"
                    print(msg)
                    syear = found_syear
                if syear not in found_yr_range:
                    msg = f"Given start year '{syear}' is not in current dataset "
                    msg += f"{case_name}, using first found year: {found_syear}\n"
                    print(msg)
                    syear = found_syear
                #End if
                if eyear is None:
                    msg = f"No given end year for {case_name}, "
                    msg += f"using last found year: {found_eyear}"
                    print(msg)
                    eyear = found_eyear
                if eyear not in found_yr_range:
                    msg = f"Given end year '{eyear}' is not in current dataset "
                    msg += f"{case_name}, using last found year: {found_eyear}\n"
                    print(msg)
                    eyear = found_eyear
                #End if
            #End if

            #Check if history file path exists:
            hist_str_case = cam_hist_str[case_idx]
            if any(cam_hist_locs):
                hist_str = hist_str_case[0]

                #Get climo years for verification or assignment if missing
                starting_location = Path(cam_hist_locs[case_idx])
                print("starting_location",starting_location,"\n")
                #ug = Path()
                print("WoOOooHAHahahAhsd",starting_location.is_dir())
                file_list = sorted(starting_location.glob('*'+hist_str+'.*.nc'))
                print("file_list",file_list)
                if len(file_list) == 0:
                    print("\tYeah, it's an empty list. Why did this not get checked before getting here. I mean come on.\n")
                else:
                    print()

                #file_list = sorted(starting_location.glob('*'+hist_str+'.*.nc'))
                #Partition string to find exactly where h-number is
                #This cuts the string before and after the `{hist_str}.` sub-string
                # so there will always be three parts:
                # before sub-string, sub-string, and after sub-string
                #Since the last part always includes the time range, grab that with last index (2)
                #NOTE: this is based off the current CAM file name structure in the form:
                #  $CASE.cam.h#.YYYY<other date info>.nc
                case_climo_yrs = [int(str(i).partition(f"{hist_str}.")[2][0:4]) for i in file_list]
                case_climo_yrs = sorted(np.unique(case_climo_yrs))
                print("case_climo_yrs",case_climo_yrs,type(case_climo_yrs))
                print(len(case_climo_yrs))
                if len(case_climo_yrs) == 0:
                    print("Yeah, it's an empty list. Why did this not get checked before getting here. I mean come on.\n")
                else:
                    print()

                case_found_syr = int(case_climo_yrs[0])
                case_found_eyr = int(case_climo_yrs[-1])

                #Check if start or end year is missing.  If so then just assume it is the
                #start or end of the entire available model data.
                if syear is None:
                    msg = f"No given start year for {case_name}, "
                    msg += f"using first found year: {case_found_syr}"
                    print(msg)
                    syear = case_found_syr
                if syear not in case_climo_yrs:
                    msg = f"Given start year '{syear}' is not in current dataset "
                    msg += f"{case_name}, using first found year: {case_climo_yrs[0]}\n"
                    print(msg)
                    syear = case_found_syr
                #End if
                if eyear is None:
                    msg = f"No given end year for {case_name}, "
                    msg += f"using last found year: {case_found_eyr}"
                    print(msg)
                    eyear = case_found_eyr
                if eyear not in case_climo_yrs:
                    msg = f"Given end year '{eyear}' is not in current dataset "
                    msg += f"{case_name}, using last found year: {case_climo_yrs[-1]}\n"
                    print(msg)
                    eyear = case_found_eyr
                #End if
            #End if

            #Update climo year lists in case anything changed
            syear = int(syear)
            eyear = int(eyear)
            #syears_fixed.append(syear)
            syears_fixed[case_name] = syear
            #eyears_fixed.append(eyear)
            eyears_fixed[case_name] = eyear

            #Update case name with provided/found years:
            case_name_long = case_name+f"_{syear}_{eyear}"

            #Set the final directory name and save it to plot_location:
            direc_name = f"{case_name_long}_vs_{data_name_long}"
            plot_loc = os.path.join(plot_dir, direc_name)
            self.__plot_location.append(plot_loc)

            #If first iteration, then save directory name for use by baseline:
            first_case_dir = ''
            if case_idx == 0:
                first_case_dir = direc_name
            #End if

            #Go ahead and make the diag plot location if it doesn't exist already
            diag_location = Path(plot_loc)
            if not diag_location.is_dir():
                print(f"\t    {diag_location} not found, making new directory")
                diag_location.mkdir(parents=True)
        #End for

        self.__syears = syears_fixed
        self.__eyears = eyears_fixed

        self.__calc_test_ts = calc_test_ts

        #Finally add baseline case (if applicable) for use by the website table
        #generator.  These files will be stored in the same location as the first
        #listed case.
        if not self.compare_obs:
            self.__plot_location.append(os.path.join(plot_dir, first_case_dir))
        #End if

        #-------------------------------------------------------------------------

        #Initialize "num_procs" variable:
        #-----------------------------------------
        temp_num_procs = self.get_basic_info('num_procs')

        if not temp_num_procs:
            #Variable not present, so set to a single processor:
            self.__num_procs = 1
        else:
            #Check if variable is a string and matches keyword:
            if isinstance(temp_num_procs, str) and \
               temp_num_procs.strip() == "*":

                #Set number of processors to total number of CPUs
                #on the node.  Please note that at some point this
                #may need to be replaced with a DASK implementation
                #instead:

                #First try to get CPUs allowed by OS for process to use:
                try:
                    self.__num_procs = len(os.sched_getaffinity(0))
                except AttributeError:
                    #Operating system doesn't support getaffinity, so try
                    #straight CPU number:
                    if os.cpu_count():
                        self.__num_procs = os.cpu_count()
                    else:
                        #Something is weird with this Operating System,
                        #so warn user and then try to run in serial mode:
                        wmsg = "WARNING!!!! ADF unable to determine how"
                        wmsg += " many processors are availble on this system,"
                        wmsg += " so defaulting to a single process/core."
                        print(wmsg)
                        self.__num_procs = 1
                    #End if
                #End except

            else:
                #If anything else, then try to convert to integer:
                try:
                    self.__num_procs = int(temp_num_procs)
                except ValueError:
                    #This variable has been set to something that
                    #can't be converted into an integer, so warn
                    #user and then try to run in serial mode:
                    wmsg = "WARNING!!!!  The 'num_procs' variable"
                    wmsg += f" has been set to '{temp_num_procs}'"
                    wmsg += " which cannot be converted to an integer."
                    wmsg += "\nThe ADF will now default to a single core"
                    wmsg += " and attempt to run."
                    print(wmsg)
                    self.__num_procs = 1
                #End except
            #End if
        #End if
        #Print number of processors being used to debug log (if requested):
        self.debug_log(f"ADF is running with {self.__num_procs} processors.")
        # -----------------------------------------

    #########
    def hist_str_to_list(self, conf_var, conf_val):
        """
        Make hist_str a nested list [ncases,nfiles] of the given value(s)
        """
        if isinstance(conf_val, list):
            hist_str = conf_val
        else:  # one case, one hist str
            hist_str = [
                conf_val
            ]
        self.__cam_climo_info[conf_var] = [hist_str]
        # -----------------------------------------

    #########

    # Create property needed to return "user" name to user:
    @property
    def user(self):
        """Return the "user" name if requested."""
        return self.__user

    # Create property needed to return "compare_obs" logical to user:
    @property
    def compare_obs(self):
        """Return the "compare_obs" logical to the user if requested."""
        return self.__compare_obs

    # Create property needed to return the number of test cases (num_cases) to user:
    @property
    def num_cases(self):
        """Return the "num_cases" integer value to the user if requested."""
        return self.__num_cases

    # Create property needed to return "diag_var_list" list to user:
    @property
    def diag_var_list(self):
        """Return a copy of the "diag_var_list" list to the user if requested."""
        #Note that a copy is needed in order to avoid having a script mistakenly
        #modify this variable, as it is mutable and thus passed by reference:
        return copy.copy(self.__diag_var_list)

    # Create property needed to return "basic_info" expanded dictionary to user:
    @property
    def basic_info_dict(self):
        """Return a copy of the "basic_info" list to the user if requested."""
        #Note that a copy is needed in order to avoid having a script mistakenly
        #modify this variable, as it is mutable and thus passed by reference:
        return copy.copy(self.__basic_info)

    # Create property needed to return "basic_info" expanded dictionary to user:
    @property
    def cam_climo_dict(self):
        """Return a copy of the "cam_climo_dict" list to the user if requested."""
        #Note that a copy is needed in order to avoid having a script mistakenly
        #modify this variable, as it is mutable and thus passed by reference:
        return copy.copy(self.__cam_climo_info)

    # Create property needed to return "basic_info" expanded dictionary to user:
    @property
    def baseline_climo_dict(self):
        """Return a copy of the "cam_bl_climo_info" list to the user if requested."""
        #Note that a copy is needed in order to avoid having a script mistakenly
        #modify this variable, as it is mutable and thus passed by reference:
        return copy.copy(self.__cam_bl_climo_info)

    # Create property needed to return "num_procs" to user:
    @property
    def num_procs(self):
        """Return the "num_procs" logical to the user if requested."""
        return self.__num_procs

    # Create property needed to return "plot_location" variable to user:
    @property
    def plot_location(self):
        """Return a copy of the '__plot_location' string list to user if requested."""
        #Note that a copy is needed in order to avoid having a script mistakenly
        #modify this variable:
        return copy.copy(self.__plot_location)

    # Create property needed to return the climo start (syear) and end (eyear) years to user:
    @property
    def climo_yrs(self):
        """Return the "syear" and "eyear" integer values to the user if requested."""
        syears = copy.copy(self.__syears) #Send copies so a script doesn't modify the original
        eyears = copy.copy(self.__eyears)
        return {"syears":syears,"eyears":eyears,
                "syear_baseline":self.__syear_baseline, "eyear_baseline":self.__eyear_baseline}


    # Create property needed to return the case nicknames to user:
    @property
    def case_nicknames(self):
        """Return the test case and baseline nicknames to the user if requested."""

        #Note that copies are needed in order to avoid having a script mistakenly
        #modify these variables, as they are mutable and thus passed by reference:
        test_nicknames = copy.copy(self.__test_nicknames)
        base_nickname = self.__base_nickname

        return {"test_nicknames":test_nicknames,"base_nickname":base_nickname}

    @property
    def hist_string(self):
        """ Return the CAM history string list to the user if requested."""
        cam_hist_strs = copy.copy(self.__hist_str)
        if self.__base_hist_str:
            base_hist_strs = copy.copy(self.__base_hist_str)
        else:
            base_hist_strs = ""
        hist_strs = {"test_hist_str":cam_hist_strs, "base_hist_str":base_hist_strs}
        return hist_strs

    @property
    def calc_climos(self):
        """ Return the history string name to the user if requested."""

        #Make list of all entries, similarly how the ADF does in various scripts
        calc_climos = []
        for key,val in self.__calc_climo_dict.items():
            if key == "test":
                for _,val2 in val.items():
                    calc_climos.append(val2)
            else: # baseline
                calc_climos.append(val)
        #The length of this list should always be the number of cases!
        #calc_climos = calc_climos + [self.__calc_bl_climo]

        return calc_climos

    @property
    def calc_climo_dict(self):
        """ Return the history string name to the user if requested."""
        return self.__calc_climo_dict

    @property
    def overwrite_climo_dict(self):
        """ Return the history string name to the user if requested."""
        return self.__overwrite_climo_dict




    @property
    def ts_done(self):
        """ Return the history string name to the user if requested."""

        #Make list of all entries, similarly how the ADF does in various scripts
        ts_done = []
        for key,val in self.__ts_done_dict.items():
            if key == "test":
                for _,val2 in val.items():
                    ts_done.append(val2)
            else: # baseline
                ts_done.append(val)
        return ts_done

    @property
    def ts_done_dict(self):
        """ Return the history string name to the user if requested."""
        return self.__ts_done_dict


    @property
    def cam_overwrite_ts(self):
        """ Return the history string name to the user if requested."""

        #Make list of all entries, similarly how the ADF does in various scripts
        overwrite_ts = []
        for key,val in self.__cam_overwrite_ts_dict.items():
            if key == "test":
                for _,val2 in val.items():
                    overwrite_ts.append(val2)
            else: # baseline
                overwrite_ts.append(val)
        return overwrite_ts

    @property
    def cam_overwrite_ts_dict(self):
        """ Return the history string name to the user if requested."""
        return self.__cam_overwrite_ts_dict


    @property
    def ts_locs_dict(self):
        """ Return the history string name to the user if requested."""
        return self.__cam_ts_locs_dict


    @property
    def climo_locs_dict(self):
        """ Return the history string name to the user if requested."""
        return self.__cam_climo_locs_dict


    # Create property needed to return whether to caculate time series files:
    @property
    def calc_ts(self):
        return {"test": copy.copy(self.__calc_test_ts),
                "baseline": copy.copy(self.__calc_baseline_ts)}


    # Create property needed to return whether to caculate time series files:
    @property
    def eyears(self):
        return {"test": copy.copy(self.__calc_test_ts),
                "baseline": copy.copy(self.__eyears_dict)}
    


    @property
    def eyears_dict(self):
        """ Return the history string name to the user if requested."""
        return self.__eyears_dict

    @property
    def syears_dict(self):
        """ Return the history string name to the user if requested."""
        return self.__syears_dict
    




    #########

    #Utility function to access expanded 'diag_basic_info' variables:
    def get_basic_info(self, var_str, required=False):
        """
        Return the config variable from 'diag_basic_info' as requested by
        the user.
        """

        return self.read_config_var(var_str,
                                    conf_dict=self.__basic_info,
                                    required=required)

    #########

    #Utility function to access expanded 'diag_cam_climo' variables:
    def get_cam_info(self, var_str, required=False):
        """
        Return the config variable from 'diag_cam_climo' as requested by
        the user.  """

        return self.read_config_var(var_str,
                                    conf_dict=self.__cam_climo_info,
                                    required=required)

    #########

    #Utility function to access expanded 'diag_cam_baseline_climo' variables:
    def get_baseline_info(self, var_str, required=False):
        """
        Return the config variable from 'diag_cam_baseline_climo' as requested by
        the user.  This function assumes that if the user is requesting it,
        then it must be required.
        """

        #Check if the cam baseline dictionary exists:
        if not self.__cam_bl_climo_info:
            #If required, then throw an error:
            if required:
                emsg = "get_baseline_info: Requested variable cannot be found"
                emsg += " because no baseline info exists.\n"
                emsg += "This is likely because an observational comparison is being done,"
                emsg += " so try adding 'required = False' to the get call."
                self.end_diag_fail(emsg)
            #End if

            #If not required, then return none:
            return None
        #End if

        #If basline dictionary exists, then search for variable like normal:
        return self.read_config_var(var_str,
                                    conf_dict=self.__cam_bl_climo_info,
                                    required=required)

    #########

    #Utility function to add a new model variable to the ADF (diag) variable list:
    def add_diag_var(self, var_str):
        """
        Adds a new variable to the ADF variable list
        """
        if var_str not in self.__diag_var_list:
            self.__diag_var_list.append(var_str)
        #End if

    #########

    # Utility function to access expanded 'diag_cvdp_info' variables
    def get_cvdp_info(self, var_str, required=False):
        """
        Return the config variable from 'diag_cvdp_info' as requested by
        the user. If 'diag_cvdp_info' is not found then try grabbing the
        variable from the top level of the YAML config file dictionary
        instead.
        """

        return self.read_config_var(
            var_str, conf_dict=self.__cvdp_info, required=required
        )

    #########

    # Utility function to access expanded 'diag_mdtf_info' variables
    def get_mdtf_info(self, var_str, required=False):
        """
        Return the config variable from 'diag_mdtf_info' as requested by
        the user. If 'diag_mdtf_info' is not found then try grabbing the
        variable from the top level of the YAML config file dictionary
        instead.
        """

        return self.read_config_var(
            var_str, conf_dict=self.__mdtf_info, required=required
        )


    #########

    # Utility function to grab climo years from pre-made time series files:
    def get_climo_yrs_from_ts(self, input_ts_loc, case_name):
        """
        Grab start and end climo years if none are specified in config file
        for pre-made time series file(s)

        Return
        ------
          - start year
          - end year
        """

        #Grab variable list
        var_list = self.diag_var_list

        #Create "Path" objects:
        input_location  = Path(input_ts_loc)

        #Check that time series input directory actually exists:
        if not input_location.is_dir():
            errmsg = f"Time series directory '{input_ts_loc}' not found.  Script is exiting."
            raise AdfError(errmsg)

        # Search for first variable in var_list to get a time series file to read
        # NOTE: it is assumed all the variables have the same dates!
        # Also, it is assumed that only h0 files should be climo-ed.
        for var in var_list:
            ts_files = sorted(input_location.glob(f"{case_name}*h0*.{var}.*nc"))
            if ts_files:
                break
            else:
                logmsg = "get years for time series:"
                logmsg = f"\tVar '{var}' not in dataset, skip to next to try and find climo years..."
                self.debug_log(logmsg)

        #Read in file(s)
        if len(ts_files) == 1:
            cam_ts_data = xr.open_dataset(ts_files[0], decode_times=True)
        else:
            cam_ts_data = xr.open_mfdataset(ts_files, decode_times=True, combine='by_coords')

        #Average time dimension over time bounds, if bounds exist:
        if 'time_bnds' in cam_ts_data:
            time_bounds_name = 'time_bnds'
        elif 'time_bounds' in cam_ts_data:
            time_bounds_name = 'time_bounds'
        else:
            time_bounds_name = None

        if time_bounds_name:
            time = cam_ts_data['time']
            #NOTE: force `load` here b/c if dask & time is cftime,
            #throws a NotImplementedError:

            time = xr.DataArray(cam_ts_data[time_bounds_name].load().mean(dim='nbnd').values,
                                dims=time.dims, attrs=time.attrs)
            cam_ts_data['time'] = time
            cam_ts_data.assign_coords(time=time)
            cam_ts_data = xr.decode_cf(cam_ts_data)

        #Extract first and last years from dataset:
        syr = int(cam_ts_data.time[0].dt.year.values)
        eyr = int(cam_ts_data.time[-1].dt.year.values)

        if eyr-syr >= 100:
            msg = f"WARNING: the found climo year range is large: {eyr-syr} years, "
            msg += "this may take a long time!"
            print(msg)

        return syr, eyr

#++++++++++++++++++++
#End Class definition
#++++++++++++++++++++