---
layout: default # Tells Jekyll to wrap this content with _layouts/default.html
title: Run SCAM nudging # Shows up as the text in the browser tab
permalink: /scam_nudge/ # Allows _layouts/default.html to find this page.
---

# Create Nudged IOP forcing using CAM for use with SCAM

This procedure will generate IOP forcing data associated with the dates and area of the SOCRATES field campaign to use with SCAM.

* The first experiment will provide initial conditions that approximates the state of the atmosphere at daily time intervals throughout the period of the SOCRATES field campaign.  The ERA5 reanalysis data will be used to nudge the thermodynamic vertical profiles via the T,U,V and Q fields for the initial CAM run.
* The second experiment will also be a full 3d cam run which uses the initial condition/restart boundary data along with the CAMIOP and windowing capability of the nudging functionality to generate CAM IOP forcing that can be used with SCAM to rerun the state of any individual column.
* The third experiment runs the single column version of CAM using the initial condition data along with the IOP forcing to rerun a specific column of the atmosphere during the SOCRATES period.

The output of the SCAM run can be compared to the same column of the initial CAM run to see how the model atmosphere evolves (away?) from the nudged (observed) atmospheric state.  Many SCAM runs can be made to analyze the physics processes and modify the parameterizations to help improve the prognosed state.

There are three scripts for making these runs under the SCAM_scripts directory:
* First: CAM6 run to create a ground truth dataset that is nudged to ERA5 reanalysis
```tcsh
SCAM_scripts/create_CAM6_ne30_Global_Nudged_SOCRATES_Jan-Feb_2018
```
* Second: CAM6 run(s) to generate CAM IOP data. CAM will nudge to ERA5 reanalysis outside the SOCRATES area
```tcsh
SCAM_scripts/create_CAM6_ne30_Window_Nudged_SOCRATES_CAMIOP_Jan-18-19_RF01
```
* Third: single column atmosphere (SCAM) run(s) using the generated CAM IOP data.
```tcsh
SCAM_scripts/create_CAM6_ne30_SCAM_RUN
```

The steps outlined below place the CAM code and COMPASS-cookbook underneath your $HOME directory.  The 3 CAM/SCAM cases that are created from the cookbook scripts are located under your scratch space on Derecho.  The CAM experiments will generated a terabyte of data which can be handled by $SCRATCH.  These initial cases are writing out a lot of data for analysis as we are fine tuning our procedures. The final requirements will be much less. Since the SCAM experiment is just a single column it always puts out much smaller data sets and can be easily run on any filesystem.

Before you begin, you may find it useful to review the [SCAM Tips]({{ site.baseurl }}/scam_tips/) page.

1. Create the collections and case directories and check out your own version of the CAM code.
```tcsh
> mkdir -p /glade/derecho/scratch/$USER/cases
> mkdir -p $HOME/collections
> cd $HOME/collections
> git clone https://github.com/jtruesdal/CAM-1 CAM_6_4_120_compass
> cd CAM_6_4_120_compass
> git checkout compass
> ./bin/git-fleximod update
```
1. If you haven't already, checkout the COMPASS cookbook repository and navigate to the SCAM_scripts dir.
```tcsh
> cd $HOME/collections
> git clone https://github.com/NCAR/INFORM-COMPASS-cookbook.git
> cd INFORM-COMPASS-cookbook/SCAM_scripts
```

1. If you changed the location of the code or case directories, edit the CESMDIR and CASEDIR variables in the following script to point to your new locations.
```tcsh
> vi create_CAM6_ne30_Global_Nudged_SOCRATES_Jan-Feb_2018
> vi create_CAM6_ne30_Window_Nudged_SOCRATES_CAMIOP_Jan-18-19_RF01
> vi create_CAM6_ne30_SCAM_RUN
```

1. The case title is set in the scripts using a combination of the model resolution, compset and other specifics about the experiment.  If you would like to rename the case for this experiment then edit the following line to set CASENAME as you wish.  The script will stop if you try to overwrite a previous case.
```tcsh
 set CASENAME=${CASETITLE}.${COMPSET}.${RES}.${CASEID}.${EXP}
```

1. If you have not already, you will need to set a PBS_ACCOUNT environment variable in your $HOME/.tcshrc file to provide your project number.
```tcsh
> vi ~/.tcshrc
  > setenv PBS_ACCOUNT "P########"
> source ~/.tcshrc (just needed the first time, will be run automatically each time you login in the future)
```

1. Run the first globally nudged experiment.
```tcsh
> cd $HOME/collections/INFORM-COMPASS-cookbook/SCAM_scripts
> qcmd -- ./create_CAM6_ne30_Global_Nudged_SOCRATES_Jan-Feb_2018
```
You can check the status of the run, or delete it, using the scam commands described on the [Configure Scam]({{ site.baseurl }}/scam/) page.

1. After the first experiment finishes, you should have output data underneath $SCRATCH/cases/your_case_name/run.  See what you have!
```tcsh
> cd /glade/derecho/scratch/$USER/cases/f.e30.cam6_4_120.FHIST_BGC.ne30_ne30_mg17.SOCRATES_nudgeUVTQsoc_full_withCOSP_tau6h_2months_inithist.100.cosp/run
> ls -al *.cam.h*
```

1. Set up the second experiment to generate the IOP data for the SCAM run.
*  Modify the following script variables to specify the dates that you want to generate IOP data for. As an example the following variables are set for the first SOCRATES flight RF01 that began Jan 15 2018.
```tcsh
set RUN_STARTDATE=2018-01-15
set STOP_OPTION=ndays
set STOP_N=3
set REST_OPTION=${STOP_OPTION}
set REST_N=${STOP_N}
set RUN_REFCASE=f.e30.cam6_4_120.FHIST_BGC.ne30_ne30_mg17.SOCRATES_nudgeUVTQsoc_full_withCOSP_tau6h_2months_inithist.100.cosp
set RUN_REFDATE=2018-01-15
set RUN_REFDIR=/glade/derecho/scratch/$USER/cases/${RUN_REFCASE}/run
set GET_REFCASE=TRUE
```

1. Run the second experiment to generate IOP data for SCAM.
```tcsh
> cd $HOME/collections/INFORM-COMPASS-cookbook/SCAM_scripts
> qcmd -- ./create_CAM6_ne30_Window_Nudged_SOCRATES_CAMIOP_Jan-18-19_RF01
```

1. See what the second experiment generated
```tcsh
> cd /glade/derecho/scratch/$USER/cases/f.e30.cam6_4_120.FHIST_BGC.ne30_ne30_mg17.SOCRATES_nudgeUVTQwindow_withCOSP_tau6h_3days_camiop.rf01.cosp/run
>
```

1. Set up for the third experiment, the SCAM run:  **NOTE**: We should change exp 2 to write out a number of days of data to the iop history file.  Here we only copy the first days worth of data to the IOP file and only run the SCAM case for 47 timesteps.  We have 3 individual days of IOP data, we could cat them all together and copy that large IOP file over to SCRATCH or just have exp 2 write a number of days worth of data in one history file.
* Copy the IOP file from exp 2 for the correct dates to $SCRATCH and modify the CAM namelist variable iopfile to point to the copied IOP file.
* modify create_CAM6_ne30_SCAM_RUN script to set REFCASE variables, paths, and dates as done for step 6.
* set PTS_LAT and PTS_LON variables in the script to point to the column you would like to simulate.  NOTE the PTS_LAT and PTS_LON should point to a column in SOCRATES area these are not.  I'll have to rerun with something along the RF01 flight path.
```tcsh
> cp /glade/derecho/scratch/$USER/cases/f.e30.cam6_4_120.FHIST_BGC.ne30_ne30_mg17.SOCRATES_nudgeUVTQwindow_withCOSP_tau6h_3days_camiop.rf01.cosp.cam.h0i.2018-01-15-00000.nc /glade/derecho/scratch/$USER/r01.IOP.nc'
> cd $HOME/collections/INFORM-COMPASS-cookbook/SCAM_scripts
> emacs create_CAM6_ne30_SCAM_RUN
```
* modify the following line to point to your iop file
```tcsh
iopfile = '/glade/derecho/scratch/jet/rf01.IOP.nc'
```

* modify PTS_LAT and PTS_LON to point to the column you want to simulate
```
set PTS_LON=276.7082039324993
set PTS_LAT=44.80320177421346
```

1. Run SCAM
```tcsh
> cd $HOME/collections/INFORM-COMPASS-cookbook/SCAM_scripts
> qcmd -- ./create_CAM6_ne30_SCAM_RUN
```

These experiments should be analyzed and improved through several iterations.  Some items for consideration:
* How long a spin up is needed to bring the CAM into a quasi equilibrated state for the SOCRATES start dates?  Our first run started 2 weeks before SOCRATES start.
* What model variables should be nudged and what nudging parameters work best to achieve a state that is close to the obs but not too far from CAM equilibrium
  -  Exp. 1 and 2 used 6 hour nudging on T,U,V, and Q. Should we just be nudging with T,U, and V? Is a 6hr Tau nudge timescale too strong?  Do we need to nudge hourly?
  -  Exp. 1 is a global nudge to bring the SOCRATES area close to the reanalysis state
  -  Exp. 2 is nudged outside the SOCRATES region using the windowing feature of CAM nudging to allow the physics parameterizations in the SOCRATES area to evolve freely
* How long is the observed state maintained inside the windowed area during the SOCRATES period?
  -  Do we need to use initial nudged IC data before each flight?  How quickly does the model initial state degrade for the various observed weather regimes?
* SCAM only works with boundary data on the dynamics grid.  The physics grid can't be used.  I interpolated Isla's ne30pg3 ERA forcing to the ne30np4 grid but it might be better to have Isla regenerate the forcing data for the ne30np4 grid directly.
* How close is the nudged CAM state to the ERA reanalysis?
