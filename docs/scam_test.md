---
layout: default # Tells Jekyll to wrap this content with _layouts/default.html
title: Run SCAMtest # Shows up as the text in the browser tab
permalink: /scam_test/ # Allows _layouts/default.html to find this page.
---

# Running SCAMtest
The scamtest scripts do a bit for bit (BFB) verification of the SCAM run. They are useful if you are running and getting unexpected results. qdif and tdif are difference from ops in the scamtest BFB comparison.

SCAM can be run from the Derecho supercomputer at NCAR. If you haven't already, [login to derecho and clone the repo]({{ site.baseurl }}/derecho/)

Login to Derecho, identify project/accounts you have access to, create a scratch directory. This and other useful SCAM information can be found on the [SCAM Configuration]({{ site.baseurl }}/scam/) page.

You will run two files (**after making the changes below**):
* First run create_scamtest.F2000.ne3_ne3_mg37.005.**new**.cold_off.derecho
* Then run create_scamtest.F2000.ne3_ne3_mg37.005.**scm.new**.cold_off.derecho

## Edit the .derecho files
Edit the two SCAM .derecho files. Add your project to the create_newcase line

> $CESMDIR/cime/scripts/create_newcase -compset $COMPSET  **--project Pxxxxxxxx** -res $RES -compiler $COMPILER -case $CASEDIR/$CASENAME  -pecount ${PES} --run-unsupported

Change the CASEDIR to point to your scratch scam dir.

> set CASEDIR=/glade/derecho/scratch/<YOUR_USERNAME>/scam

### Modify the scripts to point to the CESM collection directory
These collections are currently stored under John Truesdale's campaign dir so we are using that. Edit the SCAM .derecho scripts to set this directory.

<div style="padding: 16px; margin-bottom: 16pt; border-left: .25em solid #0969da; background-color: #f0f7ff; border-radius: 6px;">
  <p style="margin-top: 4px; margin-bottom: 0; color: #0969da;">
    <strong> Note:</strong> Need to move these to a more generic location
  </p>
</div>

> set CESMDIR=/glade/campaign/cgd/amp/jet/collections/$CAMDIRNAME

## Run the first file

<div style="padding: 16px; margin-bottom: 16pt; border-left: .25em solid #0969da; background-color: #f0f7ff; border-radius: 6px;">
  <p style="margin-top: 4px; margin-bottom: 0; color: #0969da;">
    <strong> Note:</strong> Need to describe what this file does here.
  </p>
</div>

```tcsh
> ./create_scamtest.F2000.ne3_ne3_mg37.005.new.cold_off.derecho
```

Wait for the run to complete. Instructions for checking on its progress are near the bottom of the [SCAM Tips]({{ site.baseurl }}/scam_tips/) page.

## Run the second file

<div style="padding: 16px; margin-bottom: 16pt; border-left: .25em solid #0969da; background-color: #f0f7ff; border-radius: 6px;">
  <p style="margin-top: 4px; margin-bottom: 0; color: #0969da;">
    <strong> Note:</strong> Need to describe what this file does here.
  </p>
</div>

Once the first run completes, copy the created netCDF file from the scratch scam run dir to just under scratch/<YOUR_USERNAME>. The file path is too long to access it from its original directory.
```tcsh
> cp /glade/derecho/scratch/<YOUR_USERNAME>/scam/scamtest.F2000.ne3_ne3_mg37.005.new.cold_off/run/scamtest.F2000.ne3_ne3_mg37.005.new.cold_off.cam.h1i.0001-01-01-00000.nc /glade/derecho/scratch/<YOUR_USERNAME>/.
```

Edit the iopfile dir in create_scamtest.F2000.ne3_ne3_mg37.005.scm.new.cold_off.derecho to point to the shortened location.

> iopfile ='/glade/derecho/scratch/<YOUR_USERNAME>/scamtest.F2000.ne3_ne3_mg37.005.new.cold_off.cam.h1i.0001-01-01-00000.nc'

### Run the file

```tcsh
> ./create_scamtest.F2000.ne3_ne3_mg37.005.scm.new.cold_off.derecho
```

