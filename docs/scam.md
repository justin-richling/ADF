---
layout: default # Tells Jekyll to wrap this content with _layouts/default.html
title: Configure SCAM # Shows up as the text in the browser tab
permalink: /scam/ # Allows _layouts/default.html to find this page.
---

# Running SCAM

## Introduction
SCAM, the [Single Column Atmosphere Model](https://www.cesm.ucar.edu/models/simple/scam) is a single column model version of the [Community Atmospheric Model (CAM)](https://www.cesm.ucar.edu/models/cam), a global atmosphere model developed at NSF NCAR for the weather and climate research communities.

This page describes the basic setup for running SCAM. Once you have this configured you can [Create Nudged IOP forcing using CAM for use with SCAM]({{ site.baseurl }}/scam_nudge/) or [run scamTEST]({{ site.baseurl }}/scam_test/) to ensure your setup is correct.

## Login to Derecho
SCAM can be run from the Derecho supercomputer at NCAR. If you haven't already, [login to derecho and clone the repo]({{ site.baseurl }}/derecho/)

## Identify available projects/accounts
Running SCAM on the Derecho supercomputer requires access to CPU core-hours.  You can find out which projects/accounts are linked to your login via the Systems Accounting Manager webpage

&nbsp;&nbsp;&nbsp;&nbsp;[https://sam.ucar.edu](https://sam.ucar.edu)

or from the command line while logged into Derecho by providing a bogus project to the qinteractive command
```tcsh
> qinteractive -A P99999999 @derecho
```
The command will return a list of available accounts:
```tcsh
qsub: Invalid account for CPU usage, available accounts:
Project, Status, Active
Pxxxxxxxx, Normal, True
Pxxxxxxxx, Normal, True
,etc
```
where xxxxxxxx is the actual project where you have access. If there are no projects available you will need to reach out to whomever in your group/program provides core-hours and request an allocation. Note which projects are available to you.

To get more info on how many hours are available for each account, login to SAM as described above, click on a project under reports and view it.


## Create a scratch directory
You will need to configure a scratch directory where SCAM can write files. Each Derecho user should have a /glade/derecho/scratch directory. Create a scam dir under your scratch dir.
```tcsh
> ls /glade/derecho/scratch/<YOUR_USERNAME>
> mkdir /glade/derecho/scratch/<YOUR_USERNAME>/scam
```
