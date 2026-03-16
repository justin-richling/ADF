---
layout: default # Tells Jekyll to wrap this content with _layouts/default.html
title: Scam Tips # Shows up as the text in the browser tab
permalink: /scam_tips/ # Allows _layouts/default.html to find this page.
---

# Scam Tips

## To check the status of the run:
```tcsh
> cd $SCRATCH/cases/<your_case_name>
> more CaseStatus
-or-
> qstat -u <YOUR_USERNAME>
```
[qstat documentation](https://ncar-hpc-docs.readthedocs.io/en/latest/pbs/?h=qstat#qstat)

## To delete a run
```tcsh
qdel #######
```
where ###### is the Job ID from qstat. Before you can run again, you may need to clean up the run dir:
```tcsh
> cd $SCRATCH/cases/<your_case_name>
> ./case.build --clean
> ./case.build --clean-all
```

# Errors

## Job exceeded resource walltime
The error is saying the jobs did not specify a long enough queue slot and the it ran out of time before finishing. To determine how much time you need, look in the run/ subdir of your case, which can be found under `/glade/derecho/scratch/<username>/cases`. The queue length (JOB_WALLCLOCK_TIME) is set in the `env_workflow.xml` file in the case directory. Examine the times for the daily files that were written out before the job failed. That, plus your JOB_WALLCLOCK_TIME should give you an idea of how much more time you need. You can also look at the log `run.<casename>` for additional useful info.

You have two options:

### Adjust the wallclock limit
`./xmlchange JOB_WALLCLOCK_TIME=20:00:00 --subgroup=case.run` replacing 20 with however many hours you think you need. It may take a long time to get a job into the queue that requires a long runtime. The job may wait hours or days before it gets a slot to run depending on how busy the machine is. For shorter runtimes, this can be a good option.

### Break up the run into shorter segments
For longer runtimes, a better strategy is not to ask for a longer wallclock time but to break up the run into shorter segments that don't take too long to run. The jobs usually get in the queue quicker and finish faster. Use the xmlchange commands to tell the model to run a subset of days and repeat itself until all days are run.

For example, if you want to run 2 months (60 days), you could ask the model to run 10 days and resubmit itself 5 times, for a total of 60 days.
```tcsh
> cd /glade/derecho/scratch/<username>/cases/<your_case>
> ./xmlchange STOP_OPTION=ndays
> ./xmlchange STOP_N=10
> ./xmlchange RESUBMIT=5
> ./*.submit
```

Each time the model runs 10 days it will look at the RESUBMIT variable and if it is greater than 0 the model framework will subtract 1 from the total and resubmit the job to run another 10 days, restarting from where it left off. The env_run.xml file contains the STOP and RESUBMIT parameters controlling the length of the run. The xmlchange adjusts the value of these parameters in the env_run.xml file.
