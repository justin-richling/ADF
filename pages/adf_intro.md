---
layout: default # Tells Jekyll to wrap this content with _layouts/default.html
title: Accessing Derecho # Shows up as the text in the browser tab
---

### Accessing Derecho
The Jupyter Notebooks and SCAM scripts available in this repo can be run from /glade on Derecho. If you are not already set up, follow the steps below to gain access and configure your environment.

#### Obtain access to /glade on Derecho
You will need to run the scripts from Derecho so you can access the ERA5 data stored in /glade. See the [Derecho documentation](https://ncar-hpc-docs.readthedocs.io/en/latest/compute-systems/derecho/) for instructions to obtain access.

#### Login to Derecho
See [Logging in with Duo](https://ncar-hpc-docs.readthedocs.io/en/latest/getting-started/accounts/duo/#logging-in-with-duo) for details on authenticating
```tcsh
> ssh -Y <YOUR_USERNAME>@derecho.hpc.ucar.edu
```
* When it asks for ncar-two-factor, enter your CIT password.
* Authenticate with DUO

#### Create a collections directory to work in
```tcsh
> mkdir -p $HOME/collections
```

#### Clone this repo to your /glade $HOME/collections directory
```tcsh
> cd $HOME/collections
> git clone http://github.com/NCAR/INFORM-COMPASS-cookbook
```
or if you already cloned the repo, pull the latest updates
```tcsh
> cd INFORM-COMPASS-cookbook/
> git pull
```




<h1>Intro to the ADF<a class="headerlink" href="#intro-to-the-adf" title="Permalink to this headline">#</a></h1>
<p>This package is meant to be an update/upgrade to the much used and beloved AMWG diagnostics package used by the atmospheric/chemistry CESM community. It is still under active development and is near it’s first major version.</p>
<section id="what-is-it">
<h2>What is it?<a class="headerlink" href="#what-is-it" title="Permalink to this headline">#</a></h2>
<p>The ADF is an open source, community developed Python-based set of collection of analysis (averaging), re-gridding, and plotting scripts aimed at replacing the old AMWG Diagnostics package (NCL-based). It is designed to be a multipurpose tool and has been built to be somewhat flexible for customization.</p>



