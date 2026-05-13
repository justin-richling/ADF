---
layout: default # Tells Jekyll to wrap this content with _layouts/default.html
title: Accessing Derecho # Shows up as the text in the browser tab
---

### Accessing Derecho
The Jupyter Notebooks and SCAM scripts available in this repo can be run from /glade on Derecho. If you are not already set up, follow the steps below to gain access and configure your environment.

<h1>Intro to the ADF<a class="headerlink" href="#intro-to-the-adf" title="Permalink to this headline">#</a></h1>
<p>This package is meant to be an update/upgrade to the much used and beloved AMWG diagnostics package used by the atmospheric/chemistry CESM community. It is still under active development and is near it’s first major version.</p>
<section id="what-is-it">
<h2>What is it?<a class="headerlink" href="#what-is-it" title="Permalink to this headline">#</a></h2>
<p>The ADF is an open source, community developed Python-based set of collection of analysis (averaging), re-gridding, and plotting scripts aimed at replacing the old AMWG Diagnostics package (NCL-based). It is designed to be a multipurpose tool and has been built to be somewhat flexible for customization.</p>
</section>
<section id="who-is-it-for">
<h2>Who is it for?<a class="headerlink" href="#who-is-it-for" title="Permalink to this headline">#</a></h2>
<p>The ADF is geared towards users who are running CAM or CAM-like (MPAS) simulations that are looking to compare their runs agains other CAM simulations, observations/reanalysis or model comparison sets (like CMIP/AMIP).</p>
</section>
<section id="key-features">
<h2>Key Features<a class="headerlink" href="#key-features" title="Permalink to this headline">#</a></h2>
<p> <strong>Flexible and Open Source</strong></p>
<p>   - ADF code is completely open to the public<br>
    - Users can modify their clone of the ADF to fit their needs if desired</p>
<p> <strong>Use of <a href="https://geocat-comp.readthedocs.io/en/v2023.10.1/index.html" target="_blank"><u>GeoCAT</u></a> functions (currently limited use)</strong></p>
<p> <strong>Use of <a href="https://www.redhat.com/en/topics/automation/what-is-yaml" target="_blank"><u>YAML configuration files</u></a></strong></p>
<p>   - helps to avoid changing source code</p>
<p> <strong>Option for use of multiple processors</strong></p>
<p> <strong>Installation via Git and Conda package manager</strong></p>
<p>   - CISL machines are set to run out of the box with required dependencies</p>
<p> <strong>Centralize vertical interpolation</strong></p>
<p> -&gt; Regridding and vertical interpolation script which interpolates all model variables with a vertical component onto a standard set of pressure levels</p>
<p> -&gt; Allows 3D model variable comparison against 3D observations</p>
<p>   - Assuming the observations are also on the same set of pressure levels</p>
<p> <strong>Enable interpolation on <a href="https://www.mmm.ucar.edu/models/mpas/" target="_blank"><u>MPAS</u></a> vertical coordinate</strong></p>
<p> -&gt; Checks for the MPAS height-based vertical coordinate, and if present it will enable pressure-to-pressure vertical interpolation required to get MPAS data onto the standard pressure levels used by the ADF</p>
</section>
<section id="types-of-adf-comparisons">
<h2>Types of ADF Comparisons<a class="headerlink" href="#types-of-adf-comparisons" title="Permalink to this headline">#</a></h2>
<p>There are essentially 3 types of comparisons:</p>
<ul class="simple">
<li><p>CAM vs CAM</p></li>
<li><p>CAM vs Observations/Reanalysis</p></li>
<li><p>CAM vs CMIP</p></li>
</ul>
<p> Each of which can be run as:</p>
<ul>
<li><p>Single Case Comparison - One test (experiment) case vs one baseline (control/target) case</p></li>
    <li><p>Multiple Case Comparison<a style="color:red;">**</a> - Multiple test cases vs one baseline case</p></li>
    <p style="margin-top: -5px;">&emsp;&emsp; <a style="color:red;">**</a> In progress, will not be part of this tutorial at the moment :(</p>
</ul>
<p>Each have their own requirements for the run-time config file. The requirements for each type of comparison will be addressed in their respective sections of this tutorial under the <strong>GUIDED EXAMPLES</strong> section.</p>
<div class="admonition attention">
<p class="admonition-title">Attention</p>
<p>Depending on how many years, plot types and variables you configure, the more time it takes</p>
<ul class="simple">
<li><p>Rough average for all default plots, 10-20 years for 30 variables is ~45 minutes</p></li>
</ul>
<div class="admonition note">
<p class="admonition-title">Note</p>
<p>This is an active area of development; there are potential ways of cutting time/processes</p>
</div>
</div>
</section>
<hr class="docutils" />
<section id="adf-basics">
<h2>ADF Basics<a class="headerlink" href="#adf-basics" title="Permalink to this headline">#</a></h2>
<p>Now let’s take a quick look at what the ADF actually does and the flow of it.</p>
<section id="adf-flow">
<h3>ADF Flow<a class="headerlink" href="#adf-flow" title="Permalink to this headline">#</a></h3>
<p>A simple look at the flow of the ADF<br>
 ↳ Create time series files from monthly history files<br>
  ↳ Create climatology files from either ADF generated or pre-existing time series files (ie CMIP)<br>
   ↳ Regrid Test case from Baseline case and vice-versa from climatology files<br>
    ↳ Run Analysis scripts<br>
     ↳ Run Plotting scripts<br>
      ↳ Generate Website pages</p>
<p>** Most parts of the ADF are optional and can be turned off for your desired need:</p>
<p>Potential Examples:</p>
<ul class="simple">
<li><p>If you only want regridded data, turn off plotting, analysis, website parts</p></li>
<li><p>If you want plotting and only need a few plot images, you can turn off the website part</p></li>
<li><p>If you only want time series files, you can turn off all other parts</p></li>
<li><p>If you don’t care about the statistics (AMWG) tables, turn off analysis part</p></li>
<li><p>etc.</p></li>
</ul>
</section>
<hr class="docutils" />
<section id="adf-layout">
<h3>ADF Layout<a class="headerlink" href="#adf-layout" title="Permalink to this headline">#</a></h3>
<p>A simple look at the structure of the ADF directories</p>
<div class="highlight-default notranslate"><div class="highlight"><pre><span></span>.
|-- config_amwg_default_plots.yaml
|-- config_cam_baseline_example.yaml
|-- env
|   `-- conda_environment.yaml
|-- jupyter_sample.ipynb
|-- lib
|   |-- adf_base.py
|   |-- adf_config.py
|   |-- adf_diag.py
|   |-- adf_info.py
|   |-- adf_obs.py
|   |-- adf_variable_defaults.yaml
|   |-- adf_web.py
|   |-- plotting_functions.py
|   |-- test
|   |   |-- pylintrc
|   |   `-- unit_tests
|   |       |-- pytest.ini
|   |       |-- test_adf_base.py
|   |       |-- test_adf_config.py
|   |       `-- test_files
|   |           |-- config_cam_double_nested.yaml
|   |           |-- config_cam_keywords.yaml
|   |           `-- config_cam_unset_var.yaml
|   `-- website_templates
|       |-- adf_diag.css
|       |-- NCAR.gif
|       |-- template.html
|       |-- template_index.html
|       |-- template_mean_diag.html
|       |-- template_mean_tables.html
|       |-- template_multi_case_index.html
|       |-- template_table.html
|       `-- template_var.html
|-- LICENSE
|-- README.md
|-- run_adf_diag
`-- scripts
    |-- analysis
    |   `-- amwg_table.py
    |-- averaging
    |   `-- create_climo_files.py
    |-- plotting
    |   |-- cam_taylor_diagram.py
    |   |-- global_latlon_map.py
    |   |-- global_latlon_vect_map.py
    |   |-- meridional_mean.py
    |   |-- polar_map.py
    |   |-- qbo.py
    |   |-- regional_map_multicase.py
    |   `-- zonal_mean.py
    `-- regridding
        `-- regrid_and_vert_interp.py
</pre></div>
</div>
</section>
<section id="adf-setup">
<!--<h3>ADF Setup<a class="headerlink" href="#adf-setup" title="Permalink to this headline">#</a></h3>-->
<h3>ADF Setup</h3>
<p>A simple look at the steps for using the ADF</p>
<p><a>0.</a> Run CAM/other simulation to get ADF input files <br>
   history files (ie: h0, h1, etc.)<br>
   time series files (ie: C/AMIP)</p>
<p><a>1.</a> Configure copies of yaml file(s) for the input files<br>
   config_cam_baseline_example.yaml<br>
   lib/adf_variable_defaults.yaml</p>
<p><a>2.</a> Run the ADF<br>
   simple command line call: <code class="docutils literal notranslate"><span class="pre">./run_adf_diag</span> <span class="pre">config_myadf.yaml</span></code>
   run ADF in Jupyter notebook</p>
<p><a>3.</a> View Diagnostics<br>
   Locally<br>
   Jupyter/JupyterHub<br>
   Publish to website</p>
