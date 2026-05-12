---
layout: frontpage
title: INFORM-COMPASS
banner-title: AMWG Diagnostics Framework (ADF)
banner-button-text: Download
banner-button-url: https://github.com/NCAR/ADF
---

<div markdown="0">
<a href="https://github.com/NCAR/INFORM-COMPASS-cookbook/blob/main/LICENSE" style="margin-right: 8px;"><img src="https://img.shields.io/github/license/NCAR/INFORM-COMPASS-cookbook?color=2043a1" alt="GitHub License"></a>
<a href="https://github.com/NCAR/INFORM-COMPASS-cookbook/blob/main/CHANGELOG.md" style="margin-right: 8px;"><img src="https://img.shields.io/badge/updates-CHANGELOG-4a69bd" alt="Changelog"></a>
<img src="https://img.shields.io/badge/Status-under%20development-yellow" alt="Status Badge">
<!-- Placeholder for a future DOI badge
        [![DOI](https://img.shields.io/badge/DOI-10.#####)](https://doi.org/10.#####)-->
</div>

### Introduction

<div style="padding: 16px; margin-bottom: 16pt; border-left: .25em solid #0969da; background-color: #f0f7ff; border-radius: 6px;">
  <p style="margin-top: 4px; margin-bottom: 0; color: #0969da;">
    <strong> Note:</strong> Need to assign someone to fill out this content.
  </p>
</div>

The goal of INFORM-COMPASS is to ....

This cookbook ....

#### CITATION

As this repository is under development, please reach out to us before citing this work and we will provide a DOI and guidance at that time.












<!DOCTYPE html>

<html>
  <head>
    <meta charset="utf-8" />
    <meta name="viewport" content="width=device-width, initial-scale=1.0" /><meta name="generator" content="Docutils 0.17.1: http://docutils.sourceforge.net/" />

    <title>Welcome to the AMWG Diagnostics Framework (ADF) Tutorial &#8212; ADF Tutorial</title>
    
  <!-- Loaded before other Sphinx assets -->
  <link href="_static/styles/theme.css?digest=1999514e3f237ded88cf" rel="stylesheet">
<link href="_static/styles/pydata-sphinx-theme.css?digest=1999514e3f237ded88cf" rel="stylesheet">

    
  <link rel="stylesheet"
    href="_static/vendor/fontawesome/5.13.0/css/all.min.css">
  <link rel="preload" as="font" type="font/woff2" crossorigin
    href="_static/vendor/fontawesome/5.13.0/webfonts/fa-solid-900.woff2">
  <link rel="preload" as="font" type="font/woff2" crossorigin
    href="_static/vendor/fontawesome/5.13.0/webfonts/fa-brands-400.woff2">

    <link rel="stylesheet" type="text/css" href="_static/pygments.css" />
    <link rel="stylesheet" href="_static/styles/sphinx-book-theme.css?digest=5115cc725059bd94278eecd172e13a965bf8f5a9" type="text/css" />
    <link rel="stylesheet" type="text/css" href="_static/togglebutton.css" />
    <link rel="stylesheet" type="text/css" href="_static/copybutton.css" />
    <link rel="stylesheet" type="text/css" href="_static/mystnb.css" />
    <link rel="stylesheet" type="text/css" href="_static/sphinx-thebe.css" />
    <link rel="stylesheet" type="text/css" href="_static/mine.css" />
    <link rel="stylesheet" type="text/css" href="_static/.ipynb_checkpoints/mine-checkpoint.css" />
    <link rel="stylesheet" type="text/css" href="_static/design-style.b7bb847fb20b106c3d81b95245e65545.min.css" />
    
  <!-- Pre-loaded scripts that we'll load fully later -->
  <link rel="preload" as="script" href="_static/scripts/pydata-sphinx-theme.js?digest=1999514e3f237ded88cf">

    <script data-url_root="./" id="documentation_options" src="_static/documentation_options.js"></script>
    <script src="_static/jquery.js"></script>
    <script src="_static/underscore.js"></script>
    <script src="_static/doctools.js"></script>
    <script src="_static/clipboard.min.js"></script>
    <script src="_static/copybutton.js"></script>
    <script src="_static/scripts/sphinx-book-theme.js?digest=9c920249402e914e316237a7dbc6769907cce411"></script>
    <script>let toggleHintShow = 'Click to show';</script>
    <script>let toggleHintHide = 'Click to hide';</script>
    <script>let toggleOpenOnPrint = 'true';</script>
    <script src="_static/togglebutton.js"></script>
    <script>var togglebuttonSelector = '.toggle, .admonition.dropdown, .tag_hide_input div.cell_input, .tag_hide-input div.cell_input, .tag_hide_output div.cell_output, .tag_hide-output div.cell_output, .tag_hide_cell.cell, .tag_hide-cell.cell';</script>
    <script src="_static/design-tabs.js"></script>
    <script>const THEBE_JS_URL = "https://unpkg.com/thebe@0.8.2/lib/index.js"
const thebe_selector = ".thebe,.cell"
const thebe_selector_input = "pre"
const thebe_selector_output = ".output, .cell_output"
</script>
    <script async="async" src="_static/sphinx-thebe.js"></script>
    <link rel="index" title="Index" href="genindex.html" />
    <link rel="search" title="Search" href="search.html" />
    <link rel="next" title="Intro to the ADF" href="notebooks/intro/adf.html" />
    <meta name="viewport" content="width=device-width, initial-scale=1" />
    <meta name="docsearch:language" content="None">
    

    <!-- Google Analytics -->
    
  </head>
  <body data-spy="scroll" data-target="#bd-toc-nav" data-offset="60">
<!-- Checkboxes to toggle the left sidebar -->
<input type="checkbox" class="sidebar-toggle" name="__navigation" id="__navigation" aria-label="Toggle navigation sidebar">
<label class="overlay overlay-navbar" for="__navigation">
    <div class="visually-hidden">Toggle navigation sidebar</div>
</label>
<!-- Checkboxes to toggle the in-page toc -->
<input type="checkbox" class="sidebar-toggle" name="__page-toc" id="__page-toc" aria-label="Toggle in-page Table of Contents">
<label class="overlay overlay-pagetoc" for="__page-toc">
    <div class="visually-hidden">Toggle in-page Table of Contents</div>
</label>
<!-- Headers at the top -->
<div class="announcement header-item noprint"></div>
<div class="header header-item noprint"></div>

    
    <div class="container-fluid" id="banner"></div>

    

    <div class="container-xl">
      <div class="row">
          
<!-- Sidebar -->
<div class="bd-sidebar noprint" id="site-navigation">
    <div class="bd-sidebar__content">
        <div class="bd-sidebar__top"><div class="navbar-brand-box">
    <a class="navbar-brand text-wrap" href="index.html">
      
        <!-- `logo` is deprecated in Sphinx 4.0, so remove this when we stop supporting 3 -->
        
      
      
      <img src="_static/_6e5c5f85-aa2c-401c-995a-c8a6e5da0811.jpeg" class="logo" alt="logo">
      
      
      <h1 class="site-logo" id="site-title">ADF Tutorial</h1>
      
    </a>
</div><form class="bd-search d-flex align-items-center" action="search.html" method="get">
  <i class="icon fas fa-search"></i>
  <input type="search" class="form-control" name="q" id="search-input" placeholder="Search this book..." aria-label="Search this book..." autocomplete="off" >
</form><nav class="bd-links" id="bd-docs-nav" aria-label="Main">
    <div class="bd-toc-item active">
        
        <ul class="nav bd-sidenav bd-sidenav__home-link">
            <li class="toctree-l1 current active">
                <a class="reference internal" href="#">
                    Welcome to the AMWG Diagnostics Framework (ADF) Tutorial
                </a>
            </li>
        </ul>
        <p aria-level="2" class="caption" role="heading">
 <span class="caption-text">
  ADF
 </span>
</p>
<ul class="nav bd-sidenav">
 <li class="toctree-l1">
  <a class="reference internal" href="notebooks/intro/adf.html">
   Intro to the ADF
  </a>
 </li>
 <li class="toctree-l1">
  <a class="reference internal" href="notebooks/intro/output.html">
   What do you get?
  </a>
 </li>
</ul>
<p aria-level="2" class="caption" role="heading">
 <span class="caption-text">
  Installing the ADF
 </span>
</p>
<ul class="nav bd-sidenav">
 <li class="toctree-l1">
  <a class="reference internal" href="notebooks/install/setup.html">
   Setting up ADF
  </a>
 </li>
</ul>
<p aria-level="2" class="caption" role="heading">
 <span class="caption-text">
  ADF Basics
 </span>
</p>
<ul class="nav bd-sidenav">
 <li class="toctree-l1">
  <a class="reference internal" href="notebooks/basics/adf_basics.html">
   Intro
  </a>
 </li>
 <li class="toctree-l1">
  <a class="reference internal" href="notebooks/basics/yaml_basics.html">
   YAML files
  </a>
 </li>
 <li class="toctree-l1 has-children">
  <a class="reference internal" href="notebooks/basics/run_basics.html">
   Setting up
   <strong>
    your
   </strong>
   ADF run - optional
  </a>
  <input class="toctree-checkbox" id="toctree-checkbox-1" name="toctree-checkbox-1" type="checkbox"/>
  <label for="toctree-checkbox-1">
   <i class="fas fa-chevron-down">
   </i>
  </label>
  <ul>
   <li class="toctree-l2">
    <a class="reference internal" href="notebooks/basics/run_basics/adf_run.html">
     Running
     <strong>
      your
     </strong>
     ADF
    </a>
   </li>
   <li class="toctree-l2">
    <a class="reference internal" href="notebooks/basics/run_basics/outputs.html">
     Exploring
     <strong>
      your
     </strong>
     ADF outputs
    </a>
   </li>
  </ul>
 </li>
</ul>
<p aria-level="2" class="caption" role="heading">
 <span class="caption-text">
  Guided Examples
 </span>
</p>
<ul class="nav bd-sidenav">
 <li class="toctree-l1">
  <a class="reference internal" href="notebooks/basic_examples/readme.html">
   ADF Demo Setup
  </a>
 </li>
 <li class="toctree-l1">
  <a class="reference internal" href="notebooks/basic_examples/model_model_run.html">
   CAM vs CAM
  </a>
 </li>
 <li class="toctree-l1">
  <a class="reference internal" href="notebooks/basic_examples/model_obs_run.html">
   CAM vs Obs
  </a>
 </li>
 <li class="toctree-l1 has-children">
  <a class="reference internal" href="notebooks/exercises/custom_defaults.html">
   Custom Variable Defaults For the ADF
  </a>
  <input class="toctree-checkbox" id="toctree-checkbox-2" name="toctree-checkbox-2" type="checkbox"/>
  <label for="toctree-checkbox-2">
   <i class="fas fa-chevron-down">
   </i>
  </label>
  <ul>
   <li class="toctree-l2">
    <a class="reference internal" href="notebooks/exercises/custom_defaults/single_variable.html">
     Customize Single Variable
    </a>
   </li>
   <li class="toctree-l2">
    <a class="reference internal" href="notebooks/exercises/custom_defaults/custom_file.html">
     Customized File
    </a>
   </li>
  </ul>
 </li>
 <li class="toctree-l1">
  <a class="reference internal" href="notebooks/exercises/custom_observations.html">
   Custom Set of Observations
  </a>
 </li>
 <li class="toctree-l1">
  <a class="reference internal" href="notebooks/exercises/custom_script.html">
   Adding custom script to your ADF
  </a>
 </li>
 <li class="toctree-l1">
  <a class="reference internal" href="notebooks/basic_examples/jupyter.html">
   ADF in Jupyter
  </a>
 </li>
</ul>
<p aria-level="2" class="caption" role="heading">
 <span class="caption-text">
  ADF - Reference
 </span>
</p>
<ul class="nav bd-sidenav">
 <li class="toctree-l1 has-children">
  <a class="reference internal" href="notebooks/reference/intro.html">
   ADF Basics
  </a>
  <input class="toctree-checkbox" id="toctree-checkbox-3" name="toctree-checkbox-3" type="checkbox"/>
  <label for="toctree-checkbox-3">
   <i class="fas fa-chevron-down">
   </i>
  </label>
  <ul>
   <li class="toctree-l2">
    <a class="reference internal" href="notebooks/reference/layout/lib.html">
     Libraries
    </a>
   </li>
   <li class="toctree-l2">
    <a class="reference internal" href="notebooks/reference/layout/scripts.html">
     Scripts
    </a>
   </li>
  </ul>
 </li>
 <li class="toctree-l1 has-children">
  <a class="reference internal" href="notebooks/reference/yaml.html">
   YAML Config Files
  </a>
  <input class="toctree-checkbox" id="toctree-checkbox-4" name="toctree-checkbox-4" type="checkbox"/>
  <label for="toctree-checkbox-4">
   <i class="fas fa-chevron-down">
   </i>
  </label>
  <ul>
   <li class="toctree-l2">
    <a class="reference internal" href="notebooks/reference/yaml/basics.html">
     YAML File Basics
    </a>
   </li>
   <li class="toctree-l2">
    <a class="reference internal" href="notebooks/reference/yaml/config_cam_baseline_example.html">
     Run Time - config_cam_baseline_example.yaml
    </a>
   </li>
   <li class="toctree-l2">
    <a class="reference internal" href="notebooks/reference/yaml/adf_variable_defaults.html">
     Variable Defaults - adf_variable_defaults.yaml
    </a>
   </li>
  </ul>
 </li>
 <li class="toctree-l1">
  <a class="reference internal" href="notebooks/reference/outputs.html">
   Explore ADF Outputs
  </a>
 </li>
</ul>
<p aria-level="2" class="caption" role="heading">
 <span class="caption-text">
  Troubleshooting
 </span>
</p>
<ul class="nav bd-sidenav">
 <li class="toctree-l1">
  <a class="reference internal" href="notebooks/help/contact.html">
   Contact
  </a>
 </li>
 <li class="toctree-l1">
  <a class="reference internal" href="notebooks/help/faqs.html">
   FAQS
  </a>
 </li>
</ul>

    </div>
</nav></div>
        <div class="bd-sidebar__bottom">
             <!-- To handle the deprecated key -->
            
        </div>
    </div>
    <div id="rtd-footer-container"></div>
</div>


          


          
<!-- A tiny helper pixel to detect if we've scrolled -->
<div class="sbt-scroll-pixel-helper"></div>
<!-- Main content -->
<div class="col py-0 content-container">
    
    <div class="header-article row sticky-top noprint">
        



<div class="col py-1 d-flex header-article-main">
    <div class="header-article__left">
        
        <label for="__navigation"
  class="headerbtn"
  data-toggle="tooltip"
data-placement="right"
title="Toggle navigation"
>
  

<span class="headerbtn__icon-container">
  <i class="fas fa-bars"></i>
  </span>

</label>

        
    </div>
    <div class="header-article__right">
<button onclick="toggleFullScreen()"
  class="headerbtn"
  data-toggle="tooltip"
data-placement="bottom"
title="Fullscreen mode"
>
  

<span class="headerbtn__icon-container">
  <i class="fas fa-expand"></i>
  </span>

</button>

<div class="menu-dropdown menu-dropdown-repository-buttons">
  <button class="headerbtn menu-dropdown__trigger"
      aria-label="Source repositories">
      <i class="fab fa-github"></i>
  </button>
  <div class="menu-dropdown__content">
    <ul>
      <li>
        <a href="https://github.com/justin-richling/ADF-Tutorial"
   class="headerbtn"
   data-toggle="tooltip"
data-placement="left"
title="Source repository"
>
  

<span class="headerbtn__icon-container">
  <i class="fab fa-github"></i>
  </span>
<span class="headerbtn__text-container">repository</span>
</a>

      </li>
      
      <li>
        <a href="https://github.com/justin-richling/ADF-Tutorial/issues/new?title=Issue%20on%20page%20%2FREADME.html&body=Your%20issue%20content%20here."
   class="headerbtn"
   data-toggle="tooltip"
data-placement="left"
title="Open an issue"
>
  

<span class="headerbtn__icon-container">
  <i class="fas fa-lightbulb"></i>
  </span>
<span class="headerbtn__text-container">open issue</span>
</a>

      </li>
      
    </ul>
  </div>
</div>

<div class="menu-dropdown menu-dropdown-download-buttons">
  <button class="headerbtn menu-dropdown__trigger"
      aria-label="Download this page">
      <i class="fas fa-download"></i>
  </button>
  <div class="menu-dropdown__content">
    <ul>
      <li>
        <a href="_sources/README.md"
   class="headerbtn"
   data-toggle="tooltip"
data-placement="left"
title="Download source file"
>
  

<span class="headerbtn__icon-container">
  <i class="fas fa-file"></i>
  </span>
<span class="headerbtn__text-container">.md</span>
</a>

      </li>
      
      <li>
        
<button onclick="printPdf(this)"
  class="headerbtn"
  data-toggle="tooltip"
data-placement="left"
title="Print to PDF"
>
  

<span class="headerbtn__icon-container">
  <i class="fas fa-file-pdf"></i>
  </span>
<span class="headerbtn__text-container">.pdf</span>
</button>

      </li>
      
    </ul>
  </div>
</div>

    </div>
</div>

<!-- Table of contents -->
<div class="col-md-3 bd-toc show noprint">
</div>
    </div>
    <div class="article row">
        <div class="col pl-md-3 pl-lg-5 content-container">
            <!-- Table of contents that is only displayed when printing the page -->
            <div id="jb-print-docs-body" class="onlyprint">
                <h1>Welcome to the AMWG Diagnostics Framework (ADF) Tutorial</h1>
                <!-- Table of contents -->
                <div id="print-main-content">
                    <div id="jb-print-toc">
                        
                    </div>
                </div>
            </div>
            <main id="main-content" role="main">
                
              <div>
                
  <section class="tex2jax_ignore mathjax_ignore" id="welcome-to-the-amwg-diagnostics-framework-adf-tutorial">
<h1>Welcome to the AMWG Diagnostics Framework (ADF) Tutorial<a class="headerlink" href="#welcome-to-the-amwg-diagnostics-framework-adf-tutorial" title="Permalink to this headline">#</a></h1>
<br>
<p><u><strong>AMWG</strong></u> <u><strong>D</strong></u>iagnostic <u><strong>F</strong></u>ramework</p>
<a class="reference internal image-reference" href="_images/adf_stats_image.jpeg"><img alt="_images/adf_stats_image.jpeg" class="align-center" src="_images/adf_stats_image.jpeg" style="height: 450px;" /></a>
<br>
<div></div>
<p>This material is based upon work supported by the National Center for Atmospheric Research (NCAR), which is a major facility sponsored by the National Science Foundation (NSF) under Cooperative Agreement No. 1852977. Staff time on this project was also supported by NSF award numbers <???></p>
<p>The purpose of this tutorial/demo is to give some exposure to what the ADF is and who could benefit from using it. This tutorial is meant to encompass basic to advanced uses of the ADF with guided examples to follow with sample data. It could also be used for a quick how-to for folks eager to try it with their own cases/data.</p>
<p>There is a reference section that is meant to dive deeper into aspects of the ADF as well as a Frequently Asked Questions (FAQ) section at the end for some common questions/gotchas that can arise from using the ADF.</p>
<h5>Feel free to explore this tutorial in any order you wish, but the layout is as follows:</h5>
<h3><u>Intro to the ADF</u></h3>
<ul class="simple">
<li><p>key features</p></li>
<li><p>ADF comparisons</p></li>
<li><p>ADF flow/layout/setup</p></li>
<li><p>outputs</p></li>
</ul>
<h3><u>Installation and Setup</u></h3>
<h3><u>ADF Basics</u></h3>
<ul>
<li><p>use of YAML files for configuring the ADF</p></li>
<li><p>quick ADF setup</p>
<p>-&gt; meant for users that have cases ready to try with the ADF</p>
</li>
</ul>
<h3><u>Guided Examples</u></h3>
<ul>
<li><p>step by step process to run the ADF with supplied sample data</p>
<p>-&gt; meant for users who want to try the ADF in an example setting</p>
</li>
<li><p>several examples including:</p>
<ul class="simple">
<li><p>CAM vs CAM diagnostics</p></li>
<li><p>CAM vs Obs/Reanalysis diagnostics</p></li>
<li><p>customization of conifguration files</p></li>
<li><p>adding custom plots to user’s ADF (advanced)</p></li>
<li><p>using the ADF in a Jupyter session (Jupyter Hub/Notebook)</p></li>
</ul>
</li>
</ul>
<h3><u>ADF Reference</u></h3>
<ul class="simple">
<li><p>detailed information about the ADF</p>
<ul>
<li><p>libraries and python scripts</p></li>
<li><p>config yaml files</p></li>
</ul>
</li>
</ul>
<h3><u>Troubleshooting the ADF</u></h3>
<ul class="simple">
<li><p>ways of contact for support</p></li>
<li><p>frequently asked questions</p></li>
</ul>
<div class="toctree-wrapper compound">
</div>
<div class="toctree-wrapper compound">
</div>
<div class="toctree-wrapper compound">
</div>
<div class="toctree-wrapper compound">
</div>
<div class="toctree-wrapper compound">
</div>
<div class="toctree-wrapper compound">
</div>
</section>

    <script type="text/x-thebe-config">
    {
        requestKernel: true,
        binderOptions: {
            repo: "binder-examples/jupyter-stacks-datascience",
            ref: "master",
        },
        codeMirrorConfig: {
            theme: "abcdef",
            mode: "python"
        },
        kernelOptions: {
            name: "python3",
            path: "./."
        },
        predefinedOutput: true
    }
    </script>
    <script>kernelName = 'python3'</script>

              </div>
              
            </main>
            <footer class="footer-article noprint">
                
    <!-- Previous / next buttons -->
<div class='prev-next-area'>
    <a class='right-next' id="next-link" href="notebooks/intro/adf.html" title="next page">
    <div class="prev-next-info">
        <p class="prev-next-subtitle">next</p>
        <p class="prev-next-title">Intro to the ADF</p>
    </div>
    <i class="fas fa-angle-right"></i>
    </a>
</div>
            </footer>
        </div>
    </div>
    <div class="footer-content row">
        <footer class="col footer"><p>
  
    By Justin Richling<br/>
  
</p>
        </footer>
    </div>
    
</div>


      </div>
    </div>
  
  <!-- Scripts loaded after <body> so the DOM is not blocked -->
  <script src="_static/scripts/pydata-sphinx-theme.js?digest=1999514e3f237ded88cf"></script>


  </body>
</html>