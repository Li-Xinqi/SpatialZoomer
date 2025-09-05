---
title: SpatialZoomer
layout: home
nav_order: 1
description: "SpatialZoomer uses spectral graph filter banks to extract multi-scale features from single-cell resolved spatial transcriptomics as a zoom-capable microscope."
---


# Overview
{: .fs-6 .fw-700 }

<p align="center">
  <img src="{{ '/assets/images/overview.jpg' | relative_url }}" alt="Overview" width="85%">
</p>

**SpatialZoomer** is a novel computational framework for **multi-scale feature analysis** of single-cell resolved spatial transcriptomics, functioning like a **zoom-capable microscope** that captures spatial structures from cells, niches, to tissue domains. It further **identifies critical scales automatically** by partitioning the cross-scale similarity map via dynamic programming. 

SpatialZoomer supports **diverse downstream applications**, such as identifying biologically meaningful signals, spatially dependent cell subclusters, and uncovering complex tissue architectures.

SpatialZoomer scales to datasets with millions of cells with high computational efficiency and low hardware requirements.

# Demo
{: .fs-6 .fw-700 }
<section id="demo">
Select a dataset to visualize the results of SpatialZoomer.
    <div class="dataset-buttons">
        <button type="button" name="button" class="btn btn-blue" data-dataset="1">Dataset 1</button>
        <button type="button" name="button" class="btn btn-blue" data-dataset="2">Dataset 2</button>
        <button type="button" name="button" class="btn btn-blue" data-dataset="3">Dataset 3</button>
        <button type="button" name="button" class="btn btn-blue" data-dataset="4">Dataset 4</button>
    </div>
    <div id="dataset-description">please select dataset</div>
    <div id="demo-content">
        <div class="images">
            <img id="img1" src="" alt="">
            <img id="img2" src="" alt="">
        </div>
        <img id="img3" src="" alt="">
        <p>Note: Across scales, cluster colors are aligned by an algorithm that maximizes the overlap between clusters.</p>
        <div class="params">
            <p>Please <strong>slide the scale axis</strong> to explore the results at critical scales (higher values result in more smoothed signals, identifying more global spatial structures):</p>
            <div class="param-group">
                <label for="scaleRange">Scale:</label>
                <input type="range" id="scaleRange" value="0">
                <span id="scaleValue"></span>
            </div>
            <p>Please <strong>choose the clustering resolution</strong> (a higher value results in more clusters):</p>
            <div class="param-group">
                <label for="resolutionSelect">Resolution:</label>
                <select id="resolutionSelect"></select>
            </div>
        </div>
    </div>
</section>

# Installation
{: .fs-6 .fw-700 }

### Option 1: Recommended (using `environment.yml`)
Download the `environment.yml` from [Github](https://github.com/Li-Xinqi/SpatialZoomer/blob/main/environment.yml).
```bash
# Create a virtual environment from the provided .yml file 
conda env create -f environment.yml
conda activate spatialzoomer
```


### Option 2: Manual installation
```bash
# Create and activate a virtual environment
conda create -n spatialzoomer python=3.10 -y
conda activate spatialzoomer

# Install required dependency
conda install -c conda-forge pyarrow

# Install SpatialZoomer
pip install SpatialZoomer==0.1.1
```

# Tutorial
{: .fs-6 .fw-700 }

To get started with SpatialZoomer, please refer to our tutorials available in the [GitHub repository](https://github.com/Li-Xinqi/SpatialZoomer.git).

# References
{: .fs-6 .fw-700 }

<script src="{{ '/assets/js/demo.js' | relative_url }}"></script>