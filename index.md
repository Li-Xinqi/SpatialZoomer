---
title: SpatialZoomer
layout: home
nav_order: 1
description: "SpatialZoomer uses spectral graph filter banks and dynamic programming to extract multi-scale features from spatial transcriptomics data."
---


# Overview
{: .fs-6 .fw-700 }
<div class="overview-container">
    <div class="overview-image">
    <img src="{{ '/assets/images/overview.png' | relative_url }}" alt="Overview">
    </div>
    <div class="overview-text">
    <p>SpatialZoomer uses spectral graph filter banks and dynamic programming to extract multi-scale features from spatial transcriptomics data.</p>
    </div>
</div>

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
        <div class="params">
            <div class="param-group">
                <label for="scaleRange">Scale:</label>
                <input type="range" id="scaleRange" min="0" max="6" step="1" value="0">
                <span id="scaleValue"></span>
            </div>
            <div class="param-group">
                <label for="resolutionSelect">Resolution:</label>
                <select id="resolutionSelect"></select>
            </div>
        </div>
    </div>
</section>

# Installation
{: .fs-6 .fw-700 }
```python
pip install spatialzoomer
```

# Tutorial
{: .fs-6 .fw-700 }
To get started with SpatialZoomer, check out our:

[tutorial 1](#)

[tutorial 1](#)

# References
{: .fs-6 .fw-700 }
References


<!-- <script>
// 注入 Jekyll baseurl
window.jekyllBaseurl = "{{ site.baseurl }}";
</script> -->
<script src="{{ '/assets/js/demo.js' | relative_url }}"></script>