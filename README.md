# Clinical Study Optimisation

**Cardinality-Constrained Modelling Framework for Transparent and Fair Multi-Site Recruitment in Rare Disease Research**

This repository contains the full implementation of the optimisation framework described in the manuscript:

> *Cardinality-Constrained Modelling for Transparent and Fair Multi-Site Recruitment in Rare Disease Research: Application to Recon4IMD*  
> Farid Zare¹², Ronan M.T. Fleming¹²  
> ¹School of Medicine, University of Galway, Ireland  
> ²Digital Metabolic Twin Centre, University of Galway, Ireland  

---

## 🔍 Overview

Clinical research involving rare diseases often requires multi-site collaboration, due to low disease prevalence and logistical complexity. This repository provides a unified optimisation framework for:

- Patient recruitment across multiple healthcare providers
- Sample allocation to analytical laboratories
- Platform selection based on biomarker coverage and cost
- Minimising logistics burden through local sample processing
- Incorporating clinical interest and diagnostic confidence

The framework uses **cardinality-constrained optimisation** to ensure transparent, fair, and cost-efficient clinical study design.

---

## 📁 Repository Contents
- data (Containing actual and noised data)
- graphics
- results
- An end-to-end MATLAB live script which produces all the examples and figures in the paper

---

## Requirements
- MATLAB (Only tested on 2024a version)
- COBRA Toolbox see: https://github.com/opencobra/cobratoolbox
