# Clinical Study Optimisation

**Optimising Recruitment and Sample Allocation in Multi-Site Rare Disease Studies**

This repository contains the full implementation of the optimisation framework described in the manuscript:

> *Optimising Recruitment and Sample Allocation in Multi-Site Rare Disease Studies*  
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
- 

---

## Requirements
- MATLAB (Only tested on 2024a version)
- COBRA Toolbox see: https://github.com/opencobra/cobratoolbox

---
## 🚀 Quick Start Guide

To reproduce the examples and figures from the paper, follow the steps below to set up the environment and run the live script.

### 1. 🖥️ Install MATLAB

Ensure you have MATLAB installed. This project has been tested with **MATLAB R2024a, R2024b, and R2025a**.

> 🔧 You can obtain MATLAB via your institution or from [mathworks.com](https://www.mathworks.com/).

---

### 2. 📦 Install the COBRA Toolbox

Clone and initialise the [COBRA Toolbox](https://github.com/opencobra/cobratoolbox):

```bash
git clone https://github.com/opencobra/cobratoolbox.git

In MATLAB:

```matlab
initCobraToolbox(false)  % Do not update
savepath
```
### 3. ⚙️ Install a Solver (e.g. Gurobi)

To enable optimisation, install a compatible solver such as **Gurobi**:

- Obtain a free academic license from [gurobi.com](https://www.gurobi.com/)
- Follow the installation instructions for your platform
- Use `changeCobraSolver` in MATLAB to set Gurobi as the default:

```matlab
changeCobraSolver('gurobi', 'all');
```

> 📌 Alternatively, you may use other supported solvers like CPLEX or GLPK.

### 4. 📥 Clone this Repository

Download or clone this repository:

```bash
git clone git@github.com:Digital-Metabolic-Twin-Centre/clinical-study-optimisation.git clinical_study_optimisation
```

Open MATLAB and navigate to the cloned folder:

```matlab
cd('clinical_study_optimisation/src')
```

---

### 5. ▶️ Run the Live Script

Launch the live script to run the full pipeline and reproduce results:

```matlab
open('recruitmentStrategyLiveScript.mlx')
```

Then click **Run All** to execute the entire analysis. This will generate all the results and figures presented in the paper.
