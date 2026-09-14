# Weighted Sparsity Regularization for the Inverse EEG Problem

This repository contains code for numerical experiments on **weighted sparsity regularization for EEG source localization**, with a focus on dipole sources with **fixed and limited orientations**.

The implementations are based on the methods and numerical experiments presented in:

> O. L. Elvetun and N. Sudheer,
> *Weighted sparsity regularization for solving the inverse EEG problem: A case study*,
> Biomedical Signal Processing and Control, 107 (2025), 107673.
> https://doi.org/10.1016/j.bspc.2025.107673

## Overview

The EEG inverse problem aims to identify the brain sources that generate measured electrical potentials at the scalp. Because the problem is ill-posed and non-unique, regularization is used to obtain meaningful source estimates.

The approach implemented in this work uses **weighted sparsity regularization** to recover sparse dipole sources while addressing the depth bias associated with standard regularization methods.

The regularized inverse problem is formulated as

```math
\mathbf{x}_{\alpha} = \arg\min_{\mathbf{x}} \left\{ \frac{1}{2}\|\mathbf{L}\mathbf{x}-\mathbf{d}\|_2^2 + \alpha\|\mathbf{W}\mathbf{x}\|_1 \right\}
```

where $$\mathbf{L}$$ is the lead field matrix, $$\mathbf{d}$$ represents the EEG measurements, $$\mathbf{W}$$ is the weighting matrix, and $$\alpha$$ is the regularization parameter.

## Orientation Setups

### Fixed orientation

In the fixed-orientation setup, a single dipole is associated with each source position and is oriented **perpendicular to the cortical surface**.

The simulations use randomly selected source positions and generate noisy EEG measurements from sources located outside the reconstruction source space to avoid an inverse crime.

### Limited orientation

In the limited-orientation setup, the dipole orientation is allowed to vary within a **fixed angular deviation from the direction perpendicular to the cortical surface**.

Four basis orientations are constructed at each source position using the cortical-normal and two mutually orthogonal tangential directions. This provides a representation of dipoles whose orientations lie within the specified cone.

## Numerical Experiments

The code implements numerical experiments for Fixed orientation dipoles and Limited orientation dipoles using Weighted sparsity regularization formulation and hence
evaluating the source localization accuracy.

The experiments use simulated EEG data and lead-field matrices from the **ICBM-NY (New York Head) model**.

## Performance Measures

The reconstruction results can be evaluated using measures including:

* **Dipole Localization Error (DLE)** - Euclidean distance between the true source position and the position of the dominant estimated source.
* **Spatial Dispersion (SD)** – measures the spatial spread of the estimated inverse solution.
* **Earth Mover's Distance (EMD)** – measures the spatial discrepancy between the true and estimated source distributions.
* **Depth bias** – evaluates whether estimated source depth systematically differs from the true source depth.

## Software

The study uses the SEREEGA MATLAB toolbox for generating simulated EEG data and the New York Head lead field for the numerical experiments.

## Purpose of This Repository

This repository provides an example of the computational implementation of the numerical experiments associated with the published study. 
## Reference

Elvetun, O. L., & Sudheer, N. (2025).
**Weighted sparsity regularization for solving the inverse EEG problem: A case study.**
*Biomedical Signal Processing and Control, 107*, 107673.

DOI: https://doi.org/10.1016/j.bspc.2025.107673

## Author

**Niranjana Sudheer**
PhD in Applied Mathematics
Norwegian University of Life Sciences (NMBU)
