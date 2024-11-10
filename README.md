# Finite Volume Method for 1D CFD Simulations

## Overview

Welcome to the **Finite Volume Method for 1D CFD Simulations** repository! This collection of MATLAB scripts demonstrates various numerical techniques for solving **1D steady-state heat conduction and fluid flow problems** using the **Finite Volume Method (FVM)**. These simulations are fundamental in understanding heat transfer, diffusion, and fluid dynamics, making this repository a great resource for engineers, researchers, and students.

## Table of Contents
- [Introduction](#introduction)
- [Scientific Background](#scientific-background)
- [Repository Structure](#repository-structure)
- [Methodology](#methodology)
  - [Finite Volume Method (FVM)](#finite-volume-method-fvm)
  - [Numerical Schemes](#numerical-schemes)
- [Installation and Usage](#installation-and-usage)
- [Contributing](#contributing)
- [License](#license)

## Introduction

The **Finite Volume Method (FVM)** is a powerful numerical technique used in solving partial differential equations, especially in the fields of heat transfer and fluid dynamics. This repository contains MATLAB scripts to solve various 1D problems using FVM, such as:

- Heat diffusion with and without internal heat generation.
- Convective-diffusion problems using different discretization schemes.
- SIMPLE algorithm for fluid flow in 1D.
- Temperature distribution in extended surfaces (fins).

These simulations are designed to provide a clear understanding of how numerical methods are implemented in CFD analysis.

## Scientific Background

### What is the Finite Volume Method (FVM)?

The **Finite Volume Method** is a popular numerical technique used in the discretization of partial differential equations (PDEs). It is particularly effective for conservation laws, where the integral form of the equations ensures the conservation of quantities like mass, momentum, and energy.

### Numerical Schemes Used:

1. **Central Difference Scheme (CDS)**: A second-order accurate scheme used for problems dominated by diffusion.
2. **Upwind Scheme**: A first-order accurate scheme often used for convection-dominated problems to maintain numerical stability.
3. **SIMPLE Algorithm**: A widely-used iterative method to solve pressure-velocity coupling in fluid flow problems.

## Repository Structure

Here’s an overview of the MATLAB scripts included in this repository:

### 1. `oned_SIMPLE_Algorithm.m`
- **Description**: Implements the **SIMPLE (Semi-Implicit Method for Pressure Linked Equations)** algorithm for solving 1D fluid flow problems. This script solves for velocity and pressure fields in a 1D domain using iterative pressure correction.
- **Applications**: Useful for understanding fluid flow in ducts and channels where pressure-velocity coupling is significant.

### 2. `oned_diffusion.m`
- **Description**: Solves the **1D steady-state heat diffusion equation** without any convective effects or internal heat generation using the FVM approach.
- **Applications**: Ideal for understanding pure diffusion problems, such as heat conduction in a solid wall.

### 3. `Update oned_diffusion.m`
- **Description**: An updated version of the basic diffusion script with enhanced features and optimizations for faster convergence and better accuracy.
- **Applications**: Suitable for beginners looking to see iterative improvements in numerical methods.

### 4. `oned_diffusion_convection_CDS.m`
- **Description**: Solves the **1D convection-diffusion equation** using the **Central Difference Scheme (CDS)**. This method is more accurate for problems where diffusion dominates over convection.
- **Applications**: Useful for scenarios like pollutant dispersion in air or temperature distribution in low-speed fluid flows.

### 5. `oned_diffusion_convection_UpwindScheme.m`
- **Description**: Uses the **Upwind Scheme** to solve the **1D convection-diffusion equation**. This method is particularly stable for convection-dominated flows, reducing numerical oscillations.
- **Applications**: Ideal for simulations where fluid flow velocities are high, ensuring stability in the numerical solution.

### 6. `oned_diffusion_with_heat_generation.m`
- **Description**: Extends the heat diffusion problem to include **internal heat generation**. This script is useful for analyzing scenarios like nuclear fuel rods or chemical reactors where heat is generated within the material.
- **Applications**: Essential for understanding temperature profiles in systems with internal heat sources.

### 7. `oned_fin.m`
- **Description**: Solves the **temperature distribution in a 1D fin** using the FVM. The fin equation accounts for convective heat loss along its length, making it ideal for analyzing extended surfaces used to enhance heat transfer.
- **Applications**: Commonly used in the design of heat sinks, radiators, and cooling fins in electronics.

### 8. `LICENSE`
- **Description**: The license file outlines the terms under which this code can be used, modified, and distributed.

## Methodology

### Finite Volume Method (FVM)
The core idea of FVM is to divide the domain into small control volumes and apply conservation laws to each control volume. The values of fluxes at the control volume faces are approximated to solve for the unknown variables like temperature, pressure, or velocity.

### Numerical Schemes
- **Central Difference Scheme (CDS)**: Second-order accurate, best for diffusion-dominated problems.
- **Upwind Scheme**: First-order accurate, more stable for convection-dominated problems.
- **SIMPLE Algorithm**: Iterative method for solving coupled pressure-velocity equations in fluid flow.

## Installation and Usage

### Prerequisites
- MATLAB (version R2018b or later recommended)

### Running the Scripts
1. Clone the repository:
   ```bash
   git clone https://github.com/iammohith/Finite-Volume-Method-for-1D-CFD-Simulations.git
   ```
2. Navigate to the repository directory:
   ```bash
   cd Finite-Volume-Method-for-1D-CFD-Simulations
   ```
3. Open MATLAB and run any of the scripts, for example:
   ```matlab
   oned_diffusion.m
   ```
   
## Applications

This repository can be utilized in various fields, such as:
- Thermal analysis of materials
- Heat exchanger design
- Cooling of electronic components
- Flow simulation in ducts and channels
- Academic research and teaching numerical methods

## Contributing

We welcome contributions to enhance the functionality of this repository! Please feel free to submit pull requests, report bugs, or suggest new features.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
