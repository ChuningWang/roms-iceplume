# ROMS-ICEPLUME

[![ver info](https://img.shields.io/badge/ROMS%20ver3.9-ICEPLUME%20ver1.1.1-brightgreen.svg)](https://github.com/ChuningWang/roms-iceplume)

Buoyant plume theory coupled with Rutgers ROMS.

This repository contains code of the Buoyant Plume Theory (BPT)/ROMS Coupled Model. Previously the model is developed based on Kate's ROMS branch (ROMS-ice) in order to allow coupling with sea-ice model.

However, Kate's ROMS is heavily modified and contains a lot of uniform configurations, which is potentially buggy. Since so far we haven't get the chance to couple the model with sea-ice, it is reasonable to migrate back to Rutgers ROMS to maintain a tidy version of the code.

This document is an introduction to the work in process ICEPLUME module for ROMS. It is modified from a similar package, IcePlume for the MITgcm, first developed by Dr. Tom Cowton. A detailed description of the MITgcm version is in [Cowton et al. 2015](#key-references)

## Table of Contents

- [Overview](#overview)
- [Thoery](#theory)
- [Code Structure](#code-structure)
- [Example Case](#example-case)

## Overview

In ROMS, freshwater discharge is treated as point source going into a ROMS grid, which is activated by the compiling options **LuvSrc** or **LwSrc**. It reads in a total discharge (*Q_<sub>bar</sub>*), a prescribed vertical weight function (*Q_<sub>shape</sub>*), and tracer concentrations (T, S and passive), and volume/tracers are injected into a grid cell through horizontal advection (**LuvSrc**) or vertical convergence (**LwSrc**).

This method works well for shallow estuaries with a barotropic freshwater discharge forcing. In most cases, the freshwater is advected into the grid uniformly (*Q_<sub>shape</sub>=1/N* or *Q_<sub>shape</sub>=dz/H*. The expansion in volume and dilution of tracers drives a gravitational flow, which adjusts quickly to form an estuarine circulation.

However, the freshwater discharge is not always uniform in vertical direction. At high latitude, especially fjordic systems, summer meltwater can permeate the glacier through cracks and pores, and releases at the bottom of glacier head. If the glacier extends into the ocean, the freshwater discharge is injected at depth, which is called subglacial discharge.

In fjordic systems, the depth of subglacial discharge can be a few hundreds of meters. The runoff water is cold and fresh, rises and entrains the ambient water until it reaches a neutral buoyant layer. This process is non-hydrostatic, which cannot be accurately simulated by ROMS. Therefore, a parameterization for this process is required to correctly model subglacial discharge.

This document is a technical manual to introduce the ROMS-ICEPLUME coupled model, which uses a set of parameterizations to represent the subglacial discharge driven circulation in fjords. The parameterizations for entrainment, detrainment, background melting and coupler options are described in [the Theory Section](#theory); code structure is briefly summarized in [the Code Structure Section](#code-structure); steps for model usage are given in [the Example Case](#example-case).

![schem](readme_figs/schematics.png)

## Theory
There is no way to explain the BPT in a few words in a Markdown file. A separate document is being constructed to fully explain the model.

## Example Case
An exmaple is provided for users to try out this coupled model. The files for this case is located in *./Iceplume\_Test*. It is a pure analytical case, thus no extra file is needed. The test domain is a 15x3 km channel, with an open boundary towards the east side. The first two rows of grid are masked to represent the glacier; subglacial discharge is injected from the center of the glacier into the channel. The channel is uniformly 400 m deep.

Initially the channel is stratified by a strong halocline at 50 m (set in roms\_fjord.in). Subglacial discharge of 200 m3/s is injected at 300 m depth, and a finite-line style plume of 260 m deep and 220 m long is used. Details of the plume parameters are set in *src/ana\_psource.h*.

To run this test case, set the right path in the build\_roms.bash file and compile. Different model options are provided in *src/iceplume\_test.h*; the default setup is the recommended options in terms of the ICEPLUME module.

---

## Key References
[Cowton, T., Slater, D., Sole, A., Goldberg, D., & Nienow, P. (2015). Modeling the impact of glacial runoff on fjord circulation and submarine melt rate using a new subgrid‐scale parameterization for glacial plumes. Journal of Geophysical Research: Oceans, 120(2), 796-812.](https://agupubs.onlinelibrary.wiley.com/doi/10.1002/2014JC010324)

---

Chuning Wang

chuning@marine.rutgers.edu

---

## Ver 1.1.0
This is a major update. In this update I attempt to fix and improve the **ICEPLUME_DET_AVERAGE** method. The major issue I had is to average data across tiles. Depending on the MPI method, accessing data from another tile may not be allowed. This happened when I run the code on a super computer with Intel MPI. Historically it has also caused issues for other uses; in my own applications I 'fix' it by manually set averaging span within a single tile.

However this simple 'patch' is not a good solution. It is necessary to design the code so it can either

- Automatically detect the tile boundary and restrict the average span within it; or

- Accessing data across tiles for the spatial averaging.

Therefore, two methods are applied to fix the same problem. The first solution is to iteratively find the boundary of each tile, and compare the averaging span with tile boundaries; this is easy to achieve, however, it may generate slightly difference solution when the domain is tiled differently. The second solution, which is more aggressive, uses ROMS internal function *mp_aggregate2d* and *mp_aggregate3d* to collect data from each tile and then apply the average. It guarantees the same averaging span each time, but exchanging data between tiles can be numerically expensive when the tile number is large.

In this new version, the default averaging method is the the first one, which thoeretically is 'faster'. To use the more aggressive approach, user needs to manually activate it with the cpp flag **ICEPLUME_CTIL_AVG** as well as the original flag **ICEPLUME_DET_AVERAGE**; **CTIL** stands for *cross tile*. User should run benchmark tests to decide which strategy works best for the machine.

Another improvement is that the spatial averaging calculation now skips land-masked points, and weight the average with respect to each grid size, which aims to serve realistic topography in future applications. Other modifications include but not limited to data broadcasting, tile exchange, etc.

Chuning Wang

2020-11-24

---

## Ver 1.1.1
More optimization on the averaging algorithm. Using ROMS internal function *mp_assemble* to substitute *mp_aggregate* in the coupling function. This is to avoid fetching global arrays of density at each timestep. Instead, each tile calculates the weighted SUM, and *mp_assemble* is used to get and broadcast the reduced sum from and to each tile; then the weighted average is calculated using the reduced sum.

Since this is a better scheme over the method introduced in [Ver 1.1.0](##Ver-1.1.0), I decide to remove the temporarily used CPP flag **ICEPLUME_CTIL_AVG**. From now on only the new default averaging method is available.

Other modifications - Previously **ICEPLUME_DET_AVERAGE** averages density profiles on \sigma surfaces. This is based on the assumption that the depth near glacier grounding is somewhat constant, which is fine for idealized simulation. To generalize for realistic applications, in the current version the density profile is averaged on z surfaces. This is realized by linearly interpolating the density profile from each grid.

Chuning Wang

2020-12-03
