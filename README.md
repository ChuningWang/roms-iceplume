# ROMS-ICEPLUME

[![ver info](https://img.shields.io/badge/ROMS%20ver3.9-ICEPLUME%20ver1.1.1-brightgreen.svg)](https://github.com/ChuningWang/roms-iceplume)

Buoyant plume theory coupled with Rutgers ROMS.

This repository contains code of the Buoyant Plume Theory (BPT)/ROMS Coupled Model. Previously the model is developed based on Kate's ROMS branch (ROMS-ice) in order to allow coupling with sea-ice model.

However, Kate's ROMS is heavily modified and contains a lot of uniform configurations, which is potentially buggy. Since so far we haven't get the chance to couple the model with sea-ice, it is reasonable to migrate back to Rutgers ROMS to maintain a tidy version of the code.

This document is an introduction to the work in process ICEPLUME module for ROMS. It is modified from a similar package, IcePlume for the MITgcm, first developed by Dr. Tom Cowton. A detailed description of the MITgcm version is in [Cowton et al. (2015)](#key-references).

## Table of Contents

- [Overview](#overview)
- [Thoery](#theory)
- [Code Structure](#code-structure)
- [Example Case](#example-case)

## Overview

In ROMS, freshwater discharge is treated as point source going into a ROMS grid, which is activated by the compiling options **LuvSrc** or **LwSrc**. It reads in a total discharge (*Q<sub>bar</sub>*), a prescribed vertical weight function (*Q<sub>shape</sub>*), and tracer concentrations (T, S and passive), and volume/tracers are injected into a grid cell through horizontal advection (**LuvSrc**) or vertical convergence (**LwSrc**).

This method works well for shallow estuaries with a barotropic freshwater discharge forcing. In most cases, the freshwater is advected into the grid uniformly (*Q<sub>shape</sub>=1/N* or *Q<sub>shape</sub>=dz/H*. The expansion in volume and dilution of tracers drives a gravitational flow, which adjusts quickly to form an estuarine circulation.

However, the freshwater discharge is not always uniform in vertical direction. At high latitude, especially fjordic systems, summer meltwater can permeate the glacier through cracks and pores, and releases at the bottom of glacier head. If the glacier extends into the ocean, the freshwater discharge is injected at depth, which is called subglacial discharge.

In fjordic systems, the depth of subglacial discharge can be a few hundreds of meters. The runoff water is cold and fresh, rises and entrains the ambient water until it reaches a neutral buoyant layer. This process is non-hydrostatic, which cannot be accurately simulated by ROMS. Therefore, a parameterization for this process is required to correctly model subglacial discharge.

This document is a technical manual to introduce the ROMS-ICEPLUME coupled model, which uses a set of parameterizations to represent the subglacial discharge driven circulation in fjords. The parameterizations for entrainment, detrainment, background melting and coupler options are described in [the Theory Section](#theory); code structure is briefly summarized in [the Code Structure Section](#code-structure); steps for model usage are given in [the Example Case](#example-case).

![schem](readme_figs/schematics.png)

## Theory
### Buoyant Plume Thoery

The buoyant plume theory (BPT) is a set of equations that describes the development of buoyant plume rising near an ocean/glacier boundary. It is first described using a one-dimensional model by [Jenkins (1991, 2011)](#key-references), which assumes the plume initiates from a line source and only grows in direction normal to line source. Later this model is modified by [Cowton et al. (2015)](#key-references), which assumes the plume initiates from a point source instead of a line source, and grows uniformly in horizontal direction. In reality it is likely that the development of buoyant plume falls in between the two cases [(Jackson et al., 2017)](#key-references). In this section we attempt to use a ‘generalized’ buoyant plume model to describe both cases. This generalized model gives more flexibility in modeling different geometries of buoyant plume, and is also convenient for numerical applications. 

To summarize, the development of buoyant plume is controlled by the following BPT equations: 

![eq1][1]

where *A* is the plume cross section area, *u* is the velocity along plume transport axis, *g'=g(ρ<sub>a</sub>-ρ<sub>p</sub>)/ρ<sub>0</sub>* is the reduced gravitational acceleration. The temporal derivative of mass, ![][3] is the submarine melt rate from the glacier wall. *T*, *S* and *ρ* are temperature, salinity and density. The subscripts *a*, *p* and *b* indicate the plume model component of ambient water, plume water, and boundary layer, respectively. The ambient conditions are read from ROMS at each baroclinic time step. *Γ<sub>T</sub>* and *Γ<sub>S</sub>* are non-dimensional turbulent transfer coefficients, and *C<sub>d</sub>* is the ice-plume drag coefficient.

The parameter *α* is a constant entrainment rate. It determines the growth rate buoyant plume. It is the key parameter controlling the model behavior. Historically *α=0.1* is the conventional value [(Jenkins, 2011; Cowton et al., 2015)](#key-references), which is also verified by recent numerical studies [(Ezhova et al., 2018)]($key-references).

*L<sub>m</sub>* and *L<sub>c</sub>* are length of the cross section that is in contact with the glacier ice wall and the ambient water, respectively. In this model the plume can take any arbitrary shape as long as *L<sub>m</sub>* and *L<sub>c</sub>* can be determined. For example, when *L<sub>m</sub>=2b*, *L<sub>c</sub>=πb*, and *b=√(2A/π)*, The BPT Equations degenerate into the point source model [(Cowton et al., 2015)](#key-references), where *b* is the radius of the ‘half-cone’ plume. When *L<sub>m</sub>=L<sub>c</sub>=const* and *A=D∙const*, the BPT Equations degenerate into the line source model [(Jenkins, 2011)](#key-references), where *D* is the thickness of plume. For other plume geometries, *L<sub>m</sub>* and *L<sub>c</sub>* can be determined analytically or numerically. 

The BPT equations can be easily interpreted as the conservation of mass, momentum, heat and salt. In The first equation, *Au* is the volume flux through a horizontal transect; the first term on RHS is the entrainment of volume from the ambient water; the second term on RHS is the volume of melted water from glacier wall. In the second equation, *Au<sup>2</sup>* is the momentum flux through a transect; the first term on RHS is the buoyancy force; the second term on RHS is the friction (drag) from ambient water. In Equations third and fourth equations, the first and second terms on RHS are similar to that of the first equation; the third terms are turbulent mixing of temperature and salinity between plume and boundary layer. 

The melting of plume-ice boundary layer is described by a set of three equations [(Holland & Jenkins, 1999)](#key-references)

![eq2][2]

where *c<sub>i</sub>* and *c<sub>w</sub>* are heat capacity of ice and water, *L* is the latent heat of melting/freezing, *T<sub>i</sub>* is the temperature of glacier ice, *λ<sub>1</sub>*, *λ<sub>2</sub>*, and *λ<sub>3</sub>* are constants of linearity. The first equation is the heat budget during melting/freezing; LHS is the total heat absorbed by ice; RHS is the heat transfer from plume water to the glacier wall using a bulk flux formula. The second equation is similar to the first one but for salt budget. The last equation is the linear relationship between melting temperature as a function of salinity and depth.

A list of constants and parameter values are listed in the Table below. By solving the BPT and melting Equations, it gives the volume flux and tracer concentration of the buoyant plume (*A*, *u*, *T<sub>p</sub>*, *S<sub>p</sub>*), which is then coupled with the ocean model.

| Symbol            | Name                                   | Value    | Units                             |
| ----              | ----                                   | ----     | ----                              |
| *α*               | Entrainment Rate                       | 0.1      |                                   |
| *c<sub>i</sub>*   | Ice heat capacity                      | 2009     | J∙kg<sup>-1</sup>∙°C<sup>-1</sup> |
| *c<sub>w</sub>*   | Water heat capacity                    | 3974     | J∙kg<sup>-1</sup>∙°C<sup>-1</sup> |
| *ρ<sub>0</sub>*   | Reference ambient water density        | 1020     | kg∙m<sup>-3</sup>                 |
| *ρ<sub>ice</sub>* | Reference ice density                  | 916.7    | kg∙m<sup>-3</sup>                 |
| *L*               | Latent heat of melting                 | 335000   | J∙kg<sup>-1</sup>                 |
| *Γ<sub>T</sub>*   | Thermal turbulent transfer coefficient | 0.022    |                                   |
| *Γ<sub>S</sub>*   | Salt turbulent transfer coefficient    | 0.00062  |                                   |
| *C<sub>d</sub>*   | Ice/ocean drag coefficient             | 0.065    |                                   |
| *λ<sub>1</sub>*   | Freezing point salt slope              | -0.0573  |                                   |
| *λ<sub>2</sub>*   | Freezing point offset                  | 0.0832   |                                   |
| *λ<sub>3</sub>*   | Freezing point depth slope             | 0.000761 |                                   |
| *T<sub>i</sub>*   | Ice temperature                        | -10      | °C                                |
| *u<sub>bkg</sub>* | Minimum background velocity            | 0.3      | m∙s<sup>-1</sup>                  |

The background melt rate is calculated following the same melting Equations by substituting the plume temperature, salinity and velocity *T<sub>p</sub>*, *S<sub>p</sub>* with the ambient temperature and salinity *T<sub>p</sub>*, *S<sub>p</sub>*. The melt rate is then corrected for the surface area since portion of the grid cell is covered by the buoyant plume where the melting rate is calculated already.

This background melt parameterization is known to underestimate the melt rate by at least one order of magnitude [(D. A. Sutherland et al., 2019)](#key-references). The parameterization of [Holland and Jenkins (1999)](#key-references) is first developed to estimate melt rate below ice shelves, which are horizontally aligned on top of seawater. Since the meltwater is fresher than seawater below, it forms a hydrostatically stable layer, which tends to prevent further melting. Near a vertically aligned glacier/ocean boundary, melt water is hydrostatically unstable compared to ambient water, and convection cells tend to form which accelerates melting. This may partly explain the underestimated melt rate. A temporary solution is to increase the background velocity *u<sub>bkg</sub>*, which guarantees a minimum amount of melting; however, the melt rate produced is still significantly lower than field observations.

### Outflow Parameterization

To simplify the problem, the ambient stratification is represented with a two-layer setup, where *ρ<sub>1</sub><ρ<sub>2</sub>* are average densities of the upper and lower layer, respectively; *g'=g(ρ<sub>1</sub>-ρ<sub>2</sub>)/ρ<sub>ref</sub>* is defined as the reduced gravity between two layers; *ρ<sub>ref</sub>* is a reference density and *g* is the gravitational acceleration. Subglacial discharge plume rises along the glacier wall, during which the plume properties are predicted by the BPT. The plume density is *ρ<sub>p</sub>* when the rising stage terminates. If the value of *ρ<sub>p</sub>* falls between the densities of the two layers (*ρ<sub>1</sub><ρ<sub>p</sub><ρ<sub>2</sub>*), the outflowing plume forms near the density jump; otherwise the plume outflows at surface and travels downstream as gravity current. Assuming that the plume detrains as one uniform water mass, the nose speed *U<sub>D</sub>* of the outflowing current is estimated with an empirical parameterization developed by [Noh et al. (1992); Ching et al. (1993)](#key-references). For outflow in two-layer fluid, the nose velocity is dependent on a modified Richardson number ![][4], where $l_p$

## Example Case
An exmaple is provided for users to try out this coupled model. The files for this case is located in *./Iceplume\_Test*. It is a pure analytical case, thus no extra file is needed. The test domain is a 15x3 km channel, with an open boundary towards the east side. The first two rows of grid are masked to represent the glacier; subglacial discharge is injected from the center of the glacier into the channel. The channel is uniformly 400 m deep.

Initially the channel is stratified by a strong halocline at 50 m (set in roms\_fjord.in). Subglacial discharge of 200 m3/s is injected at 300 m depth, and a finite-line style plume of 260 m deep and 220 m long is used. Details of the plume parameters are set in *src/ana\_psource.h*.

To run this test case, set the right path in the build\_roms.bash file and compile. Different model options are provided in *src/iceplume\_test.h*; the default setup is the recommended options in terms of the ICEPLUME module.

---

## Key References

[Ching, C. Y., Fernando, H. J. S., & Noh, Y. (1993). Interaction of a Negatively Buoyant Line Plume with a Density Interface. Dynamics of Atmospheres and Oceans, 19(1-4), 367-388.](https://www.sciencedirect.com/science/article/abs/pii/0377026593900426)

[Cowton, T., Slater, D., Sole, A., Goldberg, D., & Nienow, P. (2015). Modeling the impact of glacial runoff on fjord circulation and submarine melt rate using a new subgrid‐scale parameterization for glacial plumes. Journal of Geophysical Research: Oceans, 120(2), 796-812.](https://agupubs.onlinelibrary.wiley.com/doi/10.1002/2014JC010324)

[Sutherland, D. A., Jackson, R. H., Kienholz, C., Amundson, J. M., Dryer, W. P., Duncan, D., et al. (2019). Direct observations of submarine melt and subsurface geometry at a tidewater glacier. Science, 365(6451), 369.](https://science.sciencemag.org/content/365/6451/369.full)

[Ezhova, E., Cenedese, C., & Brandt, L. (2018). Dynamics of Three-Dimensional Turbulent Wall Plumes and Implications for Estimates of Submarine Glacier Melting. Journal of Physical Oceanography, 48(9), 1941-1950.](https://journals.ametsoc.org/view/journals/phoc/48/9/jpo-d-17-0194.1.xml)

[Holland, D. M., & Jenkins, A. (1999). Modeling thermodynamic ice-ocean interactions at the base of an ice shelf. Journal of Physical Oceanography, 29(8), 1787-1800.](https://journals.ametsoc.org/view/journals/phoc/29/8/1520-0485_1999_029_1787_mtioia_2.0.co_2.xml?tab_body=fulltext-display)

[Jenkins, A. (1991). A One-Dimensional Model of Ice Shelf-Ocean Interaction. Journal of Geophysical Research-Oceans, 96(C11), 20671-20677.](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/91JC01842)

[Jenkins, A. (2011). Convection-Driven Melting near the Grounding Lines of Ice Shelves and Tidewater Glaciers. Journal of Physical Oceanography, 41(12), 2279-2294.](https://journals.ametsoc.org/view/journals/phoc/41/12/jpo-d-11-03.1.xml)

[Noh, Y., Fernando, H. J. S., & Ching, C. Y. (1992). Flows Induced by the Impingement of a 2-Dimensional Thermal on a Density Interface. Journal of Physical Oceanography, 22(10), 1207-1220.](https://journals.ametsoc.org/view/journals/phoc/22/10/1520-0485_1992_022_1207_fibtio_2_0_co_2.xml)

---

## Contact Info

Chuning Wang, School of Oceanography, Shanghai Jiaotong University

wangchuning@sjtu.edu.cn

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



[1]: https://latex.codecogs.com/svg.image?\begin{aligned}\frac{d}{dz}[Au]&space;&&space;=\alpha&space;L_c&space;u&space;&plus;&space;L_m\dot{m}&space;\\\\\frac{d}{dz}[Au^2]&space;&&space;=g'A&space;&plus;&space;L_m&space;C_d&space;u^2&space;\\\\\frac{d}{dz}[AuT_p]&space;&&space;=\alpha&space;L_c&space;u&space;T_a&space;&plus;&space;L_m\dot{m}T_b&space;-&space;L_m\Gamma_T&space;C_d^{1/2}u(T_p-T_b)&space;\\\\\frac{d}{dz}[AuS_p]&space;&&space;=\alpha&space;L_c&space;u&space;S_a&space;&plus;&space;L_m\dot{m}S_b&space;-&space;L_m\Gamma_S&space;C_d^{1/2}u(S_p-S_b)\end{aligned}
[2]: https://latex.codecogs.com/svg.image?\begin{aligned}\dot{m}(c_i(T_b-T_i)&plus;L)&space;&&space;=c_w\Gamma_T&space;C_d^{1/2}u(T_p-T_b)\\\\\dot{m}S_b&space;&&space;=\Gamma_S&space;C_d^{1/2}u(S_p-S_b)\\\\T_b&space;&&space;=\lambda_1&space;S_b&space;&plus;&space;\lambda_2&space;&plus;&space;\lambda_3&space;z\end{aligned}
[3]: https://latex.codecogs.com/svg.image?\dot{m}
[4]: https://latex.codecogs.com/svg.image?Ri=\frac{g'l_P}{W_P^2}
