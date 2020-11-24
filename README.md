# ROMS-ICEPLUME
Buoyant plume theory coupled with Rutgers ROMS.

This repository contains code of the Buoyant Plume Theory (BPT)/ROMS Coupled Model. Previously the model is developed based on Kate's ROMS branch (ROMS-ice) in order to allow coupling with sea-ice model.

However, Kate's ROMS is heavily modified and contains a lot of uniform configurations, which is potentially buggy. Since so far we haven't get the chance to couple the model with sea-ice, it is reasonable to migrate back to Rutgers ROMS to maintain a tidy version of the code.

# Theory
There is no way to explain the BPT in a few words in a Markdown file. A separate document is being constructed to fully explain the model.

# Example Case
An exmaple is provided for users to try out this coupled model. The files for this case is located in *./Iceplume\_Test*. It is a pure analytical case, thus no extra file is needed. The test domain is a 15x3 km channel, with an open boundary towards the east side. The first two rows of grid are masked to represent the glacier; subglacial discharge is injected from the center of the glacier into the channel. The channel is uniformly 400 m deep.

Initially the channel is stratified by a strong halocline at 50 m (set in roms\_fjord.in). Subglacial discharge of 200 m3/s is injected at 300 m depth, and a finite-line style plume of 260 m deep and 220 m long is used. Details of the plume parameters are set in *src/ana\_psource.h*.

To run this test case, set the right path in the build\_roms.bash file and compile. Different model options are provided in *src/iceplume\_test.h*; the default setup is the recommended options in terms of the ICEPLUME module.

# Key References
Cowton, T., Slater, D., Sole, A., Goldberg, D., & Nienow, P. (2015). Modeling the impact of glacial runoff on fjord circulation and submarine melt rate using a new subgrid‐scale parameterization for glacial plumes. Journal of Geophysical Research: Oceans, 120(2), 796-812.

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
