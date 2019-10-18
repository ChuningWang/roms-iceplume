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
