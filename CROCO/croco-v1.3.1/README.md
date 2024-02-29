CROCO v1.3
Released date : 28 November 2022
Previous release : CROCO v1.2.1 (March 2022)

Reminders:
CROCO sources and CROCO_TOOLS (the follow-on of ROMS_TOOLS) are now distributed separately (for croco_tools releases, see associated tab at https://www.croco-ocean.org/download/croco-project ). ROMS_AGRIF is not maintained anymore and we strongly encourage ROMS_AGRIF users to switch to CROCO. CROCO version available directly from the git repository is an unstable development version. Standard users should use the stable one downloaded from the web site.
CROCO documentation is available at https://croco-ocean.gitlabpages.inria.fr/croco_doc
You can also subscribe at CROCO forum at https://forum.croco-ocean.org

New in CROCO v1.3
CROCO v1.3 is recommended for users of PISCES, MUSTANG and the coupling toolbox. In addition it includes correction of critical bugs, inconsistencies and light improvements

Environment

Little updates on create_config.bash : mostly comments and some little bugfixes


Coupling
full description : https://croco-ocean.gitlabpages.inria.fr/croco_doc/model/model.coupling.html and https://croco-ocean.gitlabpages.inria.fr/croco_doc/tutos/tutos.16.coupling.html

put by default exchanges of AVERAGE fields over the coupling timestep
add possibility to use coupling with WRF moving nests (coupling only on the parent domain, and interpolation in WRF between the parent and moving nest for ocean surface conditions => requires to use the last version of WRF in the WRF-CROCO github fork)
some corrections for coupling with a toy model
add the possibility to create oasis restart files from pre-existing model outputs
coupling scripts: correct some typos, comments, and directory and grid names for oasis files and oce_compile for AROME and ECMWF


Sediment
full description : https://croco-ocean.gitlabpages.inria.fr/croco_doc/model/model.modules.sediment.html

Test and reorganise flocculation capabilities. FLOCMOD module is present in MUSTANG and in SEDIMENT with equal routines. To avoid redundancies and clarify maintenance/dev on this module, an independent module has been be created, flocmod.F90 with interfaces in MUSTANG/SEDIMENT


Numerics and parametrisations
full description : https://croco-ocean.gitlabpages.inria.fr/croco_doc/model/model.numerics.html

UV_HADV_UP5: change to apply UP5 to both predictor and corrector time steps to avoid checkboard noise detected in regional simulations (thanks to J. Gula)
GLS_MIXING: add limitation in potential flow region following Larsen and Fuhrman (2018), to prevent the instable growth of turbulent viscosity for non-breaking external or internal waves


New Configurations
full description : https://croco-ocean.gitlabpages.inria.fr/croco_doc/model/model.test_cases.html

1 new test cases for sediment processes (can be used with USGS or MUSTANG sediment model) :a schematized estuary (from MARS corresponding test case). From TIDAL_FLAT, the grid is set to 200x90 then the cells sizes are fixed to dx = 600m and dy = 100m. The OBC West is the same as in TIDAL_FLAT. The bathymetry is derived from MARS code with the same parameters. A river is placed on the East side with sediment input (400m3/s - 50mg/L of MUD)


Pisces

fix various bugs and improve robustness


Runoff
Possibility to add a source point as a volume vertical influx. With this feature there is no need to take care of the position of the source with the mask and a source can be added anywhere on the grid. The outflow is applied at rho point (cpp key PSOURCE_MASS)

Miscellaneous

run_croco scripts: some corrections for adequation with prod env, comments and re-organization
NC4 for diagnostics
Some fixes when writing into existing diagnostics files (diag TS, diag Momentum, Vorticity, PV, EK, PV, VRT...)
Various fixes when using NC4PAR
Move Surface atmospheric forcing and Lateral forcing cpp-keys upper in cppdefs.h
Output in double precisions in case of PARALLEL_FILES (to be compatible with ncjoin)
Put back by default previous partit.F to handle with FRC_FRY
Changes in compilations options for ifort and gfortran


CROCO_TOOLS

Improve octave compatibility
Add global attribute in netcdf files (notably: origin date)
Tides: add the possibility to create frc files with tides only, and add a routine to create tide-only frc files for interanual simulation
Add pre-processing scripts for downlading glorys oceanic reanalysis data.
Coupling: add scripts to create weight files and smoothing for the coupler (interpolation between model grids), some corrections in pre-processing for WWIII for adequations with more recent versions of WWIII, and some other minor corrections in scripts
Biology: Pre-processing scripts have been moved to Preprocessing_tools/Bio. Update some scripts to handle PISCES-quota
Nesting_tools: some small corrections
