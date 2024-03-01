# Changelog
## [v.1.3.1] - 2023-09-21
### Added

- COUPLING
  - Add EXACT_RESTART + OA/OW_COUPLING capability : the coupler initialisation cpl_prism_define.F was wrongly done twice in this case (line 1147 in get_initial.F). It was leading to a coupler crash during the initilization phase.
  - Add the capability to cutoff the heat fluxes in case of water under the freezing point (-1.8°C) => Simplest sea ice correction to T,S fluxes.
### Fixed
- OW_COUPLING: missing if defined OW_COUPLING in forces.h for smstr (that can be updated in case of OW COUPLING). Fix issue #127

- USE_CALENDAR : 
  - Compatibility with XIOS
  - Change format of date format : 
    - in croco.in : from dd/mm/yyyy to yyyy-mm-dd
    - read origin_date in time netcdf attribute for netCDF input files.


-DIAGNOSTICS_BIO : 
  - fix netcdf parralel (NC4PAR) netcdf files writing with DIAGNOSTICS_BIO (bgc fluxes)
  - fix time and scrum_time in avergae files of croco_bio_diags_avg.nc

- WET_DRY : solve blowup problem with dry cells at boundaries

- AGRIF Compilation : fix jobcomp.bash with gcc > 10

- PISCES : fix missing local variables initialization
### Changed

### Deprecated

### Removed

### Other

## [v.1.3] - 2022-11-28
