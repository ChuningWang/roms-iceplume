# $Id: Makefile 1483 2014-03-15 17:05:10Z rblod $
#
#======================================================================
# CROCO is a branch of ROMS developped at IRD and INRIA, in France
# The two other branches from UCLA (Shchepetkin et al)
# and Rutgers University (Arango et al) are under MIT/X style license.
# CROCO specific routines (nesting) are under CeCILL-C license.
#
# CROCO website : http://www.croco-ocean.org
#======================================================================
#
# Universal machine independent makefile for ROMS model
#
#======================================================================
# Set machine dependent definitions and rules.
#======================================================================

include Makedefs
#======================================================================
# Model Configuration:
#======================================================================
# SRCS: source codes files are sorted into groups, separated by
# blanc lines:
#   1) main driving part;		 2) 2D time stepping engine;
#   3) 3D time stepping engine;		 4) non-hydrostatic engine
#   5) sea-water EOS;			 6) vertical mixing schemes;
#   7) on-fly model diagnostics;	 8) netCDF I/O routines;
#   9) main forcing;			10) Wave forcing routines;
#  11) Online surface forcing;		12) Floats routines;
#  13) Station diagnostics;		14) biology/sediment routines;
#  15) PISCES biology;			16) MPI routines;
#  17) AGRIF routines;			18) OASIS coupling interface;
#
# SRC90: additional F90 routines for PISCES code and Non-hydrostatic
#        analysis routines
#======================================================================

 SRCS = main.F		step.F		read_inp.F\
	timers_roms.F	init_scalars.F	init_arrays.F	set_weights.F\
	set_scoord.F	ana_grid.F	setup_grid1.F	setup_grid2.F\
	set_nudgcof.F	ana_initial.F	analytical.F	zonavg.F\
\
	step2d.F 	u2dbc.F		v2dbc.F		zetabc.F\
	obc_volcons.F\
\
	pre_step3d.F	step3d_t.F	step3d_uv1.F	step3d_uv2.F\
	prsgrd.F	rhs3d.F		set_depth.F	omega.F\
	uv3dmix.F	uv3dmix_spg.F	t3dmix.F	t3dmix_spg.F\
	hmix_coef.F	wetdry.F	wavedry.F\
	u3dbc.F		v3dbc.F		t3dbc.F\
\
	step3d_fast.F\
	step3d_w.F	rhs3d_w_nh.F	initial_nbq.F	grid_nbq.F\
	unbq_bc.F	vnbq_bc.F	wnbq_bc.F	rnbq_bc.F\
        w3dbc.F		nbq_bry_store.F\
\
	rho_eos.F	ab_ratio.F	alfabeta.F      alfabeta_k.F\
\
	ana_vmix.F	bvf_mix.F	lmd_vmix.F 	gls_mixing.F\
	lmd_skpp.F	lmd_bkpp.F	lmd_swfrac.F	lmd_wscale.F\
\
	diag.F		wvlcty.F	checkdims.F	grid_stiffness.F\
	bio_diag.F	setup_kwds.F    check_kwds.F	check_srcs.F\
	check_switches1.F		check_switches2.F\
	debug.F\
\
        param.F         ncscrum.F   scalars.F\
	output.F	put_global_atts.F\
	nf_fread.F	nf_fread_x.F	nf_fread_y.F	nf_read_bry.F\
	get_date.F	lenstr.F	closecdf.F	insert_node.F\
	fillvalue.F	nf_add_attribute.F		set_cycle.F\
	def_grid_2d.F   def_grid_3d.F   def_his.F       def_rst.F\
	def_diags.F	def_diagsM.F	def_bio_diags.F\
	wrt_grid.F      wrt_his.F       wrt_avg.F	wrt_rst.F\
	wrt_diags.F	wrt_diags_avg.F	wrt_diagsM.F	wrt_diagsM_avg.F\
	wrt_bio_diags.F	wrt_bio_diags_avg.F\
	set_avg.F	set_diags_avg.F	set_diagsM_avg.F\
	set_bio_diags_avg.F\
	def_diags_vrt.F	wrt_diags_vrt.F\
	set_diags_vrt.F	set_diags_vrt_avg.F	wrt_diags_vrt_avg.F\
	def_diags_ek.F	wrt_diags_ek.F\
	set_diags_ek.F	set_diags_ek_avg.F	wrt_diags_ek_avg.F\
	def_diags_pv.F	wrt_diags_pv.F\
	set_diags_pv.F	set_diags_pv_avg.F	wrt_diags_pv_avg.F\
	def_diags_eddy.F\
	set_diags_eddy_avg.F	wrt_diags_eddy_avg.F\
	def_surf.F	wrt_surf.F\
	set_surf_avg.F	wrt_surf_avg.F\
\
	get_grid.F	get_initial.F	get_vbc.F	get_wwave.F\
	get_tclima.F    get_uclima.F    get_ssh.F       get_sss.F\
	get_smflux.F    get_stflux.F    get_srflux.F    get_sst.F\
	mod_tides_mas.F tidedata.F      mas.F         get_tides.F\
        clm_tides.F     get_bulk.F      bulk_flux.F\
	get_bry.F       get_bry_bio.F	sstskin.F\
	get_psource.F   get_psource_ts.F get_btflux.F\
\
	mrl_wci.F  	wkb_wwave.F 	wkbbc.F		get_bry_wkb.F\
\
	online_bulk_var.F		online_get_bulk.F\
	online_interp.F			online_interpolate_bulk.F\
	online_set_bulk.F\
\
	init_floats.F	wrt_floats.F	step_floats.F	rhs_floats.F\
	interp_rho.F	def_floats.F	init_arrays_floats.F\
	random_walk.F	get_initial_floats.F\
\
	init_sta.F	wrt_sta.F	step_sta.F	interp_sta.F\
	def_sta.F	init_arrays_sta.F\
\
	biology.F	o2sato.F	sediment.F	bbl.F\
\
	MPI_Setup.F	MessPass2D.F	MessPass3D.F	exchange.F\
	autotiling.F\
\
	zoom.F		update2D.F	set_nudgcof_fine.F\
	zoombc_2D.F	zoombc_3D.F	uv3dpremix.F\
	t3dpremix.F     update3D.F	zoombc_3Dfast.F\
	Agrif2Model.F\
\
	send_xios_diags.F\
\
	cpl_prism_define.F	cpl_prism_put.F 	cpl_prism_init.F\
	cpl_prism_get.F 	cpl_prism_getvar.F      cpl_prism_grid.F\
\
        module_substance.F  module_MUSTANG.F\
\
        online_spectral_diags.F var3d_oa_out.F var2d_oa_out.F\
        scal0d_oa_out_loc.F scal0d_oa_out_full.F\
        output_oa.F



 SRCS90 = \
        par_pisces.F90	ocean2pisces.F90	trc.F90		sms_pisces.F90\
	p4zche.F90	p4zint.F90		p4zlys.F90	p4zflx.F90\
	p4zlim.F90	p4zsink.F90		p4zmicro.F90	p4zmeso.F90\
	p4zmort.F90	p4zopt.F90		p4zprod.F90	p4zrem.F90\
	p4zsbc.F90      p4zsed.F90	        trcwri_pisces.F90\
	p4zagg.F90      p4zfechem.F90           p4zligand.F90   p4zpoc.F90\
	p4zsms.F90      p5zlim.F90              p5zmeso.F90\
	p5zmicro.F90    p5zmort.F90             p5zprod.F90\
	p4zbio.F90      trcsms_pisces.F90	trcini_pisces.F90\
	pisces_ini.F90  oce_sed.F90             par_sed.F90     sedadv.F90\
	sedarr.F90      sedbtb.F90              sedchem.F90     sedco3.F90\
	seddiff.F90     seddsr.F90              seddta.F90      sed.F90\
	sedini.F90      sedinitrc.F90           sedinorg.F90     sedmat.F90\
	sedmbc.F90      sedmodel.F90            sed_oce.F90     sedorg.F90\
	sedrst.F90      sedsfc.F90              sedstp.F90      sedwri.F90\
	setavg_sed.F90\
	module_parameter_oa.F90 module_grd_oa.F90\
	module_oa_time.F90	module_oa_space.F90 module_oa_periode.F90\
	module_oa_variables.F90	module_oa_type.F90 module_tile_oa.F90\
        module_oa_level.F90 module_oa_interface.F90 var_oa.F90\
	tooldatosec.F90	toolsectodat.F90 tooldecompdat.F90 tooldatetosec.F90 \
        toolorigindate.F90 \
        flocmod.F90   comsubstance.F90        substance.F90\
        comMUSTANG.F90  coupler_MUSTANG.F90 sed_MUSTANG_CROCO.F90\
        sed_MUSTANG.F90 initMUSTANG.F90 plug_MUSTANG_CROCO.F90


AMRDIR = AGRIF/AGRIF_YOURFILES

#======================================================================

 TOOLS = mpc cross_matrix cppcheck srcscheck checkkwds partit ncjoin ncrename


 RCS = $(SRCS:.F=.f)

 RCS90 = $(SRCS90:.F90=.f90)

 OBJS = $(RCS:.f=.o)
 OBJS90 = $(RCS90:.f90=.o)

 SBIN = croco

 AMRRCS=$(AMRSRCS:.F=.f)

 AMROBJS=$(AMRRCS:.f=.o)


%.mod : %.o
	@touch $@


#
# Everything
# ==========
all: tools depend $(SBIN)

#
# Executable file.
# ========= =====
#
$(SBIN):   $(OBJS90) $(OBJS)
	 $(LDR) $(FFLAGS) $(LDFLAGS) -o a.out $(OBJS90) $(OBJS) $(LCDF) $(LMPI)
	 mv a.out $(SBIN)

#
# Auxiliary utility programs and List of Dependecies:
# ========= ======= ======== === ==== == ============
#

tools: $(TOOLS)

mpc: mpc.F
	$(CPP) -P $(CPPFLAGS) mpc.F > mpc_.f
	$(LDR) $(FFLAGS) $(LDFLAGS) -o mpc mpc_.f

cross_matrix: mpc cross_matrix.F
	$(CPP) -P $(CPPFLAGS) cross_matrix.F | ./mpc > cross_matrix_.f
	$(CFT) -c $(FFLAGS) cross_matrix_.f -o cross_matrix.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o cross_matrix cross_matrix.o

cppcheck:  mpc cppcheck.F
	$(CPP) -P $(CPPFLAGS) cppcheck.F | ./mpc > cppcheck_.f
	$(CFT) -c $(FFLAGS) cppcheck_.f -o cppcheck.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o cppcheck cppcheck.o

srcscheck: mpc srcscheck.F
	$(CPP) -P $(CPPFLAGS) srcscheck.F | ./mpc > srcscheck_.f
	$(CFT) -c $(FFLAGS) srcscheck_.f -o srcscheck.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o srcscheck srcscheck.o
	rm -f check_srcs.F

checkkwds: mpc checkkwds.F
	$(CPP) -P $(CPPFLAGS) checkkwds.F | ./mpc > checkkwds_.f
	$(CFT) -c $(FFLAGS) checkkwds_.f -o checkkwds.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o checkkwds checkkwds.o
	rm -f setup_kwds.F

checkdefs: check_switches1.F setup_kwds.F

check_switches1.F: cppcheck cppdefs.h
	cat cppdefs.h cppdefs_dev.h > mergcpp.txt
	./cppcheck

check_srcs.F: srcscheck Makefile
	./srcscheck

setup_kwds.F: checkkwds read_inp.F
	./checkkwds

partit: partit.F insert_node.o lenstr.o
	$(CPP) -P $(CPPFLAGS) partit.F > partit._.f
	$(CFT) -c $(FFLAGS) partit._.f -o partit.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o partit partit.o insert_node.o  lenstr.o $(LCDF)

ncjoin: ncjoin.F lenstr.o
	$(CPP) -P $(CPPFLAGS) ncjoin.F | ./mpc  > ncjoin._.f
	$(CFT) -c $(FFLAGS) ncjoin._.f -o ncjoin.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o ncjoin ncjoin.o lenstr.o $(LCDF)


ncrename: ncrename.F lenstr.o
	$(CPP) -P $(CPPFLAGS) ncrename.F | ./mpc > ncrename._.f
	$(CFT) -c $(FFLAGS) ncrename._.f -o ncrename.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o ncrename ncrename.o lenstr.o $(LCDF)

depend: checkdefs cross_matrix
	./cross_matrix  *.F90 *.F

mymodules: $(MOBJS) $(AMROBJS)
#
# Target to create tar file.
# ====== == ====== === =====
#
tarfile:
	tar cvf croco.tar Make* *.h *.F *.F90 *.in *.in.*
#
# Cleaning:
# =========
#
rmtools:
	/bin/rm -f $(TOOLS)
clean:
	/bin/rm -rf core *.o $(AMRDIR)/*.o *.i *.s *.f *.f90 *.trace *.mod ${COMP_FILES}

clobber: clean
	/bin/rm -rf $(SBIN) $(TOOLS) ./rii_files


plotter: plotter.F
	f77 -n32 -o plotter plotter.F $(LIBNCAR)

#
# Special treatment for barrier function:
# THERE SHALL BE NO OPTIMIZATION HERE!!!!
#
my_barrier.o: my_barrier.f
	$(CFT) -c -O0 my_barrier.f
#
# Include automatically generated dependency list:
#

-include Make.depend
