#include "cppdefs.h"

      MODULE MOD_TIDES_MAS
#ifdef TIDES_MAS
#include "param.h"
#include "set_global_definitions.h"
      INTEGER,PARAMETER :: nmarmax=144
      INTEGER,DIMENSION(nmarmax)                  :: nu
      INTEGER,DIMENSION(30,8)                     :: nomd
      REAL,DIMENSION(30,8)              :: fr1,fr2
      REAL,DIMENSION(30,8,GLOBAL_1D_ARRAYETA,2)    :: rouest,gouest,rest,gest
      REAL,DIMENSION(30,8,GLOBAL_1D_ARRAYXI,2)    :: rsud,gsud,rnord,gnord
      REAL,DIMENSION(130,GLOBAL_1D_ARRAYXI,2)     :: hnnord,hnsud
      REAL,DIMENSION(130,GLOBAL_1D_ARRAYETA,2)     :: hnouest,hnest
      REAL,DIMENSION(30,8):: rr,gg
      REAL,DIMENSION(130) :: hn

      CHARACTER(len=100),DIMENSION(nmarmax) :: nommar

      LOGICAL,DIMENSION(nmarmax)            :: l_presence

      INTEGER,DIMENSION(nmarmax)            :: doodson

      REAL,DIMENSION(nmarmax)     :: omega

      REAL,DIMENSION(nmarmax)     :: frequence,equitide

      INTEGER :: inivm,isa,issa,imsm,imm,imsf,imf,imfm,i2q1,isigma1
      INTEGER :: iq1,irho1,io1
      INTEGER :: ims1,imp1,ia19,im1,ikhi1,ipi1,i2ns2,i3m2s2,ioq2mnk2,i33,imns2
      INTEGER :: ip1,is1,ik1,ipsi1,iphi1,itheta1,ij1,iso1,ioo1,ikq1,i2mn2s2
      INTEGER :: imnus2,i2mk2,i2n22nm2,imu2ms2,i39,in2,inu2,iop2msk2
      INTEGER :: i43,imsk2,im2,imks2,imks22,ilambda2,i2mk3,im3,iso3
      INTEGER :: ims3,imk3,ia87
      INTEGER :: il22mn2,inkm2,it2,is2,ir22,ik2,imsn2,ikj2mkn2,i2sm2
      INTEGER :: iskm2,imq3
      INTEGER :: isp3,is3,isk3,ik3,i2nms4,i2mmus4,i2mns4,i2mnus4,i3mk4,in4
      INTEGER :: i3ms4,imn4,imnu4,i2msk4,im4,i2msn4,i2mkn4,is4,isk4,i3mnk6
      INTEGER :: i2mks4,isn4,i3mn4ml4,ink4,i2smk4,imt4,ims4,imk4,i2snm4
      INTEGER :: i3mns6,i3mnus6,i4mk6,i2nm6,i4ms6,i2mn6,i2mnu6,i3msk6,im6
      INTEGER :: i3mks6,imsn6,i4mn2ml6,imnk6,i2mt6,i3mnu8,i2msk8,im8,i4mks8
      INTEGER :: i2ms6,i2mk6,i2sn6,i3msn6,i3mkn6,i2sm6,imsk6,i2mns8,i2mn8,i3mn8
      INTEGER :: i2msn8,i3ml8,i2mnk8,i3ms8,i3mk8,i2smn8,i4msn8
      INTEGER :: imsnk8,i4mnk8,i2ms8
      INTEGER :: i4mn10,im10,i3msn10,i4ms10,i4mk10,i5msn10,i2msnk10,i3m2s10
      INTEGER :: i2msk82,i5mns10,i3m2n10
      
      INTEGER,PARAMETER               :: tide_nbharmcp=Ntides
   
!tidepredic
      integer                   :: nol,nd0,ir0,ir1,ir2,nd,jourmem=-1
      integer,save :: nbouc_restart=0
      integer,dimension(8),save :: num
      integer,dimension(144)    :: rr0,rr1,rr2
      integer,dimension(144),PARAMETER :: ide=(/
     &   55555,  56555,  57555,  63655,  65455,  73555,  75555, 125755, 
     &  127555, 135655, 137455, 145555, 146555, 147555, 153655, 155655, 
     &  157455, 162556, 163555, 164555, 165555, 166554, 167555, 173655,
     &  175455, 183555, 185555, 195455, 209655, 217755, 219555, 225655,
     &  225855, 227655, 229455, 235555, 235755, 237555, 238554, 245655,
     &  247455, 253555, 253755, 254556, 255555, 256554, 257555, 263655,
     &  265455, 265655, 272556, 273555, 274554, 275555, 283455, 285455,
     &  291555, 293555, 335655, 345555, 355555, 363555, 364555, 365555,
     &  375555, 381555, 382555, 383555, 385555, 417755, 419555, 427655,
     &  429455, 435555, 435755, 437555, 445655, 447455, 453555, 455555,
     &  457555, 463655, 465455, 465655, 471555, 472556, 473555, 475555,
     &  481655, 483455, 485455, 491555, 493555, 625655, 627655, 629455,
     &  635555, 635755, 637555, 645655, 647455, 653555, 655555, 657555,
     &  663655, 665455, 665655, 672556, 673555, 675555, 681655, 683455,
     &  685455, 691555, 693555, 827655, 835755, 845655, 847455, 853555,
     &  855555, 857555, 863655, 865455, 865655, 873555, 875555, 881655,
     &  883455, 883655, 885455, 891555, 893555,1027655,1035755,1045655,
     &  1055555,1063655,1073555,1075555,1083455,1083655,1091555,  85455/)

       real                                 :: tempis,t_0
       real                                 :: tide_coef

      data rr1/ 
     & 0,0,0,-72,-66,72,0,189,189,189,188,189,0,0,227,0,0,0,0,0,  20,3*0
     & ,0,0,0,0,-112,-75,-112,279,0,-75,-75,279,-37,-37,0,-37,-37,260,0,0
     & ,-37,0,-37,-37,-37,-37,0,2,0,0,-37,-37,0,0,152, 61,-56,189,0,0,0,0
     & ,0, 20, 20,-112,-112,-112,-112,186,-75,-112,-75,-75,223,-75,-75,-37,
     & -112,-37,261,-37,-37,-37,-37,-75,-75,0,0,
     & 148,-149,-149,148,-112,-149,-112,-112,186,-112,-112,-75,-149,-75,
     & -75,-75,-75,-37,-112,-112,-37,-37,
     & -112,-149,-149,-149,148,-149,-149,-112,-149,-112,-112,-112,-75,
     & -149,-75,-149,-75,-75,
     & -224,-186,-186,-186,-149,-149,-149,-186,-112,-112,41/
!
      data rr2/ 
     &  3*0,-65,-65,-64,414, 
     &  6*0,218,0,-200,-219,3*0,-135,3*0,-198,-195,-640,-640,
     &  17*0,-477,298,2*0,260,3*0,298,-37,260,-37,260,
     &  0,20,3 *0,135,439,2*0,-135,-433,
     &  11*0,298,0,-37,298,3*0,298,-37,-37,261,0,298,
     &  10*0,297,0,-37,297,2*0,297,0,-37,260,0,297,
     &  6*0,297,2*0,297,0,297,0,-37,297,260,0,297,
     &  6*0,297,-37,297,0,0/
!
      data rr0/ 
     &  13*1000,4*-1000,3*1000,8*-1000,19*1000,2*-1000,14*1000,
     &  2*-1000,2*1000,2*-1000,54*1000,-1000,19*1000,1000/
!

#endif
      END MODULE MOD_TIDES_MAS
