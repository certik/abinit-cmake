!{\src2tex{textfont=tt}}
!!****f* ABINIT/wrtout
!! NAME
!! wrtout
!!
!! FUNCTION
!! Organizes the sequential or parallel version of the write intrinsic
!! Also allows to treat correctly the write operations for
!! Unix (+DOS) and MacOS.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  unit=unit number for writing
!!  mode_paral='COLL' if all procs are calling the routine with the same message
!!           to be written once only
!!          or 'PERS' if the procs are calling the routine with different mesgs
!!           each to be written, or if one proc is calling the routine
!!
!! OUTPUT
!!  (only writing)
!!
!! SIDE EFFECTS
!!  message=(character(len=500)) message to be written
!!
!! PARENTS
!!      abi_etsf_electrons_put,abi_etsf_geo_put,abi_etsf_init,abinit,acfd_dyson
!!      acfd_intexact,afterscfloop,anaddb,anascr,append_cml,append_cml2,asria9
!!      asrif9,asrprs,assemblychi0sf,assemblychi0sfq0,atm2fft
!!      bands_classification,berryphase,berryphase_new,besjm,bfactor,bigbx9
!!      bldgrp,blok8,bonds_lgth_angles,bound,bound_new,brdmin,calc_cs
!!      calc_density,calc_efg,calc_fc,calc_rpa_functional,calc_sig_ppm
!!      calc_vHxc_braket,calc_wf_qp,calcdensph,canat9,ccfft,ccgradvnl
!!      ccgradvnl_ylm,cchi0,cchi0q0,cgwf,cgwf3,chkdilatmx,chkdpr,chkexi,chkgrp
!!      chki8,chkilwf,chkin9,chkinp,chkint,chkneu,chknm8,chkorthsy,chkpawovlp
!!      chkph3,chkprimit,chkr8,chkrp9,chkvars,chneu9,cigfft,clnup1,clnup2
!!      clsopn,cmpar8,completeperts,constrf,contract_dp_ge_val
!!      contract_int_ge_val,contract_int_le_val,contract_int_list,cppm1par
!!      cppm2par,cppm3par,cppm4par,cprj_utils,crho,crystal_methods,csigme
!!      ctocprj,cutoff_cylinder,cutoff_sphere,cutoff_surface,cvxclda,d3output
!!      ddkten,debug_tools,delocint,der_int,diel9,dielmt,dielmt2,dieltcel
!!      diisrelax,distrb2,dos_hdr_write,dotprod_vn,dotprodm_v,dotprodm_vn
!!      driver,drivexc,dsksta,dtsetcopy,dump_wfgw,dyson_sc,eig1fixed,elast9
!!      electrooptic,elphon,elpolariz,eltxccore,energy,entropyrec,eps1_tc,etot3
!!      ewald,ewald3,ewald4,ewald9,extraprho,fconv,fermi,fftpac,fftw,fftwfn
!!      filterpot,findggp,findk,findmin,findnq,findq,findqg0,fixsym,forstr
!!      forstrnps,fourdp,fourwf,free_energyrec,fsumrule,ftgam,ftgkk,ftiaf9
!!      fxphas,gath3,gensymshub,gensymshub4,gensymspgr,get_all_gkq
!!      get_bands_sym_GW,get_fs_kpts,get_full_gsphere,get_full_kgrid,get_g_tiny
!!      get_gkk_qpt_tr,get_irredg,get_ng0sh,get_tetra,getattribute,getcprj
!!      getcut,getdim_nloc,getfreqsus,getgh1c,getghc,getgsc,getkgrid,getlambda
!!      getmpw,getnel,getng,getshell,getspinrot,green_kernel,gstate,gtblk9
!!      gtdyn9,gw_tools,gwcompleteness,handle_err_netcdf,hartre,hartre1
!!      hartrestr,hdr_check,hdr_init,hdr_io,hdr_io_etsf,hdr_io_netcdf
!!      hdr_vs_dtset,herald,hermit,hessinit,identk,identq,importcml,inarray
!!      incomprs,ingeo,ingeobld,ini_wf_etsf,init8,initang,initberry,initberry3
!!      initmpi_respfn,initmv,initrhoij,initro,initwf,inkpts,inprep8,instr9
!!      instrng,insy3,int2char,int2char4,intagm,integrate_gamma
!!      integrate_gamma_tr,inupper,invars0,invars1,invars1m,invars2,invars2m
!!      invars9,invcb,inwffil,inwffil3,ioarr,ioddb8,iofn1,iofn2,ioppm,iosig
!!      irrzg,isfile,jellium,joint_dos,klocal,kpgio,kpgsph,kpgstr,kramerskronig
!!      kxc_alda,kxc_eok,ladielmt,lattice,lavnl,leave_new,leave_test
!!      lifetime_bn,lifetime_rpa,linemin,listkk,lobpcgIIwf,lobpcgccIIwf
!!      lobpcgccwf,lobpcgwf,loop3dte,loper3,lwf,m_commutator_vkbr,m_errors
!!      m_gwdefs,m_special_funcs,make_epsm1_driver,matcginv,mati3inv,matrginv
!!      matrixelmt_g,mean_fftr,meanvalue_g,memana,memerr,memkss,memorf,memory
!!      merge_kgirr,metcon,metric,metstr,mka2f,mka2fQgrid,mka2f_tr,mkcor3
!!      mkcore,mkdenpos,mkeuler,mkffnl,mkfilename,mkfskgrid,mkifc9,mkkpg,mklocl
!!      mklocl_realspace,mklocl_recipspace,mklocl_wavelets,mknesting,mknormpath
!!      mkph_linwid,mkphdos,mkqptequiv,mkrho,mkrho3,mkvxc3,mkvxcstr3,mlwfovlp
!!      mlwfovlp_proj,mlwfovlp_pw,mlwfovlp_radial,mlwfovlp_setup,moddiel,moldyn
!!      move,mrgddb,mrggkk,mrgscr,mv_3dte,nanal9,nderiv_gen,newfermie1,newkpt
!!      newocc,newrho,newsp,newvtr,newvtr3,nhatgrid,nmsq_gam,nmsq_gam_sumfs
!!      nmsq_pure_gkk,nmsq_pure_gkk_sumfs,nonlinear,nonlop,nonlop_pl,nonlop_ylm
!!      normsq_gkq,normv,nres2vres,nselt3,nstdy3,nstwf3,occeig,occred,odamix
!!      old_setmesh,operat,opernl2,opernl3,opernl4b,optics_paw,out1dm,outelph
!!      outeps,outkss,outqmc,outscfcv,outvars,outwant,outwf,overlap,overlap_g
!!      pareigocc,parsefile,partial_dos_fractions_paw,paw_mkrhox,paw_mkrhox_spl
!!      paw_onsite_vxc,paw_symcprj,paw_tools,pawaccrhoij,pawdenpot,pawdij
!!      pawfgrtab_utils,pawgrnl,pawgylm,pawinit,pawlsylm,pawmknhat,pawmknhat3
!!      pawmkrhoij,pawmkrhoij3,pawnabla_init,pawprt,pawpupot,pawpuxinit,pawr
!!      pawsushat,pawuenergy,pawxc,pawxcm,pawxcsph,pawxenergy,pawxpot,pclock
!!      ph1d3d,phfrq3,piezo9,pl_deriv,plm_coeff,plm_d2theta,plm_dphi,plm_dtheta
!!      polcart,ppmodel_methods,prcref,prcref_PMA,prctfvw1,prctfvw2,precon
!!      precon2,prep_fourwf,prep_getghc,prep_nonlop,print_ij,print_ngfft
!!      print_paw_ij,print_psps,printbxsf,printxsf,projbd,prt_cml,prt_cml2
!!      prteigrs,prtene,prtene3,prtocc,prtph3,prtrhomxmn,prtspgroup,prttagm
!!      prtxf,prtxvf,psddb8,psolver_hartree,psolver_kernel,psolver_rhohxc
!!      psp10in,psp10nl,psp1cc,psp1in,psp1nl,psp2in,psp2lo,psp3in,psp3nl,psp4cc
!!      psp5in,psp5nl,psp6in,psp7in,psp7nl,psp8cc,psp8in,psp9in,pspatm,pspini
!!      pstate,psxml2ab,q0dy3,ramansus,randac,rchkgsheader,rdddb9,rdkss,rdlda
!!      rdm,rdnpw,rdqps,rdscr,read_gkk,readeig,recursion,relaxpol,respfn
!!      respfunc_methods,rhofermi3,rhohxc,rhotov3,rotmat,rrho,rsiaf9,rwwan,rwwf
!!      scalewf_nonlop,scfcge,scfcv,scfcv3,scfeig,scfopt,scprqt,screening
!!      setmesh,setnoccmmp,setshells,setup1,setup2,setup_FFT_rotation
!!      setup_G_rotation,setup_G_rotation_old,setup_coulombian,setup_kmesh
!!      setup_little_group,setup_ppmodel,setup_qmesh,setvtr,sg_ctrig,sg_fft
!!      sg_fftpad,sg_fftpx,sg_fftrisc,sg_fftx,sg_ffty,sg_fftz,sg_fourwf,sigma
!!      sizefft,smallprim,smatrix,smatrix_paw,smatrix_pawinit,smpbz,solve_Dyson
!!      spectral,sphere,sphereboundary,sphericaldens,spin_current,split_work
!!      split_work2,status,stress,subdiago,suscep,suscep_dyn,suscep_kxc_dyn
!!      suscep_stat,suskmm,suskmm_dyn,suskmm_kxc_dyn,sym_gkk,symanal,symatm
!!      symbrav,symdet,symdij,symdm9,symf12,symfind,symg,symkchk,symkpt
!!      symmultsg,symph3,symq3,symrelrot,symrhg,symrhoij,symrhoij3,symspgr
!!      tddft,testkgrid,testlda,testscr,tetrahedron,thm9,timab,timana,time_accu
!!      timein,transgrid,uderiv,vdotw,vlocalstr,vtorho,vtorho3,vtorhorec
!!      vtorhotf,vtowfk,vtowfk3,wannier,wf_info,wfconv,wffclose,wffile,wffopen
!!      wffreadnpwrec,wfkfermi3,wfsinp,wght9,wrqps,wrscr,wrtloctens
!!      wvl_free_type,wvl_init_type_proj,wvl_init_type_wfs,wvl_memory,wvl_mkrho
!!      wvl_newvtr,wvl_nl_gradient,wvl_rwwf,wvl_setboxgeometry,wvl_setngfft
!!      wvl_tail_corrections,wvl_vtorho,wvl_wfsinp_disk,wvl_wfsinp_reformat
!!      wvl_wfsinp_scratch,xc_kernel,xcacfd,xcden,xchcth,xchelu,xcpbe,xcpot
!!      xcpzca,xcspol,xctetr,xcwign,xcxalp,xfpack,xredxcart,ylmc,ylmcd,zprecon3
!!
!! CHILDREN
!!      mpi_comm_rank,wrtout_myproc
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine wrtout(unit,message,mode_paral)

 use defs_basis
#if defined MPI && defined MPI2
 use mpi
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
!End of the abilint section

 implicit none
#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 integer,intent(in) :: unit
 character(len=4),intent(in) :: mode_paral
 character(len=500),intent(inout) :: message

!Local variables-------------------------------
!no_abirules
#if defined MPI 
 integer :: rtnpos
 character(len=7) :: tag
 character(len=500) :: messtmp
 character(len=550) :: string
          !Variables introduced for MPI version
           integer, save :: master=0
           integer :: ierr,me
           character(len=12) :: form,strfmt
#endif

!******************************************************************
!BEGIN EXECUTABLE SECTION

!Be careful with the coding  of the parallel case ...
#if defined MPI 
          !Determine who I am
           call MPI_COMM_RANK(MPI_COMM_WORLD,me,ierr)
           if(mode_paral=='COLL') then
            if(me==master) then
#endif

 call wrtout_myproc(unit, message)

#if defined MPI 
             end if
           elseif(mode_paral=='PERS') then
             if(me<10) then
               write(tag,'("-P-000",i1)') me
             elseif(me<100) then
               write(tag,'("-P-00",i2)') me
             elseif(me<1000) then
               write(tag,'("-P-0",i3)') me
             elseif(me<10000) then
               write(tag,'("-P-",i4)') me
             else
               tag=' ######'
             end if
             rtnpos=index(message,ch10)
             do while(rtnpos/=0)
               write(string, "(A,A,A)") tag, ' ', message(1:rtnpos-1)
               write(unit, '(a)' ) trim(string)
               message=message(rtnpos+1:len(message))
               rtnpos=index(message,ch10)
             end do
             write(string, "(A,A,A)") tag, ' ', message
             write(unit, '(a)' ) trim(string)
           elseif(mode_paral=='INIT') then
             master=unit
           end if
#endif

end subroutine wrtout
!!***
