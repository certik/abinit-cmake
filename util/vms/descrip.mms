#*****************************************************************************
#                                                                            *
#       Author : J.Jansen                                                    *
#       Version : 2.4                                                        *
#       Date : 14 February 2007                                              *
#       Purpose : Abinit make file for OpenVMS                               *
#                                                                            *
#*****************************************************************************
#

MAINS=aim,anaddb,cut3d,mrgddb,newsp,lwf,conducti

# Archives for Src* directories
AR_PARAL=[.src.21paral_md]lib21paral_md.olb
AR_DRIVE=[.src.21drive]lib21drive.olb
AR_CUT3D=[.src.19cut3d]lib19cut3d.olb
AR_SEQ=[.src.18seqpar]lib18abinis.olb
AR_PAR=[.src.18seqpar]lib18abinip.olb
AR_DDB=[.src.17ddb]lib17ddb.olb
AR_LWF=[.src.17lwf]lib17lwf.olb
AR_SUSCEP=[.src.17suscep]lib17suscep.olb
AR_RESPONSE=[.src.16response]lib16response.olb
AR_COMMON=[.src.15common]lib15common.olb
AR_IOWFDENPOT=[.src.14iowfdenpot]lib14iowfdenpot.olb
AR_WFS=[.src.14wfs]lib14wfs.olb
AR_WVL_WFS=[.src.14wvl_wfs]lib14wvl_wfs.olb
AR_RECURSION=[.src.15recursion]lib15recursion.olb
AR_RSPRC=[.src.15rsprc]lib15rsprc.olb
AR_GW=[.src.15gw]lib15gw.olb
AR_IOVARS=[.src.13iovars]lib13iovars.olb
AR_IONETCDF=[.src.13ionetcdf]lib13ionetcdf.olb
AR_PAW=[.src.13paw]lib13paw.olb
AR_RECIPSPACE=[.src.13recipspace]lib13recipspace.olb
AR_XC=[.src.13xc]lib13xc.olb
AR_XML=[.src.13xml]lib13xml.olb
AR_BADER=[.src.14bader]lib14bader.olb
AR_NONLOCAL=[.src.13nonlocal]lib13nonlocal.olb
AR_FFTS=[.src.12ffts]lib12ffts.olb
AR_PSP=[.src.13psp]lib13psp.olb
AR_GEOMETRY=[.src.12geometry]lib12geometry.olb
AR_PARSER=[.src.12parser]lib12parser.olb
AR_SPACEPAR=[.src.12spacepar]lib12spacepar.olb
AR_MGMPI=[.src.01manage_mpi]lib01manage_mpi.olb
AR_CONTRACT=[.src.11contract]lib11contract.olb
AR_CG=[.src.lib01cg]liblib01cg.olb
AR_HIDEMPI=[.src.lib01hidempi]liblib01hidempi.olb
AR_UTIL=[.src.11util]lib11util.olb
AR_BASIS=[.src.00basis]lib00basis.olb
AR_BASIP=[.src.00basis]lib00basip.olb
AR_DEFS=[.src.defs]libdefs.olb
AR_GEOMOPTIM=[.src.16geomoptim]lib16geomoptim.olb
AR_OCCEIG=[.src.14occeig]lib14occeig.olb
AR_NLSTRAIN=[.src.12nlstrain]lib12nlstrain.olb
AR_POISSON=[.src.14poisson]lib14poisson.olb
AR_IOMPI=[.src.13io_mpi]lib13io_mpi.olb

AR_NUMERIC=[.lib.numeric]libnumeric.olb
AR_NUMERICF90=[.lib.numericf90]libnumericf90.olb
AR_FFTNEW=[.src.lib01fftnew]liblib01fftnew.olb

allseq : version,crea_descrip_mms.exe,vms_prepare_input.exe,prepare_input
	$(MMS) numeric,abinis,$(MAINS)

# Create a one-line version definition
version : [.src.defs]defs_infos.F90
	@ write sys$output "version file created"

[.src.defs]defs_infos.F90 : [.src.defs]defs_infos.F90_vms
	copy [.src.defs]defs_infos.F90_vms [.src.defs]defs_infos.F90

crea_descrip_mms.exe : crea_descrip_mms.obj
	link crea_descrip_mms

crea_descrip_mms.obj : [.vms]crea_descrip_mms.f90
	cc/define=(VMS)/comment=as_is/prep=[.vms]crea_descrip_mms.f_\
	[.vms]crea_descrip_mms.f90
	f90/source=free/convert=big_endian/name=lowercase\
	/opti=(tune=host,nopipe)/arch=host/debug\
	/check=(noarg,bo,form,fp,noout,over,powe,under)\
	[.vms]crea_descrip_mms.f_
	delete [.vms]crea_descrip_mms.f_;*

vms_prepare_input.exe : vms_prepare_input.obj
	link vms_prepare_input

vms_prepare_input.obj : [.vms]vms_prepare_input.f90
	cc/define=(VMS)/comment=as_is/prep=[.vms]vms_prepare_input.f_\
	[.vms]vms_prepare_input.f90
	f90/source=free/convert=big_endian/name=lowercase\
	/opti=(tune=host,nopipe)/arch=host/debug\
	/check=(noarg,bo,form,fp,noout,over,powe,under)\
	[.vms]vms_prepare_input.f_
	delete [.vms]vms_prepare_input.f_;*

numeric :
	$(MMS) [.lib.numeric]descrip.mms
	set default [.lib.numeric]
	$(MMS)
	set default [--]

[.src.lib01fftnew]liblib01fftnew.olb :
	$(MMS) [.src.lib01fftnew]descrip.mms
	set default [.src.lib01fftnew]
	$(MMS)
	set default [--]

[.lib.numericf90]libnumericf90.olb :
	$(MMS) [.lib.numericf90]descrip.mms
	set default [.lib.numericf90]
	$(MMS)
	set default [--]

DEP_ABINIS=$(AR_DRIVE) $(AR_PARAL) $(AR_SEQ) $(AR_SUSCEP) \
 $(AR_RESPONSE) $(AR_COMMON) $(AR_IOWFDENPOT) $(AR_WFS) $(AR_WVL_WFS) \
 $(AR_GW) $(AR_RECURSION) $(AR_RSPRC)  \
 $(AR_IOVARS) $(AR_IOMPI) $(AR_IONETCDF) $(AR_PAW) $(AR_RECIPSPACE) $(AR_XC) \
 $(AR_XML) $(AR_NONLOCAL) $(AR_FFTS) $(AR_PSP) $(AR_GEOMETRY) $(AR_PARSER) \
 $(AR_SPACEPAR) \
 $(AR_UTIL) $(AR_CONTRACT) $(AR_MGMPI) $(AR_BASIS) $(AR_DEFS) $(AR_FFTNEW)\
 $(AR_NUMERIC) $(AR_NUMERICF90) $(AR_GEOMOPTIM) $(AR_OCCEIG) $(AR_NLSTRAIN)\
 $(AR_POISSON) $(AR_CG) $(AR_HIDEMPI)

abinis : abinis.exe
	write sys$output " abinis.exe has been made "
	write sys$output " "

abinis.exe : [.src.main]abinis.obj $(DEP_ABINIS)
	write sys$output " abinis.exe will be made "
	write sys$output " "
	open/write optf abinit.opt
	write optf "$(AR_DRIVE)/lib"
	write optf "$(AR_PARAL)/lib"
	write optf "$(AR_SEQ)/lib"
	write optf "$(AR_SUSCEP)/lib"
	write optf "$(AR_RESPONSE)/lib"
	write optf "$(AR_GEOMOPTIM)/lib"
	write optf "$(AR_COMMON)/lib"
	write optf "$(AR_GW)/lib"
	write optf "$(AR_RECURSION)/lib"
	write optf "$(AR_RSPRC)/lib"
	write optf "$(AR_IOWFDENPOT)/lib"
	write optf "$(AR_WFS)/lib"
	write optf "$(AR_WVL_WFS)/lib"
	write optf "$(AR_OCCEIG)/lib"
	write optf "$(AR_IOVARS)/lib"
	write optf "$(AR_IOMPI)/lib"
	write optf "$(AR_IONETCDF)/lib"
	write optf "$(AR_PAW)/lib"
	write optf "$(AR_RECIPSPACE)/lib"
	write optf "$(AR_XC)/lib"
	write optf "$(AR_XML)/lib"
	write optf "$(AR_NONLOCAL)/lib"
	write optf "$(AR_FFTS)/lib"
	write optf "$(AR_PSP)/lib"
	write optf "$(AR_NLSTRAIN)/lib"
	write optf "$(AR_POISSON)/lib"
	write optf "$(AR_GEOMETRY)/lib"
	write optf "$(AR_PARSER)/lib"
	write optf "$(AR_SPACEPAR)/lib"
	write optf "$(AR_UTIL)/lib"
	write optf "$(AR_CONTRACT)/lib"
	write optf "$(AR_MGMPI)/lib"
	write optf "$(AR_FFTNEW)/lib"
	write optf "$(AR_CG)/lib"
	write optf "$(AR_HIDEMPI)/lib"
	write optf "$(AR_BASIS)/lib"
	write optf "$(AR_DEFS)/lib"
	write optf "$(AR_NUMERIC)/lib"
	write optf "$(AR_NUMERICF90)/lib"
	write optf "sys$library:libnetcdf/lib"
	write optf "sys$library:libxmlf90/lib"
	close optf
	link/exec=abinis.exe [.src.main]abinis.obj,[]abinit.opt/opt

[.src.main]abinis.obj : [.src.main]abinit.F90 [.src.defs]libdefs.olb
	set default [.src.main]
	cc/define=(VMS,HAVE_NETCDF,HAVE_XMLF90)/comment=as_is/prep=abinit.f_ abinit.F90
	f90/source=free/convert=big_endian/name=lowercase\
	/opti=(tune=host,nopipe)/arch=host/debug/object=abinis.obj\
	/check=(noarg,bo,form,fp,noout,over,powe,under)\
	/include=([--],[--.src.defs]) abinit.f_
	delete abinit.f_;*
	set default [--]

DEP_AIM= $(AR_BADER) $(AR_GEOMETRY) $(AR_IOWFDENPOT) $(AR_UTIL) $(AR_BASIS)\
	$(AR_DEFS) $(AR_MGMPI) $(AR_NUMERIC) $(AR_PARSER) $(AR_NUMERICF90)\
	$(AR_UTIL) $(AR_IOMPI) $(AR_UTIL)

DEP_AIM_LIB=$(AR_BADER)/lib,$(AR_GEOMETRY)/lib,$(AR_IOWFDENPOT)/lib,\
	$(AR_PARSER)/lib,$(AR_UTIL)/lib,$(AR_MGMPI)/lib,$(AR_BASIS)/lib,\
	$(AR_DEFS)/lib,$(AR_NUMERIC)/lib,$(AR_NUMERICF90)/lib,$(AR_UTIL)/lib,\
	$(AR_IOMPI)/lib,$(AR_UTIL)/lib

aim : aim.exe
	write sys$output " aim.exe has been made "
	write sys$output " "

aim.exe : [.src.main]aim.obj $(DEP_AIM)
	write sys$output " aim.exe will be made "
	write sys$output " "
	link/exec=aim.exe [.src.main]aim.obj,$(DEP_AIM_LIB)

[.src.main]aim.obj : [.src.main]aim.F90 [.src.defs]libdefs.olb
	set default [.src.main]
	cc/define=(VMS,HAVE_NETCDF,HAVE_XMLF90)/comment=as_is/prep=aim.f_ aim.F90
	f90/source=free/convert=big_endian/name=lowercase\
	/opti=(tune=host,nopipe)/arch=host/debug/object=aim.obj\
	/check=(noarg,bo,form,fp,noout,over,powe,under)\
	/include=([--],[--.src.defs]) aim.f_
	delete aim.f_;*
	set default [--]

DEP_ANADDB=$(AR_DDB) $(AR_RESPONSE) $(AR_WFS) $(AR_COMMON) $(AR_IOWFDENPOT)\
 $(AR_RECIPSPACE) $(AR_GEOMETRY) $(AR_NONLOCAL) $(AR_PARSER) $(AR_SPACEPAR) \
 $(AR_UTIL) $(AR_MGMPI) $(AR_BASIS) $(AR_DEFS) $(AR_NUMERIC) $(AR_NUMERICF90) \
 $(AR_UTIL) $(AR_OCCEIG) $(AR_HIDEMPI) $(AR_IOMPI) $(AR_WVL_WFS)

DEP_ANADDB_LIB=$(AR_DDB)/lib,$(AR_RESPONSE)/lib,$(AR_WFS)/lib,\
	$(AR_WVL_WFS)/lib,\
	$(AR_COMMON)/lib,$(AR_IOWFDENPOT)/lib,$(AR_OCCEIG)/lib,\
	$(AR_RECIPSPACE)/lib,$(AR_IOMPI)/lib,$(AR_HIDEMPI)/lib,\
	$(AR_GEOMETRY)/lib,$(AR_NONLOCAL)/lib,$(AR_PARSER)/lib,\
	$(AR_SPACEPAR)/lib,$(AR_UTIL)/lib,$(AR_MGMPI)/lib,$(AR_BASIS)/lib,\
	$(AR_DEFS)/lib,$(AR_NUMERIC)/lib,$(AR_NUMERICF90)/lib,$(AR_UTIL)/lib

anaddb : anaddb.exe
	write sys$output " anaddb.exe has been made "
	write sys$output " "

anaddb.exe : [.src.main]anaddb.obj $(DEP_ANADDB)
	write sys$output " anaddb.exe will be made "
	write sys$output " "
	link/exec=anaddb.exe [.src.main]anaddb.obj,$(DEP_ANADDB_LIB),\
	sys$library:libnetcdf/lib

[.src.main]anaddb.obj : [.src.main]anaddb.F90 [.src.defs]libdefs.olb
	set default [.src.main]
	cc/define=(VMS,HAVE_NETCDF,HAVE_XMLF90)/comment=as_is/prep=anaddb.f_ anaddb.F90
	f90/source=free/convert=big_endian/name=lowercase\
	/opti=(tune=host,nopipe)/arch=host/debug/object=anaddb.obj\
	/check=(noarg,bo,form,fp,noout,over,powe,under)\
	/include=([--],[--.src.defs]) anaddb.f_
	delete anaddb.f_;*
	set default [--]

DEP_CUT3D=$(AR_CUT3D) $(AR_SEQ) $(AR_COMMON) $(AR_WFS) $(AR_IOWFDENPOT) \
	$(AR_GEOMETRY) $(AR_FFTS) $(AR_SPACEPAR) $(AR_RECIPSPACE) $(AR_UTIL) \
	$(AR_CONTRACT) $(AR_MGMPI) $(AR_BASIS) $(AR_DEFS) $(AR_FFTNEW)\
	$(AR_NUMERIC) $(AR_NONLOCAL) $(AR_NUMERICF90) $(AR_UTIL)\
	$(AR_OCCEIG) $(AR_PARSER) $(AR_HIDEMPI) $(AR_IOMPI) $(AR_WVL_WFS)

DEP_CUT3D_LIB=$(AR_CUT3D)/lib,$(AR_SEQ)/lib,$(AR_COMMON)/lib,$(AR_WFS)/lib,\
	$(AR_WVL_WFS)/lib,\
	$(AR_IOWFDENPOT)/lib,$(AR_OCCEIG)/lib,$(AR_GEOMETRY)/lib,\
	$(AR_FFTS)/lib,\
	$(AR_SPACEPAR)/lib,$(AR_RECIPSPACE)/lib,$(AR_UTIL)/lib,\
	$(AR_PARSER)/lib,$(AR_IOMPI)/lib,$(AR_HIDEMPI)/lib,\
	$(AR_NONLOCAL)/lib,$(AR_CONTRACT)/lib,$(AR_MGMPI)/lib,$(AR_BASIS)/lib,\
	$(AR_DEFS)/lib,$(AR_FFTNEW)/lib,$(AR_NUMERIC)/lib,\
	$(AR_NUMERICF90)/lib,$(AR_UTIL)/lib,sys$library:libnetcdf/lib

cut3d : cut3d.exe
	write sys$output " cut3d.exe has been made "
	write sys$output " "

cut3d.exe : [.src.main]cut3d.obj $(DEP_CUT3D)
	write sys$output " cut3d.exe will be made "
	write sys$output " "
	link/exec=cut3d.exe [.src.main]cut3d.obj,$(DEP_CUT3D_LIB)

[.src.main]cut3d.obj : [.src.main]cut3d.F90 [.src.defs]libdefs.olb
	set default [.src.main]
	cc/define=(VMS,HAVE_NETCDF,HAVE_XMLF90)/comment=as_is/prep=cut3d.f_ cut3d.F90
	f90/source=free/convert=big_endian/name=lowercase\
	/opti=(tune=host,nopipe)/arch=host/debug/object=cut3d.obj\
	/check=(noarg,bo,form,fp,noout,over,powe,under)\
	/include=([--],[--.src.defs]) cut3d.f_
	delete cut3d.f_;*
	set default [--]

DEP_MRGDDB=$(AR_DDB) $(AR_RESPONSE) $(AR_UTIL) $(AR_MGMPI) $(AR_BASIS)\
	$(AR_DEFS) $(AR_NUMERIC)

DEP_MRGDDB_LIB=$(AR_DDB)/lib,$(AR_RESPONSE)/lib,$(AR_UTIL)/lib,\
	$(AR_MGMPI)/lib,$(AR_BASIS)/lib,$(AR_DEFS)/lib,$(AR_NUMERIC)/lib

mrgddb : mrgddb.exe
	write sys$output " mrgddb.exe has been made "
	write sys$output " "

mrgddb.exe : [.src.main]mrgddb.obj $(DEP_MRGDDB)
	write sys$output " mrgddb.exe will be made "
	write sys$output " "
	link/exec=mrgddb.exe [.src.main]mrgddb.obj,$(DEP_MRGDDB_LIB)

[.src.main]mrgddb.obj : [.src.main]mrgddb.F90 [.src.defs]libdefs.olb
	set default [.src.main]
	cc/define=(VMS,HAVE_NETCDF,HAVE_XMLF90)/comment=as_is/prep=mrgddb.f_ mrgddb.F90
	f90/source=free/convert=big_endian/name=lowercase\
	/opti=(tune=host,nopipe)/arch=host/debug/object=mrgddb.obj\
	/check=(noarg,bo,form,fp,noout,over,powe,under)\
	/include=([--],[--.src.defs]) mrgddb.f_
	delete mrgddb.f_;*
	set default [--]

DEP_NEWSP=$(AR_SEQ) $(AR_COMMON) $(AR_WFS) $(AR_IOWFDENPOT) $(AR_IOVARS) \
	$(AR_IONETCDF) $(AR_HIDEMPI) $(AR_IOMPI) $(AR_WVL_WFS) \
	$(AR_RECIPSPACE) $(AR_XML) $(AR_NONLOCAL)     \
	$(AR_FFTS) $(AR_GEOMETRY) $(AR_PARSER) $(AR_SPACEPAR)         \
	$(AR_UTIL) $(AR_CONTRACT) $(AR_MGMPI) $(AR_BASIS) $(AR_DEFS)\
	$(AR_FFTNEW) $(AR_NUMERIC) $(AR_NUMERICF90) $(AR_UTIL)

DEP_NEWSP_LIB=$(AR_SEQ)/lib,$(AR_COMMON)/lib,$(AR_WFS)/lib,$(AR_WVL_WFS)/lib,\
	$(AR_IOWFDENPOT)/lib,$(AR_IOVARS)/lib,$(AR_IONETCDF)/lib,$(AR_RECIPSPACE)/lib,\
	$(AR_XML)/lib,$(AR_NONLOCAL)/lib,$(AR_FFTS)/lib,$(AR_GEOMETRY)/lib,\
	$(AR_PARSER)/lib,$(AR_SPACEPAR)/lib,$(AR_UTIL)/lib,$(AR_CONTRACT)/lib,\
	$(AR_IOMPI)/lib,$(AR_HIDEMPI)/lib,\
	$(AR_MGMPI)/lib,$(AR_BASIS)/lib,$(AR_DEFS)/lib,$(AR_FFTNEW)/lib,\
	$(AR_NUMERIC)/lib,$(AR_NUMERICF90)/lib,$(AR_UTIL)/lib,\
	sys$library:libnetcdf/lib

newsp : newsp.exe
	write sys$output " newsp.exe has been made "
	write sys$output " "

newsp.exe : [.src.main]newsp.obj $(DEP_NEWSP)
	write sys$output " newsp.exe will be made "
	write sys$output " "
	link/exec=newsp.exe [.src.main]newsp.obj,$(DEP_NEWSP_LIB)

[.src.main]newsp.obj : [.src.main]newsp.F90 [.src.defs]libdefs.olb
	set default [.src.main]
	cc/define=(VMS,HAVE_NETCDF,HAVE_XMLF90)/comment=as_is/prep=newsp.f_ newsp.F90
	f90/source=free/convert=big_endian/name=lowercase\
	/opti=(tune=host,nopipe)/arch=host/debug/object=newsp.obj\
	/check=(noarg,bo,form,fp,noout,over,powe,under)\
	/include=([--],[--.src.defs]) newsp.f_
	delete newsp.f_;*
	set default [--]

DEP_LWF=$(AR_LWF) $(AR_PARSER) $(AR_UTIL) $(AR_MGMPI) $(AR_BASIS) $(AR_DEFS)\
	$(AR_NUMERIC)

DEP_LWF_LIB=$(AR_LWF)/lib,$(AR_PARSER)/lib,$(AR_UTIL)/lib,$(AR_MGMPI)/lib,\
	$(AR_BASIS)/lib,$(AR_DEFS)/lib,$(AR_NUMERIC)/lib

lwf : lwf.exe
	write sys$output " lwf.exe has been made "
	write sys$output " "

lwf.exe : [.src.main]lwf.obj $(DEP_LWF)
	write sys$output " lwf.exe will be made "
	write sys$output " "
	link/exec=lwf.exe [.src.main]lwf.obj,$(DEP_LWF_LIB)

[.src.main]lwf.obj : [.src.main]lwf.F90 [.src.defs]libdefs.olb
	set default [.src.main]
	cc/define=(VMS,HAVE_NETCDF,HAVE_XMLF90)/comment=as_is/prep=lwf.f_ lwf.F90
	f90/source=free/convert=big_endian/name=lowercase\
	/opti=(tune=host,nopipe)/arch=host/debug/object=lwf.obj\
	/check=(noarg,bo,form,fp,noout,over,powe,under)\
	/include=([--],[--.src.defs]) lwf.f_
	delete lwf.f_;*
	set default [--]

DEP_CONDUCTI=$(AR_COMMON) $(AR_IOWFDENPOT) $(AR_GEOMETRY) $(AR_UTIL)\
	$(AR_MGMPI) $(AR_BASIS) $(AR_DEFS) $(AR_NUMERIC) $(AR_NUMERICF90)\
	$(AR_UTIL) $(AR_OCCEIG) $(AR_IOMPI) $(AR_WVL_WFS)

DEP_CONDUCTI_LIB=$(AR_COMMON)/lib,$(AR_IOWFDENPOT)/lib,$(AR_OCCEIG)/lib,\
	$(AR_GEOMETRY)/lib,$(AR_IOMPI)/lib,$(AR_WVL_WFS)/lib,\
	$(AR_UTIL)/lib,$(AR_MGMPI)/lib,$(AR_BASIS)/lib,$(AR_DEFS)/lib,\
	$(AR_NUMERIC)/lib,$(AR_NUMERICF90)/lib,$(AR_UTIL)/lib,\
	sys$library:libnetcdf/lib

conducti : conducti.exe
	write sys$output " conducti.exe has been made "
	write sys$output " "

conducti.exe : [.src.main]conducti.obj $(DEP_CONDUCTI)
	write sys$output " conducti.exe will be made "
	write sys$output " "
	link/exec=conducti.exe [.src.main]conducti.obj,$(DEP_CONDUCTI_LIB)

[.src.main]conducti.obj : [.src.main]conducti.F90 [.src.defs]libdefs.olb
	set default [.src.main]
	cc/define=(VMS,HAVE_NETCDF,HAVE_XMLF90)/comment=as_is/prep=conducti.f_ conducti.F90
	f90/source=free/convert=big_endian/name=lowercase\
	/opti=(tune=host,nopipe)/arch=host/debug/object=conducti.obj\
	/check=(noarg,bo,form,fp,noout,over,powe,under)\
	/include=([--],[--.src.defs]) conducti.f_
	delete conducti.f_;*
	set default [--]

[.src.defs]libdefs.olb :
	$(MMS) [.src.defs]descrip.mms
	set default [.src.defs]
	$(MMS) m_pseudo_types.obj
	$(MMS) defs_basis.obj
	$(MMS) defs_datatypes.obj
	$(MMS) interfaces_12spacepar.obj
	$(MMS) interfaces_13recipspace.obj
	$(MMS) interfaces_11util.obj
	$(MMS) interfaces_lib01cg.obj
	$(MMS) interfaces_13xc.obj
	$(MMS) interfaces_12ffts.obj
	$(MMS)
	set default [--]

[.src.21drive]lib21drive.olb :
	$(MMS) [.src.21drive]descrip.mms
	set default [.src.21drive]
	$(MMS)
	set default [--]

[.src.21paral_md]lib21paral_md.olb :
	$(MMS) [.src.21paral_md]descrip.mms
	set default [.src.21paral_md]
	$(MMS)
	set default [--]

[.src.16geomoptim]lib16geomoptim.olb :
	$(MMS) [.src.16geomoptim]descrip.mms
	set default [.src.16geomoptim]
	$(MMS)
	set default [--]

[.src.14occeig]lib14occeig.olb :
	$(MMS) [.src.14occeig]descrip.mms
	set default [.src.14occeig]
	$(MMS)
	set default [--]

[.src.12nlstrain]lib12nlstrain.olb :
	$(MMS) [.src.12nlstrain]descrip.mms
	set default [.src.12nlstrain]
	$(MMS)
	set default [--]

[.src.14poisson]lib14poisson.olb :
	$(MMS) [.src.14poisson]descrip.mms
	set default [.src.14poisson]
	$(MMS)
	set default [--]

[.src.19cut3d]lib19cut3d.olb :
	$(MMS) [.src.19cut3d]descrip.mms
	set default [.src.19cut3d]
	$(MMS)
	set default [--]

[.src.18seqpar]lib18abinis.olb :
	$(MMS) [.src.18seqpar]descrip.mms
	set default [.src.18seqpar]
	$(MMS)
	copy lib18seqpar.olb lib18abinis.olb
	set default [--]

[.src.17suscep]lib17suscep.olb :
	$(MMS) [.src.17suscep]descrip.mms
	set default [.src.17suscep]
	$(MMS)
	set default [--]

[.src.17ddb]lib17ddb.olb :
	$(MMS) [.src.17ddb]descrip.mms
	set default [.src.17ddb]
	$(MMS)
	set default [--]

[.src.17lwf]lib17lwf.olb :
	$(MMS) [.src.17lwf]descrip.mms
	set default [.src.17lwf]
	$(MMS)
	set default [--]

[.src.16response]lib16response.olb :
	$(MMS) [.src.16response]descrip.mms
	set default [.src.16response]
	$(MMS)
	set default [--]

[.src.15common]lib15common.olb :
	$(MMS) [.src.15common]descrip.mms
	set default [.src.15common]
	$(MMS)
	set default [--]

[.src.14iowfdenpot]lib14iowfdenpot.olb :
	$(MMS) [.src.14iowfdenpot]descrip.mms
	set default [.src.14iowfdenpot]
	$(MMS)
	set default [--]

[.src.14wfs]lib14wfs.olb :
	$(MMS) [.src.14wfs]descrip.mms
	set default [.src.14wfs]
	$(MMS)
	set default [--]

[.src.14wvl_wfs]lib14wvl_wfs.olb :
	$(MMS) [.src.14wvl_wfs]descrip.mms
	set default [.src.14wvl_wfs]
	$(MMS)
	set default [--]

[.src.15gw]lib15gw.olb :
	$(MMS) [.src.15gw]descrip.mms
	set default [.src.15gw]
	$(MMS)
	set default [--]

[.src.15recursion]lib15recursion.olb :
	$(MMS) [.src.15recursion]descrip.mms
	set default [.src.15recursion]
	$(MMS)
	set default [--]

[.src.15rsprc]lib15rsprc.olb :
	$(MMS) [.src.15rsprc]descrip.mms
	set default [.src.15rsprc]
	$(MMS)
	set default [--]

[.src.lib01cg]liblib01cg.olb :
	$(MMS) [.src.lib01cg]descrip.mms
	set default [.src.lib01cg]
	$(MMS)
	set default [--]

[.src.lib01hidempi]liblib01hidempi.olb :
	$(MMS) [.src.lib01hidempi]descrip.mms
	set default [.src.lib01hidempi]
	$(MMS)
	set default [--]

[.src.13iovars]lib13iovars.olb :
	$(MMS) [.src.13iovars]descrip.mms
	set default [.src.13iovars]
	$(MMS)
	set default [--]

[.src.13io_mpi]lib13io_mpi.olb :
	$(MMS) [.src.13io_mpi]descrip.mms
	set default [.src.13io_mpi]
	$(MMS)
	set default [--]

[.src.13ionetcdf]lib13ionetcdf.olb :
	$(MMS) [.src.13ionetcdf]descrip.mms
	set default [.src.13ionetcdf]
	$(MMS)
	set default [--]

[.src.13paw]lib13paw.olb :
	$(MMS) [.src.13paw]descrip.mms
	set default [.src.13paw]
	$(MMS)
	set default [--]

[.src.13recipspace]lib13recipspace.olb :
	$(MMS) [.src.13recipspace]descrip.mms
	set default [.src.13recipspace]
	$(MMS)
	set default [--]

[.src.13xc]lib13xc.olb :
	$(MMS) [.src.13xc]descrip.mms
	set default [.src.13xc]
	$(MMS)
	set default [--]

[.src.13xml]lib13xml.olb :
	$(MMS) [.src.13xml]descrip.mms
	set default [.src.13xml]
	$(MMS)
	set default [--]

[.src.13nonlocal]lib13nonlocal.olb :
	$(MMS) [.src.13nonlocal]descrip.mms
	set default [.src.13nonlocal]
	$(MMS)
	set default [--]

[.src.12ffts]lib12ffts.olb :
	$(MMS) [.src.12ffts]descrip.mms
	set default [.src.12ffts]
	$(MMS)
	set default [--]

[.src.13psp]lib13psp.olb :
	$(MMS) [.src.13psp]descrip.mms
	set default [.src.13psp]
	$(MMS) smoothvlocal.obj
	$(MMS)
	set default [--]

[.src.12geometry]lib12geometry.olb :
	$(MMS) [.src.12geometry]descrip.mms
	set default [.src.12geometry]
	$(MMS)
	set default [--]

[.src.12parser]lib12parser.olb :
	$(MMS) [.src.12parser]descrip.mms
	set default [.src.12parser]
	$(MMS)
	set default [--]

[.src.12spacepar]lib12spacepar.olb :
	$(MMS) [.src.12spacepar]descrip.mms
	set default [.src.12spacepar]
	$(MMS)
	set default [--]

[.src.14bader]lib14bader.olb :
	$(MMS) [.src.14bader]descrip.mms
	set default [.src.14bader]
	$(MMS)
	set default [--]

[.src.11contract]lib11contract.olb :
	$(MMS) [.src.11contract]descrip.mms
	set default [.src.11contract]
	$(MMS)
	set default [--]

[.src.01manage_mpi]lib01manage_mpi.olb :
	$(MMS) [.src.01manage_mpi]descrip.mms
	set default [.src.01manage_mpi]
	$(MMS)
	set default [--]

[.src.11util]lib11util.olb :
	$(MMS) [.src.11util]descrip.mms
	set default [.src.11util]
	$(MMS)
	set default [--]

[.src.00basis]lib00basis.olb :
	$(MMS) [.src.00basis]descrip.mms
	set default [.src.00basis]
	$(MMS)
	set default [--]

[.src.lib01fftnew]descrip.mms : crea_descrip_mms.exe
	set def [.src.lib01fftnew]
	run [--]crea_descrip_mms
	set def [--]

[.lib.numeric]descrip.mms : crea_descrip_mms.exe
	set def [.lib.numeric]
	run [--]crea_descrip_mms
	set def [--]

[.lib.numericf90]descrip.mms : crea_descrip_mms.exe
	set def [.lib.numericf90]
	run [--]crea_descrip_mms
	set def [--]

[.src.00basis]descrip.mms : crea_descrip_mms.exe
	set def [.src.00basis]
	run [--]crea_descrip_mms
	set def [--]

[.src.01manage_mpi]descrip.mms : crea_descrip_mms.exe
	set def [.src.01manage_mpi]
	run [--]crea_descrip_mms
	set def [--]

[.src.11contract]descrip.mms : crea_descrip_mms.exe
	set def [.src.11contract]
	run [--]crea_descrip_mms
	set def [--]

[.src.11util]descrip.mms : crea_descrip_mms.exe
	set def [.src.11util]
	run [--]crea_descrip_mms
	set def [--]

[.src.14bader]descrip.mms : crea_descrip_mms.exe
	set def [.src.14bader]
	run [--]crea_descrip_mms
	set def [--]

[.src.12ffts]descrip.mms : crea_descrip_mms.exe
	set def [.src.12ffts]
	run [--]crea_descrip_mms
	set def [--]

[.src.12geometry]descrip.mms : crea_descrip_mms.exe
	set def [.src.12geometry]
	run [--]crea_descrip_mms
	set def [--]

[.src.13nonlocal]descrip.mms : crea_descrip_mms.exe
	set def [.src.13nonlocal]
	run [--]crea_descrip_mms
	set def [--]

[.src.12parser]descrip.mms : crea_descrip_mms.exe
	set def [.src.12parser]
	run [--]crea_descrip_mms
	set def [--]

[.src.13psp]descrip.mms : crea_descrip_mms.exe
	set def [.src.13psp]
	run [--]crea_descrip_mms
	set def [--]

[.src.12spacepar]descrip.mms : crea_descrip_mms.exe
	set def [.src.12spacepar]
	run [--]crea_descrip_mms
	set def [--]

[.src.15gw]descrip.mms : crea_descrip_mms.exe
	set def [.src.15gw]
	run [--]crea_descrip_mms
	set def [--]

[.src.15recursion]descrip.mms : crea_descrip_mms.exe
	set def [.src.15recursion]
	run [--]crea_descrip_mms
	set def [--]

[.src.15rsprc]descrip.mms : crea_descrip_mms.exe
	set def [.src.15rsprc]
	run [--]crea_descrip_mms
	set def [--]

[.src.lib01cg]descrip.mms : crea_descrip_mms.exe
	set def [.src.lib01cg]
	run [--]crea_descrip_mms
	set def [--]

[.src.lib01hidempi]descrip.mms : crea_descrip_mms.exe
	set def [.src.lib01hidempi]
	run [--]crea_descrip_mms
	set def [--]

[.src.13iovars]descrip.mms : crea_descrip_mms.exe
	set def [.src.13iovars]
	run [--]crea_descrip_mms
	set def [--]

[.src.13io_mpi]descrip.mms : crea_descrip_mms.exe
	set def [.src.13io_mpi]
	run [--]crea_descrip_mms
	set def [--]

[.src.13ionetcdf]descrip.mms : crea_descrip_mms.exe
	set def [.src.13ionetcdf]
	run [--]crea_descrip_mms
	set def [--]

[.src.13paw]descrip.mms : crea_descrip_mms.exe
	set def [.src.13paw]
	run [--]crea_descrip_mms
	set def [--]

[.src.13recipspace]descrip.mms : crea_descrip_mms.exe
	set def [.src.13recipspace]
	run [--]crea_descrip_mms
	set def [--]

[.src.13xc]descrip.mms : crea_descrip_mms.exe
	set def [.src.13xc]
	run [--]crea_descrip_mms
	set def [--]

[.src.13xml]descrip.mms : crea_descrip_mms.exe
	set def [.src.13xml]
	run [--]crea_descrip_mms
	set def [--]

[.src.14iowfdenpot]descrip.mms : crea_descrip_mms.exe
	set def [.src.14iowfdenpot]
	run [--]crea_descrip_mms
	set def [--]

[.src.14wvl_wfs]descrip.mms : crea_descrip_mms.exe
	set def [.src.14wvl_wfs]
	run [--]crea_descrip_mms
	set def [--]

[.src.15common]descrip.mms : crea_descrip_mms.exe
	set def [.src.15common]
	run [--]crea_descrip_mms
	set def [--]

[.src.16response]descrip.mms : crea_descrip_mms.exe
	set def [.src.16response]
	run [--]crea_descrip_mms
	set def [--]

[.src.17ddb]descrip.mms : crea_descrip_mms.exe
	set def [.src.17ddb]
	run [--]crea_descrip_mms
	set def [--]

[.src.17lwf]descrip.mms : crea_descrip_mms.exe
	set def [.src.17lwf]
	run [--]crea_descrip_mms
	set def [--]

[.src.17suscep]descrip.mms : crea_descrip_mms.exe
	set def [.src.17suscep]
	run [--]crea_descrip_mms
	set def [--]

[.src.18seqpar]descrip.mms : crea_descrip_mms.exe
	set def [.src.18seqpar]
	run [--]crea_descrip_mms
	set def [--]

[.src.19cut3d]descrip.mms : crea_descrip_mms.exe
	set def [.src.19cut3d]
	run [--]crea_descrip_mms
	set def [--]

[.src.21drive]descrip.mms : crea_descrip_mms.exe
	set def [.src.21drive]
	run [--]crea_descrip_mms
	set def [--]

[.src.21paral_md]descrip.mms : crea_descrip_mms.exe
	set def [.src.21paral_md]
	run [--]crea_descrip_mms
	set def [--]

[.src.12nlstrain]descrip.mms : crea_descrip_mms.exe
	set def [.src.12nlstrain]
	run [--]crea_descrip_mms
	set def [--]

[.src.14poisson]descrip.mms : crea_descrip_mms.exe
	set def [.src.14poisson]
	run [--]crea_descrip_mms
	set def [--]

[.src.16geomoptim]descrip.mms : crea_descrip_mms.exe
	set def [.src.16geomoptim]
	run [--]crea_descrip_mms
	set def [--]

[.src.14occeig]descrip.mms : crea_descrip_mms.exe
	set def [.src.14occeig]
	run [--]crea_descrip_mms
	set def [--]

[.src.defs]descrip.mms : crea_descrip_mms.exe
	set def [.src.defs]
	run [--]crea_descrip_mms
	set def [--]

prepare_input :
	set def [.tests.tutorial.Input]
	mcr [---]vms_prepare_input t1x.files
	mcr [---]vms_prepare_input t2x.files
	mcr [---]vms_prepare_input t3x.files
	mcr [---]vms_prepare_input t4x.files
	set def [---]

# Build-in tests
tests : test1 test2 test3 test4 test5 test6
	write sys$output "Performed test 1-6"

test12 : test1 test2 
	write sys$output "Performed test 1-2"

test123 : test1 test2 test3 
	write sys$output "Performed test 1-3"

test124 : test1 test2 test4 
	write sys$output "Performed test 1,2,4"

test125 : test1 test2 test5 
	write sys$output "Performed test 1,2,5"

test1 :
	@ @[.tests.Scripts]run-basic-tests.com [.tests.built-in] 1

test2 :
	@ @[.tests.Scripts]run-basic-tests.com [.tests.built-in] 2

test3 :
	@ @[.tests.Scripts]run-basic-tests.com [.tests.built-in] 3

test4 :
	@ @[.tests.Scripts]run-basic-tests.com [.tests.built-in] 4

test5 :
	@ @[.tests.Scripts]run-basic-tests.com [.tests.built-in] 5

test6 :
	@ @[.tests.Scripts]run-basic-tests.com [.tests.built-in] 6
