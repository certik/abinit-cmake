!{\src2tex{textfont=tt}}
!!****f* ABINIT/append_cml
!! NAME
!! append_cml
!!
!! FUNCTION
!! Translate the data from a CML string (string_cml), of length lenstr_cml,
!! and add it at the end of the usual ABINIT input data string (string),
!! taking into account the dtset (dtset_char)
!!
!! COPYRIGHT
!! Copyright (C) 2002-2008 ABINIT group (XG).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtset_char*2=possible dtset label
!!  lenstr_cml=actual number of characters in string
!!  strln=maximal number of characters of string, as declared in the calling routine
!!  string_cml*(strln)=string of characters from the CML file
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  lenstr=actual number of characters in string
!!  string*(strln)=string of characters  (upper case) to which the CML data are appended
!!
!! PARENTS
!!
!! CHILDREN
!!      findmarkup,leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine append_cml (dtset_char,lenstr,lenstr_cml,string,string_cml,strln)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_13xml, except_this_one => append_cml
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lenstr_cml,strln
 integer,intent(inout) :: lenstr
 character(len=2),intent(in) :: dtset_char
 character(len=strln),intent(in) :: string_cml
 character(len=strln),intent(inout) :: string

!Local variables-------------------------------
 character :: blank=' '
!scalars
 integer :: iatom,index_blank,index_lower,index_lower_trial,index_upper,isym
 integer :: lenstr_new,lenstr_old,markuplen,mu,natom,nsym,nu,tspacegroup
 character(len=2) :: string2
 character(len=20) :: string20
 character(len=3) :: string3
 character(len=5) :: string5
 character(len=500) :: message
 character(len=fnlen) :: builtin
!arrays
 integer :: indices_atomArray(3),indices_crystal(3),indices_markup(3)
 integer :: indices_molecule(3),indices_symmetry(3)
 integer,allocatable :: symrel(:,:,:)
 real(dp),allocatable :: tnons(:,:),xred(:,:)

!************************************************************************

!DEBUG
!write(6,*)' append_cml : enter , lenstr=',lenstr
!write(6,*)trim(string(1:lenstr))
!string(lenstr+1:lenstr+5)=' TEST'
!write(6,*)trim(string(1:lenstr+5))
!stop
!ENDDEBUG

 lenstr_new=lenstr

!Find the critical indices of the first 'molecule' mark-up in CML string
 builtin=blank
 index_lower=1
 index_upper=lenstr_cml
 markuplen=8
 call findmarkup(builtin,index_lower,index_upper,indices_molecule,&
& 'molecule',markuplen,strln,string_cml)

!write(6,*)string_cml(indices_molecule(1):indices_molecule(2))

 if(indices_molecule(1)>0)then

  write(message,'(a)') ' Identified CML markup <molecule>'
  call wrtout(6,message,'COLL')
  call wrtout(ab_out,message,'COLL')

! ---------------------------------------------------------------------------

! Find the critical indices of the 'crystal' mark-up,
! inside the first 'molecule' block.
  builtin=blank
  index_lower=indices_molecule(2)
  index_upper=indices_molecule(3)
  markuplen=7
  call findmarkup(builtin,index_lower,index_upper,indices_crystal,&
&  'crystal',markuplen,strln,string_cml)

  if(indices_crystal(1)>0)then

   write(message,'(a)') ' Identified CML markup <crystal>'
   call wrtout(6,message,'COLL')
   call wrtout(ab_out,message,'COLL')

   index_lower=indices_crystal(2)
   index_upper=indices_crystal(3)

!  Find <string builtin="spacegroup"> , and create adequate append
   tspacegroup=0
   builtin='spacegroup'
   markuplen=6
   call findmarkup(builtin,index_lower,index_upper,indices_markup,&
&   'string',markuplen,strln,string_cml)

   if(indices_markup(1)>0)then

    write(message,'(a)') ' Identified CML markup <string builtin="spacegroup">'
    tspacegroup=1
    call wrtout(6,message,'COLL')
    call wrtout(ab_out,message,'COLL')
    lenstr_old=lenstr_new
    lenstr_new=lenstr_new+9+len_trim(dtset_char)+1+(indices_markup(3)-indices_markup(2)-1)
    string(lenstr_old+1:lenstr_new)=&
&    " _SPGROUP"//trim(dtset_char)//blank//string_cml(indices_markup(2)+1:indices_markup(3)-1)

   end if ! Found <string builtin="spacegroup">

!  Find <float builtin="aCell">,<float builtin="bCell">,<float builtin="cCell"> ,
!  and create adequate append
!  WARNING : suppose all are given ! also suppose that angstroms are used !
   builtin='acell'
   markuplen=5
   call findmarkup(builtin,index_lower,index_upper,indices_markup,&
&   'float',markuplen,strln,string_cml)

   if(indices_markup(1)>0)then

    write(message,'(a)') ' Identified CML markup <float builtin="acell">'
    call wrtout(6,message,'COLL')
    call wrtout(ab_out,message,'COLL')

    lenstr_old=lenstr_new
    lenstr_new=lenstr_new+7+len_trim(dtset_char)+1+(indices_markup(3)-indices_markup(2)-1)
    string(lenstr_old+1:lenstr_new)=&
&    " _ACELL"//trim(dtset_char)//blank//string_cml(indices_markup(2)+1:indices_markup(3)-1)

    builtin='bcell'
    call findmarkup(builtin,index_lower,index_upper,indices_markup,&
&    'float',markuplen,strln,string_cml)

    if(indices_markup(1)>0)then

     write(message,'(a)') ' Identified CML markup <float builtin="bcell">'
     call wrtout(6,message,'COLL')
     call wrtout(ab_out,message,'COLL')

     lenstr_old=lenstr_new
     lenstr_new=lenstr_new+1+(indices_markup(3)-indices_markup(2)-1)
     string(lenstr_old+1:lenstr_new)=&
&     blank//string_cml(indices_markup(2)+1:indices_markup(3)-1)

    else

     write(message,'(8a)')ch10,&
&     ' append_cml : ERROR -',ch10,&
&     '  Could not identify <float builtin="bcell">,',ch10,&
&     '  while <float builtin="acell"> is present.',ch10,&
&     '  Action : check your CML file ; it is likely that ABINIT is not yet able to read it.'
     call wrtout(6,message,'COLL')
     call leave_new('COLL')

    end if ! Found <float builtin="beta">

    builtin='ccell'
    call findmarkup(builtin,index_lower,index_upper,indices_markup,&
&    'float',markuplen,strln,string_cml)

    if(indices_markup(1)>0)then

     write(message,'(a)') ' Identified CML markup <float builtin="ccell">'
     call wrtout(6,message,'COLL')
     call wrtout(ab_out,message,'COLL')

     lenstr_old=lenstr_new
     lenstr_new=lenstr_new+1+(indices_markup(3)-indices_markup(2)-1)+9
     string(lenstr_old+1:lenstr_new)=&
&     blank//string_cml(indices_markup(2)+1:indices_markup(3)-1)//" ANGSTROM"

    else

     write(message,'(8a)')ch10,&
&     ' append_cml : ERROR -',ch10,&
&     '  Could not identify <float builtin="ccell">,',ch10,&
&     '  while <float builtin="acell"> is present.',ch10,&
&     '  Action : check your CML file ; it is likely that ABINIT is not yet able to read it.'
     call wrtout(6,message,'COLL')
     call leave_new('COLL')

    end if ! Found <float builtin="ccell">

   end if ! Found <float builtin="acell">

!  Find <float builtin="alpha">,<float builtin="beta">,<float builtin="gamma"> ,
!  and create adequate append
!  WARNING : suppose all are given ! also suppose that degrees are used !
   builtin='alpha'
   markuplen=5
   call findmarkup(builtin,index_lower,index_upper,indices_markup,&
&   'float',markuplen,strln,string_cml)

   if(indices_markup(1)>0)then

    write(message,'(a)') ' Identified CML markup <float builtin="alpha">'
    call wrtout(6,message,'COLL')
    call wrtout(ab_out,message,'COLL')

    lenstr_old=lenstr_new
    lenstr_new=lenstr_new+8+len_trim(dtset_char)+1+(indices_markup(3)-indices_markup(2)-1)
    string(lenstr_old+1:lenstr_new)=&
&    " _ANGDEG"//trim(dtset_char)//blank//string_cml(indices_markup(2)+1:indices_markup(3)-1)

    builtin='beta'
    call findmarkup(builtin,index_lower,index_upper,indices_markup,&
&    'float',markuplen,strln,string_cml)

    if(indices_markup(1)>0)then

     write(message,'(a)') ' Identified CML markup <float builtin="beta">'
     call wrtout(6,message,'COLL')
     call wrtout(ab_out,message,'COLL')

     lenstr_old=lenstr_new
     lenstr_new=lenstr_new+1+(indices_markup(3)-indices_markup(2)-1)
     string(lenstr_old+1:lenstr_new)=&
&     blank//string_cml(indices_markup(2)+1:indices_markup(3)-1)

    else

     write(message,'(8a)')ch10,&
&     ' append_cml : ERROR -',ch10,&
&     '  Could not identify <float builtin="beta">,',ch10,&
&     '  while <float builtin="alpha"> is present.',ch10,&
&     '  Action : check your CML file ; it is likely that ABINIT is not yet able to read it.'
     call wrtout(6,message,'COLL')
     call leave_new('COLL')

    end if ! Found <float builtin="beta">

    builtin='gamma'
    call findmarkup(builtin,index_lower,index_upper,indices_markup,&
&    'float',markuplen,strln,string_cml)

    if(indices_markup(1)>0)then

     write(message,'(a)') ' Identified CML markup <float builtin="gamma">'
     call wrtout(6,message,'COLL')
     call wrtout(ab_out,message,'COLL')

     lenstr_old=lenstr_new
     lenstr_new=lenstr_new+1+(indices_markup(3)-indices_markup(2)-1)
     string(lenstr_old+1:lenstr_new)=&
&     blank//string_cml(indices_markup(2)+1:indices_markup(3)-1)

    else

     write(message,'(8a)')ch10,&
&     ' append_cml : ERROR -',ch10,&
&     '  Could not identify <float builtin="gamma">,',ch10,&
&     '  while <float builtin="alpha"> is present.',ch10,&
&     '  Action : check your CML file ; it is likely that ABINIT is not yet able to read it.'
     call wrtout(6,message,'COLL')
     call leave_new('COLL')

    end if ! Found <float builtin="gamma">

   end if ! Found <float builtin="alpha">

  end if ! Found a crystal markup

! ---------------------------------------------------------------------------

  if(tspacegroup==1)then

   write(message,'(a)') ' Spacegroup being defined, skip <symmetry> section'
   call wrtout(6,message,'COLL')
   call wrtout(ab_out,message,'COLL')

  else ! tspacegroup==0

!  Find the critical indices of the 'symmetry' mark-up,
!  inside the first 'molecule' block.
   builtin=blank
   index_lower=indices_molecule(2)
   index_upper=indices_molecule(3)
   markuplen=8
   call findmarkup(builtin,index_lower,index_upper,indices_symmetry,&
&   'symmetry',markuplen,strln,string_cml)

   if(indices_symmetry(1)>0)then

    write(message,'(a)') ' Identified CML markup <symmetry>'
    call wrtout(6,message,'COLL')
    call wrtout(ab_out,message,'COLL')

    index_lower=indices_symmetry(2)
    index_upper=indices_symmetry(3)

!   Count the number of <floatMatrix builtin="symOp">
    builtin='symOp'
    markuplen=11
    nsym=0

    do
     call findmarkup(builtin,index_lower,index_upper,indices_markup,&
&     'floatMatrix',markuplen,strln,string_cml)
     if(indices_markup(1)>0)then
      nsym=nsym+1
      index_lower=indices_markup(3)
     else
      exit
     end if
    end do

    if(nsym/=0)then

     write(message,'(a,i5,a)') ' Found',nsym,' symmetry operations ; translate them.'
     call wrtout(6,message,'COLL')
     call wrtout(ab_out,message,'COLL')

     write(string3,'(i3)')nsym
     lenstr_old=lenstr_new
     lenstr_new=lenstr_new+6+len_trim(dtset_char)+1+3
     string(lenstr_old+1:lenstr_new)=&
&     " _NSYM"//trim(dtset_char)//blank//string3

     allocate(symrel(3,3,nsym),tnons(3,nsym))

     index_lower=indices_symmetry(2)

!    Read symrel and tnons from CML string
     do isym=1,nsym
      call findmarkup(builtin,index_lower,index_upper,indices_markup,&
&      'floatMatrix',markuplen,strln,string_cml)
      index_lower=indices_markup(3)
      read(string_cml(indices_markup(2)+1:indices_markup(3)-1),*) &
&      symrel(1,1:3,isym),tnons(1,isym),&
&      symrel(2,1:3,isym),tnons(2,isym),&
&      symrel(3,1:3,isym),tnons(3,isym)
     end do

!    Write symrel
     lenstr_old=lenstr_new
     lenstr_new=lenstr_new+8+len_trim(dtset_char)
     string(lenstr_old+1:lenstr_new)=" _SYMREL"//trim(dtset_char)
     do isym=1,nsym
      do mu=1,3
       do nu=1,3
        write(string2,'(i2)')symrel(mu,nu,isym)
        lenstr_old=lenstr_new
        lenstr_new=lenstr_new+3
        string(lenstr_old+1:lenstr_new)=blank//string2
       end do
      end do
     end do

!    Write tnons
     lenstr_old=lenstr_new
     lenstr_new=lenstr_new+7+len_trim(dtset_char)
     string(lenstr_old+1:lenstr_new)=" _TNONS"//trim(dtset_char)
     do isym=1,nsym
      do mu=1,3
       write(string20,'(f20.12)')tnons(mu,isym)
       lenstr_old=lenstr_new
       lenstr_new=lenstr_new+21
       string(lenstr_old+1:lenstr_new)=blank//string20
      end do
     end do

     deallocate(symrel,tnons)

    end if ! nsym/=0

   end if ! Found <symmetry>

  end if ! tspacegroup

! ---------------------------------------------------------------------------

! Find the critical indices of the 'atomArray' mark-up,
! inside the first 'molecule' block.
  builtin=blank
  index_lower=indices_molecule(2)
  index_upper=indices_molecule(3)
  markuplen=9
  call findmarkup(builtin,index_lower,index_upper,indices_atomArray,&
&  'atomArray',markuplen,strln,string_cml)

  if(indices_atomArray(1)>0)then

   write(message,'(a)') ' Identified CML markup <atomArray>'
   call wrtout(6,message,'COLL')
   call wrtout(ab_out,message,'COLL')

   index_lower=indices_atomArray(2)
   index_upper=indices_atomArray(3)

!  Find <stringArray builtin="elementType"> , and create adequate append
   builtin='elementType'
   markuplen=11
   call findmarkup(builtin,index_lower,index_upper,indices_markup,&
&   'stringArray',markuplen,strln,string_cml)

   if(indices_markup(1)>0)then

    write(message,'(a)') ' Identified CML markup <stringArray builtin="elementType">'
    call wrtout(6,message,'COLL')
    call wrtout(ab_out,message,'COLL')

!   Count the number of atoms
    natom=0
    index_lower_trial=indices_markup(2)+1
    if(index_lower_trial<=indices_markup(3)-1)then
     do
      index_blank=index(trim(string_cml(index_lower_trial:indices_markup(3)-1)),blank)
      if(index_blank/=0)then
       if(index_blank/=1)natom=natom+1
       index_lower_trial=index_lower_trial+index_blank
      else
       exit
      end if
     end do
     natom=natom+1
    end if

    write(message,'(a,i5,a)') ' Found',natom,' atoms'
    call wrtout(6,message,'COLL')
    call wrtout(ab_out,message,'COLL')

    write(string5,'(i5)')natom
    lenstr_old=lenstr_new
    lenstr_new=lenstr_new+7+len_trim(dtset_char)+1+5
    string(lenstr_old+1:lenstr_new)=" _NATOM"//trim(dtset_char)//blank//string5

!   Transfer the content of the elementType section
    lenstr_old=lenstr_new
    lenstr_new=lenstr_new+7+len_trim(dtset_char)+1+(indices_markup(3)-indices_markup(2)-1)+4
    string(lenstr_old+1:lenstr_new)=&
&    " _TYPAX"//trim(dtset_char)//blank//string_cml(indices_markup(2)+1:indices_markup(3)-1)//" XX "

!   Find <stringArray builtin="xFract"> , and create adequate append
    builtin='xFract'
    markuplen=11
    call findmarkup(builtin,index_lower,index_upper,indices_markup,&
&    'stringArray',markuplen,strln,string_cml)

    if(indices_markup(1)>0)then

     write(message,'(a)') ' Identified CML markup <stringArray builtin="xFract">'
     call wrtout(6,message,'COLL')
     call wrtout(ab_out,message,'COLL')

     allocate(xred(3,natom))

     read(string_cml(indices_markup(2)+1:indices_markup(3)-1),*) xred(1,1:natom)

     builtin='yFract'
     call findmarkup(builtin,index_lower,index_upper,indices_markup,&
&     'stringArray',markuplen,strln,string_cml)

     if(indices_markup(1)>0)then

      write(message,'(a)') ' Identified CML markup <stringArray builtin="yFract">'
      call wrtout(6,message,'COLL')
      call wrtout(ab_out,message,'COLL')

      read(string_cml(indices_markup(2)+1:indices_markup(3)-1),*) xred(2,1:natom)

     else

      write(message,'(8a)')ch10,&
&      ' append_cml : ERROR -',ch10,&
&      '  Could not identify <stringArray builtin="yFract">,',ch10,&
&      '  while <stringArray builtin="xFract"> is present.',ch10,&
&      '  Action : check your CML file ; it is likely that ABINIT is not yet able to read it.'
      call wrtout(6,message,'COLL')
      call leave_new('COLL')

     end if ! Found <stringArray builtin="yFract">

     builtin='zFract'
     call findmarkup(builtin,index_lower,index_upper,indices_markup,&
&     'stringArray',markuplen,strln,string_cml)

     if(indices_markup(1)>0)then

      write(message,'(a)') ' Identified CML markup <stringArray builtin="zFract">'
      call wrtout(6,message,'COLL')
      call wrtout(ab_out,message,'COLL')

      read(string_cml(indices_markup(2)+1:indices_markup(3)-1),*) xred(3,1:natom)

     else

      write(message,'(8a)')ch10,&
&      ' append_cml : ERROR -',ch10,&
&      '  Could not identify <stringArray builtin="zFract">,',ch10,&
&      '  while <stringArray builtin="xFract"> is present.',ch10,&
&      '  Action : check your CML file ; it is likely that ABINIT is not yet able to read it.'
      call wrtout(6,message,'COLL')
      call leave_new('COLL')

     end if ! Found <stringArray builtin="zFract">

     lenstr_old=lenstr_new
     lenstr_new=lenstr_new+6+len_trim(dtset_char)+1
     string(lenstr_old+1:lenstr_new)=" _XRED"//trim(dtset_char)//blank

     do iatom=1,natom
      do mu=1,3
       write(string20,'(f20.12)')xred(mu,iatom)
       lenstr_old=lenstr_new
       lenstr_new=lenstr_new+20
       string(lenstr_old+1:lenstr_new)=string20
      end do
     end do

     deallocate(xred)

    end if ! Found <stringArray builtin="xFract">

   end if ! Found <stringArray builtin="elementType">

  end if ! Found <atomArray>

! ---------------------------------------------------------------------------

 end if ! Found <molecule>

!Check the length of the string
 if(lenstr_new>strln)then
  write(message,'(6a)')ch10,&
&  ' append_cml : BUG -',ch10,&
&  '  The maximal size of the input variable string has been exceeded.',ch10,&
&  '  The use of a CML file is more character-consuming than the usual input file. Sorry.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

!Update the length of the string
 lenstr=lenstr_new

!DEBUG
!write(6,*)' append_cml : exit , lenstr=',lenstr
!write(6,*)trim(string(1:lenstr))
!stop
!ENDDEBUG

end subroutine append_cml
!!***
