module defs_berry
 implicit none

 interface
  subroutine initberry(dtefield,dtfil,dtset,gmet,kg,mband,mkmem,mpi_enreg,&  
&  mpw,nkpt,npwarr,nsppol,nsym,occ,pwind,pwind_alloc,pwnsfac,&  
&  rprimd,symrec)
   use defs_basis
   use defs_datatypes
   integer,intent(in) :: mband
   integer,intent(in) :: mkmem
   integer,intent(in) :: mpw
   integer,intent(in) :: nkpt
   integer,intent(in) :: nsppol
   integer,intent(in) :: nsym
   integer,intent(out) :: pwind_alloc
   type(efield_type),intent(out) :: dtefield
   type(datafiles_type),intent(in) :: dtfil
   type(dataset_type),intent(inout) :: dtset
   type(MPI_type),intent(inout) :: mpi_enreg
   real(dp),intent(in) :: gmet(3,3)
   integer,intent(in) :: kg(3,mpw*mkmem)
   integer,intent(in) :: npwarr(nkpt)
   real(dp),intent(in) :: occ(mband*nkpt*nsppol)
   integer,pointer :: pwind(:,:,:)
   real(dp),pointer :: pwnsfac(:,:)
   real(dp),intent(in) :: rprimd(3,3)
   integer,intent(in) :: symrec(3,3,nsym)
  end subroutine initberry
 end interface

end module defs_berry
