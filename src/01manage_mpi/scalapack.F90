!{\src2tex{textfont=tt}}
!!****f* ABINIT/scalapack
!! NAME
!! scalapack
!!
!! FUNCTION
!! This module contains functions and subroutine using ScaLAPACK library.
!! The code have to be compiled with the HAVE_SCALAPACK CPP flags.
!!
!! COPYRIGHT
!! Copyright (C) 2001-2008 ABINIT group (CS,GZ,FB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/Infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#if defined HAVE_SCALAPACK

!-------------------------------------------------------
! set up of a processor grid for ScaLAPACK
! as a function of the total number of processors attributed to the grid
!-------------------------------------------------------

SUBROUTINE build_grid_scalapack(grid,nbprocs, communicator)

  use defs_scalapack
  use defs_basis

  implicit none

  TYPE(grid_scalapack),INTENT(out)     :: grid
  INTEGER,INTENT(in)                     :: nbprocs
  INTEGER, INTENT(in)                    :: communicator

  INTEGER  :: i

  grid%nbprocs=nbprocs

  ! Search for a rectangular grid of processors 
  i=INT(SQRT(float(nbprocs)))
  DO WHILE (MOD(nbprocs,i) /= 0)
     i = i-1
  END DO

  grid%dims(1) = i
  grid%dims(2) = INT(nbprocs/i)

  grid%ictxt = communicator

  CALL BLACS_GRIDINIT(grid%ictxt,'R',grid%dims(1),grid%dims(2))           

END SUBROUTINE build_grid_scalapack
!!***


!-------------------------------------------------------
! Build of the data related to one processor in a grid
!-------------------------------------------------------

SUBROUTINE build_processor_scalapack(processor,grid,myproc, comm)

  use defs_scalapack
  use defs_basis

  implicit none

  TYPE(processor_scalapack),INTENT(out)  :: processor
  TYPE(grid_scalapack),INTENT(in)       :: grid
  INTEGER,INTENT(in)                      :: myproc
  INTEGER,INTENT(in)                      :: comm

  processor%grid= grid

  processor%myproc = myproc

  processor%comm = comm

  CALL BLACS_GRIDINFO(grid%ictxt,processor%grid%dims(1), &
       &              processor%grid%dims(2),processor%coords(1), &
       &              processor%coords(2))

  ! These values are the same as those computed by BLACS_GRIDINFO
  ! except in the case where the mmyproc argument is not the
  ! local proc
  processor%coords(1) = INT((myproc) / grid%dims(2))
  processor%coords(2) = MOD((myproc), grid%dims(2))


END SUBROUTINE build_processor_scalapack

!-------------------------------------------------------
! initialisation generale de ScaLAPACK
!-------------------------------------------------------

SUBROUTINE init_scalapack(processor,communicator)

  use defs_scalapack
  use defs_basis
#if defined MPI && defined MPI2
 use mpi
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi, except_this_one => init_scalapack
!End of the abilint section

 implicit none
#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

  TYPE(processor_scalapack),INTENT(out)    :: processor
  INTEGER, INTENT(in)                       :: communicator

  TYPE(grid_scalapack)                    :: grid
  INTEGER                                   :: nbproc,myproc
  INTEGER                                   :: ierr


  CALL MPI_COMM_SIZE(communicator, nbproc, ierr)
  CALL MPI_COMM_RANK(communicator, myproc, ierr)

  CALL build_grid_scalapack(grid, nbproc, communicator)
  CALL build_processor_scalapack(processor, grid, myproc, communicator)


END SUBROUTINE init_scalapack


!-------------------------------------------------------
! General closing of ScaLAPACK
!-------------------------------------------------------

SUBROUTINE end_scalapack(processor)

  use defs_scalapack
  use defs_basis

  implicit none

  TYPE(processor_scalapack),INTENT(inout)    :: processor

  CALL BLACS_GRIDEXIT(processor%grid%ictxt)

  !CALL BLACS_EXIT(0)

END SUBROUTINE end_scalapack

!-------------------------------------------------------
! Initialisation of a SCALAPACK matrix (each proc initialize its own part of the matrix)
!-------------------------------------------------------

SUBROUTINE init_matrix_scalapack(matrix,nbli_global, &
     &                                      nbco_global,processor,tbloc)

  use defs_scalapack
  use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi, except_this_one => init_matrix_scalapack
!End of the abilint section

  implicit none

  TYPE(matrix_scalapack),INTENT(out)    :: matrix
  TYPE(processor_scalapack),INTENT(in),TARGET  :: processor
  INTEGER,INTENT(in)                     :: nbli_global,nbco_global
  INTEGER,INTENT(in),OPTIONAL            :: tbloc

  INTEGER, PARAMETER                :: SIZE_BLOCS = 40
  INTEGER             :: info,sizeb

  INTEGER :: NUMROC
  EXTERNAL NUMROC

  IF (PRESENT(tbloc)) THEN
     sizeb = tbloc
  ELSE
     sizeb  = SIZE_BLOCS
  END IF

  ! Records of the matrix type :
  matrix%processor => processor
  matrix%sizeb_blocs(1) = MIN(sizeb,nbli_global)
  matrix%sizeb_blocs(2) = MIN(sizeb,nbco_global)
  matrix%sizeb_global(1) = nbli_global
  matrix%sizeb_global(2) = nbco_global


  ! Size of the local buffer
  matrix%sizeb_local(1) = NUMROC(nbli_global,matrix%sizeb_blocs(1), &
       &                            processor%coords(1),0, &
       &                            processor%grid%dims(1))

  matrix%sizeb_local(2) = NUMROC(nbco_global,matrix%sizeb_blocs(2), &
       &                            processor%coords(2),0, &
       &                            processor%grid%dims(2))

  CALL idx_loc(matrix,matrix%sizeb_global(1),matrix%sizeb_global(2), &
       &       matrix%sizeb_local(1),matrix%sizeb_local(2))

  ! Initialisation of the SCALAPACK description of the matrix
  CALL DESCINIT(matrix%descript%tab, nbli_global, nbco_global, &
       &        matrix%sizeb_blocs(1), matrix%sizeb_blocs(2), 0,0 , &
       &        processor%grid%ictxt, MAX(1,matrix%sizeb_local(1)), &
       &        info)

  IF (info /= 0) THEN
     PRINT *,processor%myproc,'error initialisation matrix scalapack',info
  END IF

  ALLOCATE(matrix%buffer(matrix%sizeb_local(1),matrix%sizeb_local(2)))

  matrix%buffer(:,:) = (0._DP,0._DP)

END SUBROUTINE init_matrix_scalapack

!-------------------------------------------------------
! Copy of a Scalapack matrix
!-------------------------------------------------------
SUBROUTINE matrix_copy(source,destination)

  use defs_scalapack
  use defs_basis

  implicit none

  TYPE(matrix_scalapack),INTENT(in)  :: source
  TYPE(matrix_scalapack),INTENT(out) :: destination

  destination%processor => source%processor
  destination%sizeb_global = source%sizeb_global
  destination%sizeb_local = source%sizeb_local
  destination%sizeb_blocs = source%sizeb_blocs
  destination%descript%tab = source%descript%tab

  ALLOCATE(destination%buffer(destination%sizeb_local(1),destination%sizeb_local(2)))

  destination%buffer = source%buffer

  IF(ASSOCIATED(source%ipiv)) THEN
     ALLOCATE(destination%ipiv(SIZE(source%ipiv)))
     destination%ipiv = source%ipiv
  END IF
!    destination% = source%

END SUBROUTINE matrix_copy

!-------------------------------------------------------
! Destruction of the records of a SCALAPACK matrix
!-------------------------------------------------------
SUBROUTINE destruction_matrix_scalapack(matrix)

  use defs_scalapack
  use defs_basis

  implicit none

  TYPE(matrix_scalapack),INTENT(inout)    :: matrix

  NULLIFY(matrix%processor)
  matrix%sizeb_global = 0
  IF (ASSOCIATED(matrix%buffer)) THEN
     DEALLOCATE(matrix%buffer)
  ENDIF
  IF (ASSOCIATED(matrix%ipiv)) THEN
     DEALLOCATE(matrix%ipiv)
  ENDIF
  matrix%sizeb_blocs = 0
  matrix%sizeb_local = 0
  matrix%descript%tab = 0

END SUBROUTINE destruction_matrix_scalapack


!-------------------------------------------------------
!-------------------------------------------------------
! routines:
! access to the components of a matrix
!-------------------------------------------------------
!-------------------------------------------------------



!-------------------------------------------------------
! Access to a component thanks to its local indices
!-------------------------------------------------------

FUNCTION matrix_get_local(matrix,i,j)

  use defs_scalapack
  use defs_basis

  implicit none

  TYPE(matrix_scalapack),INTENT(in)    :: matrix
  INTEGER, INTENT(in)                   :: i,j
  COMPLEX(dp)                           :: matrix_get_local

  matrix_get_local = matrix%buffer(i,j)

END FUNCTION matrix_get_local

!-------------------------------------------------------
! Positioning of a component of a matrix thanks to its local indices
!-------------------------------------------------------

SUBROUTINE matrix_set_local(matrix,i,j,value)

  use defs_scalapack
  use defs_basis

  implicit none

  TYPE(matrix_scalapack),INTENT(out)   :: matrix

  INTEGER, INTENT(in)                   :: i,j
  COMPLEX(dp)                           :: value

  matrix%buffer(i,j) = value

END SUBROUTINE matrix_set_local

!-------------------------------------------------------
! Add a value to a component of the matrix thanks ot its local indices
!-------------------------------------------------------

SUBROUTINE matrix_add_local(matrix,i,j,value)

  use defs_scalapack
  use defs_basis

  implicit none

  TYPE(matrix_scalapack),INTENT(out)   :: matrix
  INTEGER, INTENT(in)                   :: i,j
  COMPLEX(dp)                           :: value

  matrix%buffer(i,j) = matrix%buffer(i,j) + value

END SUBROUTINE matrix_add_local

!-------------------------------------------------------
! Determins whether a term given by its local indices is locally stored
!-------------------------------------------------------

FUNCTION idx_processor_is_local(matrix,i,j)

  use defs_scalapack
  use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi, except_this_one => idx_processor_is_local
!End of the abilint section

  implicit none

  TYPE(matrix_scalapack),INTENT(in)    :: matrix
  INTEGER, INTENT(in)                   :: i,j

  LOGICAL ::idx_processor_is_local
  idx_processor_is_local = (matrix%processor%coords(1) == &
       &                      idx_processor_concerned(matrix,i,1)) &
       &                    .AND. (matrix%processor%coords(2) == &
       &                           idx_processor_concerned(matrix,j,2))

END FUNCTION idx_processor_is_local

!-------------------------------------------------------
! Determine the number of the column/row  of targeted processors thanks
! to a global column/row number.
!-------------------------------------------------------

FUNCTION idx_processor_concerned(matrix,idx,lico)

  use defs_scalapack
  use defs_basis

  implicit none

  TYPE(matrix_scalapack),INTENT(in)    :: matrix
  INTEGER, INTENT(in)                   :: idx,lico !lico= row or column

  INTEGER :: idx_processor_concerned
  INTEGER :: INDXG2P
  EXTERNAL INDXG2P


  idx_processor_concerned = INDXG2P(idx,matrix%sizeb_blocs(lico),0,0, &
       &                               matrix%processor%grid%dims(lico))
END FUNCTION idx_processor_concerned

!-------------------------------------------------------
! Determination of the local indices of a term of a matrix with respect
! to its lobal indices independently of the proc 
!-------------------------------------------------------

SUBROUTINE idx_loc(matrix,i,j,iloc,jloc)

  use defs_scalapack
  use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi, except_this_one => idx_loc
!End of the abilint section

  implicit none

  TYPE(matrix_scalapack),INTENT(in)    :: matrix
  INTEGER, INTENT(in)                   :: i,j
  INTEGER, INTENT(out)                  :: iloc,jloc

  INTEGER :: NUMROC
  EXTERNAL NUMROC
  iloc = glob_loc(matrix,i,1)

  jloc = glob_loc(matrix,j,2)

END SUBROUTINE idx_loc

!-------------------------------------------------------
! 
!-------------------------------------------------------
FUNCTION glob_loc(matrix,idx,lico)

  use defs_scalapack
  use defs_basis

  implicit none

  TYPE(matrix_scalapack),INTENT(in)    :: matrix
  INTEGER, INTENT(in)                   :: idx, lico

  INTEGER :: glob_loc

  INTEGER :: NUMROC
  EXTERNAL NUMROC

  glob_loc = NUMROC(idx,matrix%sizeb_blocs(lico), &
       &        matrix%processor%coords(lico),0, &
       &        matrix%processor%grid%dims(lico))


END FUNCTION glob_loc

!-------------------------------------------------------
! Determination of the global indices of a term of the matrix with respect
! to its local indices
!-------------------------------------------------------

SUBROUTINE idx_glob(matrix,iloc,jloc,i,j)

  use defs_scalapack
  use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi, except_this_one => idx_glob
!End of the abilint section

  implicit none

  TYPE(matrix_scalapack),INTENT(in)    :: matrix
  INTEGER, INTENT(out)                  :: i,j
  INTEGER, INTENT(in)                   :: iloc,jloc

  INTEGER :: nbcycli,nbcycco,resteli,resteco,nblocsli,nblocsco
  i = loc_glob(matrix,matrix%processor,iloc,1)
  j = loc_glob(matrix,matrix%processor,jloc,2)

END SUBROUTINE idx_glob

!-------------------------------------------------------
! Determination of the global index from a local index (row or column)
! as a function of a given processor
!-------------------------------------------------------
FUNCTION loc_glob(matrix,proc,idx,lico)

  use defs_scalapack
  use defs_basis

  implicit none

  TYPE(matrix_scalapack),INTENT(in)    :: matrix
  TYPE(processor_scalapack),INTENT(in) :: proc
  INTEGER, INTENT(in)                   :: idx,lico

  INTEGER :: loc_glob

  INTEGER :: nbcyc,reste,nblocs

  nbcyc = INT((idx-1)/matrix%sizeb_blocs(lico))
  reste = MOD(idx-1,matrix%sizeb_blocs(lico))
  nblocs = nbcyc*proc%grid%dims(lico)+ &
       & proc%coords(lico)

  loc_glob = nblocs * matrix%sizeb_blocs(lico) + reste + 1

END FUNCTION loc_glob

!-------------------------------------------------------
! Query of a term of the matrix from its global indices
! WARNING : this is valid only if the calling processor
! effectively possesses locally the term.
!-------------------------------------------------------

FUNCTION matrix_get_global(matrix,i,j)

  use defs_scalapack
  use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi, except_this_one => matrix_get_global
!End of the abilint section

  implicit none

  TYPE(matrix_scalapack),INTENT(in)    :: matrix
  INTEGER, INTENT(in)                   :: i,j
  COMPLEX(dp)                           :: matrix_get_global

  INTEGER :: iloc,jloc
  CALL idx_loc(matrix,i,j,iloc,jloc)

  matrix_get_global = matrix_get_local(matrix,iloc,jloc)

END FUNCTION matrix_get_global

!-------------------------------------------------------
! Localication of the value of a term of a matrix thanks to its global indices.
! WARNING : this is valid only if the calling processor 
! effectively possesses locally the term.
!-------------------------------------------------------

SUBROUTINE matrix_set_global(matrix,i,j,value)

  use defs_scalapack
  use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi, except_this_one => matrix_set_global
!End of the abilint section

  implicit none

  TYPE(matrix_scalapack),INTENT(inout) :: matrix
  INTEGER, INTENT(in)                   :: i,j
  COMPLEX(dp)                           :: value

  INTEGER :: iloc,jloc

  CALL idx_loc(matrix,i,j,iloc,jloc)

  CALL matrix_set_local(matrix,iloc,jloc,value)

END SUBROUTINE matrix_set_global

!-------------------------------------------------------
! Add a value to a term of the matrix thanks to its global indices.
! WARNING : this is valid only if the calling processor 
! effectively possesses locally the term.
!-------------------------------------------------------

SUBROUTINE matrix_add_global(matrix,i,j,value)

  use defs_scalapack
  use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi, except_this_one => matrix_add_global
!End of the abilint section

  implicit none

  TYPE(matrix_scalapack),INTENT(inout) :: matrix
  INTEGER, INTENT(in)                   :: i,j
  COMPLEX(dp)                           :: value

  INTEGER :: iloc,jloc

  CALL idx_loc(matrix,i,j,iloc,jloc)

  CALL matrix_add_local(matrix,iloc,jloc,value)

END SUBROUTINE matrix_add_global


!-------------------------------------------------------
! Checking routine of a SCALAPACK matrix with respect to a complete matrix
!-------------------------------------------------------

SUBROUTINE matrix_check_global(matrix,reference)

  use defs_scalapack
  use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi, except_this_one => matrix_check_global
!End of the abilint section

  implicit none

  TYPE(matrix_scalapack),INTENT(in)   :: matrix
  COMPLEX(dp),DIMENSION(:,:)           :: reference

  INTEGER :: i,j,iglob,jglob
  COMPLEX(dp) :: diff,z
  REAL(dp) :: err
  INTEGER :: cptr
  err = 0._DP
  cptr = 0

  DO i=1,matrix%sizeb_local(1)
     DO j=1,matrix%sizeb_local(2)
        CALL idx_glob(matrix,i,j,iglob,jglob)
        z = matrix_get_local(matrix,i,j)
        diff = z-reference(iglob,jglob)
        err = err + ABS(diff)/MAX(ABS(z),MAX(1.E-35,ABS(reference(iglob,jglob))))
        cptr = cptr + 1
     END DO
  END DO

  IF (cptr /= 0) THEN
     PRINT *,matrix%processor%myproc,"error Linf matrix scalapack", &
          &  err,"on",cptr,"terms"
  END IF
END SUBROUTINE matrix_check_global

!-------------------------------------------------------
! Check a segment of SCALAPACK matrix with respect to a vector
!-------------------------------------------------------

SUBROUTINE matrix_check_global_vector(matrix,reference,decli,nom)

  use defs_scalapack
  use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi, except_this_one => matrix_check_global_vector
!End of the abilint section

  implicit none

  TYPE(matrix_scalapack),INTENT(in)   :: matrix
  COMPLEX(dp),DIMENSION(:)             :: reference
  INTEGER, OPTIONAL                    :: decli
  CHARACTER(len=*),OPTIONAL            :: nom

  INTEGER :: i,j,deb
  COMPLEX(dp) :: diff,z
  REAL(dp) :: err
  INTEGER :: cptr
  CHARACTER(len=150)            :: nomvec
  IF (PRESENT(decli)) THEN
     deb = decli
  ELSE
     deb = 0
  END IF
  IF (PRESENT(nom)) THEN
     nomvec=nom
  ELSE
     nomvec='(vector)'
  END IF

  err = 0._DP
  cptr = 0

  j=1
  IF (matrix%processor%coords(2) == &
       & idx_processor_concerned(matrix,j,2)) THEN
     DO i=1,MIN(UBOUND(reference,1),matrix%sizeb_global(1)-deb)
        IF (matrix%processor%coords(1) == &
             & idx_processor_concerned(matrix,deb+i,1)) THEN
           z = matrix_get_global(matrix,deb+i,j)
           diff = z-reference(i)
           err = err + ABS(diff)/MAX(ABS(z),MAX(1.E-35,ABS(reference(i))))
           !print *,'error',i,j
           cptr = cptr + 1
        END IF
     END DO
  END IF
  IF (cptr /= 0) THEN
     PRINT *,matrix%processor%myproc,"error L1 ",TRIM(nomvec), &
          &  err/cptr,'on',cptr,'stored terms'
  END IF
END SUBROUTINE matrix_check_global_vector


!-----------------------------------------------------------------
! Routine to fill a SCALAPACK matrix with respect to a full matrix
!-----------------------------------------------------------------

SUBROUTINE matrix_from_global(matrix,reference)

  use defs_scalapack
  use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi, except_this_one => matrix_from_global
!End of the abilint section

  implicit none

  TYPE(matrix_scalapack),INTENT(inout)  :: matrix
  REAL(dp),DIMENSION(:)                  :: reference
  COMPLEX(dp)::val

  INTEGER :: i,j,iglob,jglob,ind
  REAL(dp) :: err
  INTEGER :: cptr

!    err = 0._DP
!    cptr = 0

  DO i=1,matrix%sizeb_local(1)
     DO j=1,matrix%sizeb_local(2)
        CALL idx_glob(matrix,i,j,iglob,jglob)

        ind = jglob*(jglob-1)+2*iglob-1
        val = dcmplx(reference(ind),reference(ind+1))
        CALL matrix_set_local(matrix,i,j,val)

!          cptr = cptr + 1
     END DO
  END DO

!    IF (cptr /= 0) THEN
!       PRINT *,matrix%processor%myproc,"error Linf matrix scalapack", &
!            &  err,"on",cptr,"terms"
!    END IF

END SUBROUTINE matrix_from_global

!-------------------------------------------------------------------
! Routine to fill a full matrix with respect to a SCALAPACK matrix
!-------------------------------------------------------------------
SUBROUTINE matrix_to_reference(matrix,reference)

  use defs_scalapack
  use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi, except_this_one => matrix_to_reference
!End of the abilint section

  implicit none

  TYPE(matrix_scalapack),INTENT(in)        :: matrix
  REAL(dp),DIMENSION(:,:),INTENT(inout)     :: reference

  INTEGER  :: i,j,iglob,jglob,ind
  REAL(dp) :: err
  INTEGER  :: cptr

!    err = 0._DP
!    cptr = 0

  DO i=1,matrix%sizeb_local(1)
     DO j=1,matrix%sizeb_local(2)
        CALL idx_glob(matrix,i,j,iglob,jglob)

        ind=(iglob-1)*2+1
        reference(ind,jglob)   = REAL(matrix%buffer(i,j))
        reference(ind+1,jglob) = IMAG(matrix%buffer(i,j))

!          cptr = cptr + 1
     END DO
  END DO

!    IF (cptr /= 0) THEN
!       PRINT *,matrix%processor%myproc,"error Linf matrix scalapack", &
!            &  err,"on",cptr,"terms"
!    END IF

END SUBROUTINE matrix_to_reference

!-------------------------------------------------------
! Storage of a vector in a column of a SCALAPACK matrix
!-------------------------------------------------------
SUBROUTINE matrix_storage_local_vector(matrix,array,decli,decco)

  use defs_scalapack
  use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi, except_this_one => matrix_storage_local_vector
!End of the abilint section

  implicit none

  COMPLEX(dp),DIMENSION(:),INTENT(in)       :: array
  TYPE(matrix_scalapack),INTENT(inout)     :: matrix
  INTEGER, INTENT(in)                       :: decli ! line shift in the matrix
  INTEGER, INTENT(in)                       :: decco ! column shift in the matrix 

  INTEGER :: numli


  ! If the term is stored locally
  IF (matrix%processor%coords(2) == &
       & idx_processor_concerned(matrix,decco+1,2)) THEN

     DO numli = 1, UBOUND(array,1)
        IF (matrix%processor%coords(1) == &
             & idx_processor_concerned(matrix,decli+numli,1)) THEN

           CALL matrix_set_global(matrix,decli+numli,decco+1, &
                &                  array(numli))

        END IF
     END DO
  END IF

END SUBROUTINE matrix_storage_local_vector

!-------------------------------------------------------
! Extract the components of a column of the matrix to a vector.
! The components that a re no locally stored by the matrix are left unchanged.
!-------------------------------------------------------
SUBROUTINE matrix_extract_vector(matrix,array,decli,decco,finli)

  use defs_scalapack
  use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi, except_this_one => matrix_extract_vector
!End of the abilint section

  implicit none

  COMPLEX(dp),DIMENSION(:),INTENT(out)      :: array
  TYPE(matrix_scalapack),INTENT(in)        :: matrix
  INTEGER, INTENT(in)                       :: decli ! line shift in the matrix 
  INTEGER, INTENT(in)                       :: decco ! column shift in the matrix 
  INTEGER, INTENT(in), OPTIONAL             :: finli

  INTEGER :: i,iglob,jglob,maxli,debli,debco
  COMPLEX(dp) :: x
  IF (PRESENT(finli)) THEN
     maxli = finli
  ELSE
     maxli = matrix%sizeb_global(1)
  END IF

  i = decli+1

  DO WHILE (.NOT. matrix%processor%coords(1) == &
       & idx_processor_concerned(matrix,i,1))
     i=i+1
  END DO

  IF (i .GT. maxli) RETURN

  CALL idx_loc(matrix,i,decco+1,debli,debco)

  ! If one is concerned by the jglob column
  IF (matrix%processor%coords(2) == &
       & idx_processor_concerned(matrix,decco+1,2)) THEN

     DO i = debli,matrix%sizeb_local(1)
        x = matrix_get_local(matrix,i,debco)
        CALL idx_glob(matrix,i,debco,iglob,jglob)

        IF (iglob .GT. maxli) THEN
           RETURN
        END IF
        array(iglob-decli) = x
     END DO
  END IF

END SUBROUTINE matrix_extract_vector

!-------------------------------------------------------
! A=LU factorization of a (square) matrix
!-------------------------------------------------------
SUBROUTINE matrix_pzgetrf(matrix)

  use defs_scalapack
  use defs_basis

  implicit none

  TYPE(matrix_scalapack),INTENT(inout)        :: matrix

  INTEGER                                        :: info, nipiv

  EXTERNAL PZGETRF

  IF (.NOT. ASSOCIATED(matrix%ipiv)) THEN
     nipiv = matrix%sizeb_global(1) + matrix%sizeb_blocs(1)
     ALLOCATE(matrix%ipiv(nipiv))
     matrix%ipiv(:) = 0
  END IF

  CALL PZGETRF(matrix%sizeb_global(1), matrix%sizeb_global(2), &
       &       matrix%buffer(1,1),1,1,matrix%descript%tab, &
       &       matrix%ipiv, info)

  IF (info /= 0) THEN
     PRINT *,matrix%processor%myproc,'error pzgetrf',info
  END IF
END SUBROUTINE matrix_pzgetrf

!-------------------------------------------------------
! Determination of the solution of a linear system whose matrix is factorized
!-------------------------------------------------------
SUBROUTINE matrix_pzgetrs(matrix,vector)

  use defs_scalapack
  use defs_basis

  implicit none

  TYPE(matrix_scalapack),INTENT(in)        :: matrix
  TYPE(matrix_scalapack),INTENT(inout)        :: vector

  INTEGER :: info

  CALL PZGETRS('N',matrix%sizeb_global(1),vector%sizeb_global(2), &
       &       matrix%buffer,1,1,matrix%descript%tab,matrix%ipiv, &
       &       vector%buffer,1,1,vector%descript%tab,info)

  IF (info /= 0) THEN
     PRINT *,matrix%processor%myproc,'error pzgetrs',info
  END IF
END SUBROUTINE matrix_pzgetrs

!-------------------------------------------------------
! A=LLT factorisation of a (hermitian) matrix
!-------------------------------------------------------
SUBROUTINE matrix_pzpotrf(matrix)

  use defs_scalapack
  use defs_basis

  implicit none

  TYPE(matrix_scalapack),INTENT(inout)        :: matrix

  INTEGER                                        :: info, nipiv

  EXTERNAL PZPOTRF

  CALL PZPOTRF('U',matrix%sizeb_global(1), &
       &       matrix%buffer(1,1),1,1,matrix%descript%tab,info)

  IF (info /= 0) THEN
     PRINT *,matrix%processor%myproc,'error pzpotrf',info
  END IF
END SUBROUTINE matrix_pzpotrf

!-------------------------------------------------------
! Determination of the solution of a linear system whose (hermitian) matrix is factorized
!-------------------------------------------------------
SUBROUTINE matrix_pzpotrs(matrix,vector)

  use defs_scalapack
  use defs_basis

  implicit none

  TYPE(matrix_scalapack),INTENT(in)        :: matrix
  TYPE(matrix_scalapack),INTENT(inout)        :: vector

  INTEGER :: info

  CALL PZPOTRS('U',matrix%sizeb_global(1),vector%sizeb_global(2), &
       &       matrix%buffer,1,1,matrix%descript%tab, &
       &       vector%buffer,1,1,vector%descript%tab,info)

  IF (info /= 0) THEN
     PRINT *,matrix%processor%myproc,'error pzpotrs',info
  END IF
END SUBROUTINE matrix_pzpotrs

!-------------------------------------------------------
! Extended matrix*matrix product
! C := alpha*A*B - beta*C
!
! For a simple matrix vector product, one can simply pass 
! alpha = (1.,0.) and beta (0.,0.)
!-------------------------------------------------------
SUBROUTINE matrix_pzgemm(matrix1,alpha,matrix2,beta,results)

  use defs_scalapack
  use defs_basis

  implicit none

  TYPE(matrix_scalapack),INTENT(in)        :: matrix1,matrix2
  TYPE(matrix_scalapack),INTENT(inout)     :: results
  COMPLEX(dp), intent(in)                   :: alpha, beta

  CALL PZGEMM('N','N',matrix1%sizeb_global(1),matrix2%sizeb_global(2),&
       &      matrix1%sizeb_global(2),alpha,matrix1%buffer,1,1, &
       &      matrix1%descript%tab,matrix2%buffer,1,1, &
       &      matrix2%descript%tab,beta,results%buffer,1,1, &
       &      results%descript%tab)

END SUBROUTINE matrix_pzgemm


!-------------------------------------------------------
!  Calculation of eigenvalues and eigenvectors
!-------------------------------------------------------
SUBROUTINE matrix_pzheevx(processor,matrix,results,eigen,communicator)

  use defs_scalapack
  use defs_basis

#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none
#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

  TYPE(processor_scalapack),INTENT(in)       :: processor
  TYPE(matrix_scalapack),INTENT(in)          :: matrix
  TYPE(matrix_scalapack),INTENT(inout)       :: results
  DOUBLE PRECISION,DIMENSION(:),INTENT(inout) :: eigen
  INTEGER,INTENT(in)  :: communicator

  INTEGER            :: LRWORK,LIWORK,LWORK,INFO
 
  INTEGER         , dimension(3) :: IWORK_tmp
  DOUBLE PRECISION, dimension(3) :: RWORK_tmp
  COMPLEX(dp)     , dimension(3) :: WORK_tmp

  INTEGER         , allocatable  :: IWORK(:)
  DOUBLE PRECISION, allocatable  :: RWORK(:)
  COMPLEX(dp)     , allocatable  :: WORK(:)

  INTEGER,          allocatable :: ICLUSTR(:)
  INTEGER,          allocatable :: IFAIL(:)
  DOUBLE PRECISION, allocatable :: GAP(:)

  DOUBLE PRECISION, PARAMETER :: ABSTOL=-1.D+0,ORFAC=-1.D+0
  INTEGER,          PARAMETER :: IZERO=0

  INTEGER ::  M,NZ,IA,JA,IZ,JZ,ierr,TWORK_tmp(3),TWORK(3)


  INFO   = 0   

  ! Allocation of the variables for the results of the calculations
  allocate(IFAIL(matrix%sizeb_global(2)))
  allocate(ICLUSTR(2*processor%grid%dims(1)*processor%grid%dims(2)))
  allocate(GAP(processor%grid%dims(1)*processor%grid%dims(2)))
  
  ! Get the size of the work arrays
  CALL PZHEEVX('V','A','U',&
       &      matrix%sizeb_global(2),&
       &      matrix%buffer,1,1,matrix%descript%tab, &
       &      ZERO,ZERO,IZERO,IZERO,ABSTOL,&
       &      m,nz,eigen,ORFAC, &
       &      results%buffer,1,1,results%descript%tab, &
       &      WORK_tmp,-1,RWORK_tmp,-1,IWORK_tmp,-1,&
       &      IFAIL,ICLUSTR,GAP,INFO)

  TWORK_tmp(1) = IWORK_tmp(1)
  TWORK_tmp(2) = INT(RWORK_tmp(1))
  TWORK_tmp(3) = INT(WORK_tmp(1))

 !! Get the maximum of the size of the work arrays processor%comm
  CALL MPI_ALLREDUCE(TWORK_tmp,TWORK,3,MPI_INTEGER,MPI_MAX,communicator,ierr)

  LIWORK = TWORK(1)
  LRWORK = TWORK(2)
  LWORK  = TWORK(3)

 ! Allocation of the work arrays
  allocate(IWORK(LIWORK))
  allocate(WORK(LWORK))
  allocate(RWORK(LRWORK))

  ! Call the calculation routine
  CALL PZHEEVX('V','A','U',&
       &      matrix%sizeb_global(2),&
       &      matrix%buffer,1,1,matrix%descript%tab, &
       &      ZERO,ZERO,IZERO,IZERO,ABSTOL,&
       &      m,nz,eigen,ORFAC, &
       &      results%buffer,1,1,results%descript%tab, &
       &      WORK,LWORK,RWORK,LRWORK,IWORK,LIWORK,&
       &      IFAIL,ICLUSTR,GAP,INFO)


  deallocate(IFAIl,ICLUSTR,GAP,WORK,RWORK,IWORK)

END SUBROUTINE matrix_pzheevx

!-------------------------------------------------------
!  Calculation of eigenvalues and eigenvectors
!  A * X = lambda * B * X 
!-------------------------------------------------------
SUBROUTINE matrix_pzhegvx(processor,matrix1,matrix2,results,eigen,communicator)

  use defs_scalapack
  use defs_basis

#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none
#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

  TYPE(processor_scalapack),INTENT(in)       :: processor
  TYPE(matrix_scalapack),INTENT(in)          :: matrix1,matrix2
  TYPE(matrix_scalapack),INTENT(inout)       :: results
  DOUBLE PRECISION,DIMENSION(:),INTENT(inout) :: eigen

  INTEGER,INTENT(in)  :: communicator

  INTEGER            :: LRWORK,LIWORK,LWORK,INFO
 
  INTEGER         , dimension(3) :: IWORK_tmp
  DOUBLE PRECISION, dimension(3) :: RWORK_tmp
  COMPLEX(dp)     , dimension(3) :: WORK_tmp

  INTEGER         , allocatable  :: IWORK(:)
  DOUBLE PRECISION, allocatable  :: RWORK(:)
  COMPLEX(dp)     , allocatable  :: WORK(:)


  INTEGER,          allocatable :: ICLUSTR(:)
  INTEGER,          allocatable :: IFAIL(:)
  DOUBLE PRECISION, allocatable :: GAP(:)

  DOUBLE PRECISION, PARAMETER :: ABSTOL=-1.D+0,ORFAC=-1.D+0
  INTEGER         , PARAMETER :: IZERO=0

  INTEGER ::  M,NZ,IA,JA,IZ,JZ,ierr,TWORK_tmp(3),TWORK(3)


  INFO   = 0   

  ! Allocate the arrays for the results of the calculation
  allocate(IFAIL(matrix1%sizeb_global(2)))
  allocate(ICLUSTR(2*processor%grid%dims(1)*processor%grid%dims(2)))
  allocate(GAP(processor%grid%dims(1)*processor%grid%dims(2)))

  ! Get the size of the work arrays
  CALL PZHEGVX(1,'V','A','U',&
       &      matrix1%sizeb_global(2),&
       &      matrix1%buffer,1,1,matrix1%descript%tab, &
       &      matrix2%buffer,1,1,matrix2%descript%tab, &
       &      ZERO,ZERO,IZERO,IZERO,ABSTOL,&
       &      m,nz,eigen,ORFAC, &
       &      results%buffer,1,1,results%descript%tab, &
       &      WORK_tmp,-1,RWORK_tmp,-1,IWORK_tmp,-1,&
       &      IFAIL,ICLUSTR,GAP,INFO)


  TWORK_tmp(1) = IWORK_tmp(1)
  TWORK_tmp(2) = INT(RWORK_tmp(1))
  TWORK_tmp(3) = INT(WORK_tmp(1))
 
 ! Get the maximum of sizes of the work arrays processor%comm 
  CALL MPI_ALLREDUCE(TWORK_tmp,TWORK,3,MPI_INTEGER,MPI_MAX,communicator,ierr)

  LIWORK = TWORK(1)
  LRWORK = TWORK(2)
  LWORK  = TWORK(3)

 ! Allocate the work arrays
  allocate(IWORK(LIWORK))
  allocate(WORK(LWORK))
  allocate(RWORK(LRWORK))

  ! Call the calculation routine
  CALL PZHEGVX(1,'V','A','U',&
       &      matrix1%sizeb_global(2),&
       &      matrix1%buffer,1,1,matrix1%descript%tab, &
       &      matrix2%buffer,1,1,matrix2%descript%tab, &
       &      ZERO,ZERO,IZERO,IZERO,ABSTOL,&
       &      m,nz,eigen,ORFAC, &
       &      results%buffer,1,1,results%descript%tab, &
       &      WORK,LWORK,RWORK,LRWORK,IWORK,LIWORK,&
       &      IFAIL,ICLUSTR,GAP,INFO)

  deallocate(IFAIl,ICLUSTR,GAP,IWORK,WORK,RWORK)

END SUBROUTINE matrix_pzhegvx

#else
   SUBROUTINE NO_SCALAPACK


    implicit none
   END SUBROUTINE NO_SCALAPACK
#endif
!!***
