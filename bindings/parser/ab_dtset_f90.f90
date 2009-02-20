module ab_dtset

  use defs_basis
  use defs_datatypes
  use interfaces_13iovars

  implicit none

  private

  ! We store here a list of dtset arrays to be able to
  ! parse several ABINIT files without freeing it.
  ! The simplest portable way to do it, is to create
  ! a list of dtsets arrays and to use the list index
  ! as an identifier that can be given to the other languages.
  type, private :: dtsets_list
     type(dtsets_list),  pointer :: next => null()
     type(dataset_type), pointer :: dtsets(:)
  end type dtsets_list
  type(dtsets_list), pointer :: my_dtsets => null()
  integer :: nb_dtsets = 0

  ! Error codes
  integer, parameter, public :: AB_NO_ERROR         = 0
  integer, parameter, public :: AB_ERROR_DTSET_OBJ  = 1
  integer, parameter, public :: AB_ERROR_DTSET_ATT  = 2
  integer, parameter, public :: AB_ERROR_DTSET_ID   = 3
  integer, parameter, public :: AB_ERROR_DTSET_SIZE = 4

  logical, private, parameter :: AB_DBG = .false.

  include "dtset_f90.inc"

  public :: ab_dtset_new
  public :: ab_dtset_new_from_string
  public :: ab_dtset_free

  public :: ab_dtset_get_ndtset
  public :: ab_dtset_get_integer
  public :: ab_dtset_get_real
  public :: ab_dtset_get_shape
  public :: ab_dtset_get_integer_array
  public :: ab_dtset_get_real_array

contains

  subroutine ab_dtset_new_from_string(dtsetsId, instr, len)
    integer, intent(out) :: dtsetsId
    integer, intent(in) :: len
    character(len = len), intent(in) :: instr

    character(len = strlen) :: string
    integer :: lenstr, ndtset
    integer :: marr, tread
    character(len = 30) :: token
    integer :: intarr(1)
    real(dp) :: dprarr(1)
    character(len=500) :: message

    if (len > strlen) then
       dtsetsId = 0
       return
    end if

    write(string,*) instr

    !To make case-insensitive, map characters of string to upper case:
    call inupper(string(1:len))

    !Might import data from CML file(s) into string
    !Need string_raw to deal properly with CML filenames
    lenstr = len
    call importcml(lenstr, instr, string, len)

    !6) Take ndtset from the input string, then allocate
    !the arrays whose dimensions depends only on ndtset and msym.

    ndtset=0 ; marr=1
    token = 'ndtset'
    call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),token,tread,'INT')
    if(tread==1) ndtset=intarr(1)
    !Check that ndtset is not negative
    if (ndtset<0 .or. ndtset>99) then
       write(message, '(a,a,a,a,i12,a,a,a,a)' ) ch10,&
            &  ' abinit : ERROR -',ch10,&
            &  '  Input ndtset must be non-negative and < 100, but was ',ndtset,ch10,&
            &  '  This is not allowed.  ',ch10,&
            &  '  Action : modify ndtset in the input file.'
       call wrtout(06,  message,'COLL')
       call leave_new('COLL')
    end if

    call ab_dtset_load(dtsetsId, string, lenstr, ndtset)
  end subroutine ab_dtset_new_from_string

  subroutine ab_dtset_new(dtsetsId, filename, n)
    integer, intent(out) :: dtsetsId
    integer, intent(in) :: n
    character(len = n), intent(in) :: filename

    character(len = strlen) :: string
    integer :: lenstr, ndtset

    if (AB_DBG) write(0,*) "AB module: read '", filename, "' to string."
    call parsefile(filename, lenstr, ndtset, string)
    if (AB_DBG) write(0,*) "AB module: read OK, string length ", lenstr

    call ab_dtset_load(dtsetsId, string, lenstr, ndtset)
  end subroutine ab_dtset_new

  subroutine ab_dtset_load(dtsetsId, string, lenstr, ndtset)
    integer, intent(out) :: dtsetsId
    character(len = strlen), intent(inout) :: string
    integer, intent(in) :: lenstr, ndtset

    type(dtsets_list), pointer :: token
    type(pspheader_type), pointer :: pspheads(:)
    integer :: istatr, istatshft, jdtset
    integer :: ndtset_alloc
    type(MPI_type) :: mpi_enreg
    integer :: npsp, ii, idtset, msym, usepaw, dmatpuflag
    integer :: mxnatom, mxntypat, mxlpawu, mxmband_upper, mxnatpawu, &
         & mxnatsph, mxncenter, mxnconeq, mxnkptgw, mxnatvshift, &
         & mxnkpt, mxnorb, mxnnos, mxnqptdm, mxnsppol, mxnsym, mxnspinor
    integer,allocatable :: bravais_(:,:),mband_upper_(:)
    real(dp) :: zion_max

    ! We allocate a new list token and prepend it.
    if (AB_DBG) write(0,*) "AB module: allocate a new object."
    allocate(token)
    token%next => my_dtsets
    my_dtsets => token
    nb_dtsets = nb_dtsets + 1
    dtsetsId = nb_dtsets

    ndtset_alloc=ndtset ; if(ndtset==0)ndtset_alloc=1
    allocate(token%dtsets(0:ndtset_alloc))
    do idtset = 0, size(token%dtsets) - 1, 1
       call dtsetInit(token%dtsets(idtset))
    end do
    if (AB_DBG) write(0,*) "AB module: allocation OK at ", dtsetsId

    if (AB_DBG) write(0,*) "AB module: call invars0()."
    mpi_enreg%nproc = 1
    mpi_enreg%me = 0
    !7) Continue to analyze the input string, get upper dimensions,
    !and allocate the remaining arrays.
    call invars0(token%dtsets,ab_out,istatr,istatshft,lenstr,&
         & mxnatom,mxntypat,ndtset,ndtset_alloc,npsp,string)

    !8) Finish to read the "file" file completely, as npsp is known,
    !and also initialize pspheads, that contains the important information
    !from the pseudopotential headers, as well as the psp filename
    usepaw=0
    nullify(pspheads)
    allocate(pspheads(npsp))
    ! No psp files are given, we put default values into pspheads.
    pspheads(:)%zionpsp = 1
    pspheads(:)%pspxc   = token%dtsets(1)%ixc
    pspheads(:)%pspso   = 0
    pspheads(:)%xccc    = 0
    do idtset=0,ndtset_alloc
       token%dtsets(idtset)%usepaw=usepaw
    end do

    !Take care of other dimensions, and part of the content of dtsets
    !that is or might be needed early.
    !zion_max=maxval(pspheads(1:npsp)%zionpsp) ! This might not work properly with HP compiler
    zion_max=pspheads(1)%zionpsp
    do ii=1,npsp
       zion_max=max(pspheads(ii)%zionpsp,zion_max)
    end do
    if (AB_DBG) write(0,*) "AB module: OK."

    allocate(bravais_(11,0:ndtset_alloc))
    allocate(mband_upper_ (  0:ndtset_alloc))

    if (AB_DBG) write(0,*) "AB module: call invars1m()."
    !Here, the default msym.
    msym=384
    call invars1m(bravais_,dmatpuflag,token%dtsets,ab_out,lenstr,mband_upper_,mpi_enreg,&
         & msym,mxlpawu,mxmband_upper,mxnatom,mxnatpawu,mxnatsph,mxnatvshift,mxncenter,mxnconeq,mxnkptgw,&
         & mxnkpt,mxnorb,mxnnos,mxnqptdm,mxnspinor,mxnsppol,mxnsym,ndtset,ndtset_alloc,string,&
         & zion_max)
    if (AB_DBG) write(0,*) "AB module: OK."

    !9) Provide defaults for the variables that have not yet been initialized.
    if (AB_DBG) write(0,*) "AB module: call indefo()."
    call indefo(token%dtsets,mpi_enreg,ndtset_alloc)

    !If all the pseudopotentials have the same pspxc, override the default
    !value for dtsets 1 to ndtset
    !if(minval(abs((pspheads(1:npsp)%pspxc-pspheads(1)%pspxc)))==0)then
    !   token%dtsets(1:ndtset_alloc)%ixc=pspheads(1)%pspxc
    !end if
    if (AB_DBG) write(0,*) "AB module: OK."

    !10) Call the main input routine.
    if (AB_DBG) write(0,*) "AB module: call invars2()."
    do idtset = 1, ndtset_alloc, 1
       jdtset=token%dtsets(idtset)%jdtset ; if(ndtset==0)jdtset=0
       call invars2(bravais_(:, idtset),token%dtsets(idtset),ab_out,jdtset,lenstr,&
            &  mband_upper_(idtset),msym,npsp,pspheads,string,usepaw,pspheads(1:npsp)%zionpsp)
    end do
!!$    call invars2m(bravais_,dtsets,ab_out,lenstr,&
!!$         & mband_upper_,mpi_enreg,msym,ndtset,ndtset_alloc,npsp,pspheads,string)
    deallocate(bravais_, mband_upper_)
    if (AB_DBG) write(0,*) "AB module: OK."

    deallocate(pspheads)

  end subroutine ab_dtset_load

  subroutine ab_dtset_get_from_id(token, id)
    type(dtsets_list), pointer :: token
    integer, intent(in) :: id

    type(dtsets_list), pointer :: tmpLst
    integer :: i

    if (AB_DBG) write(0,*) "AB module: request list element ", id, "from ", nb_dtsets
    if (id > nb_dtsets .or. id <= 0) then
       nullify(token)
       if (AB_DBG) write(0,*) " | but not in list, nullify token."
       return
    end if
    ! List element are prepended so element id is at (nb - id) position.
    tmpLst => my_dtsets
    do i = 1, nb_dtsets - id, 1
       tmpLst => tmpLst%next
    end do
    if (associated(tmpLst%dtsets)) then
       token => tmpLst
    else
       nullify(token)
    end if
  end subroutine ab_dtset_get_from_id

  subroutine ab_dtset_free(dtsetsId)
    integer, intent(in) :: dtsetsId

    type(dtsets_list), pointer :: token
    integer :: idtset

    if (AB_DBG) write(0,*) "AB module: Free request on dataset array ", dtsetsId
    nullify(token)
    call ab_dtset_get_from_id(token, dtsetsId)
    if (associated(token)) then
       if (associated(token%dtsets)) then
          if (AB_DBG) write(0,*) " | ", size(token%dtsets), "dtsets found."
          do idtset = 0, size(token%dtsets) - 1, 1
             if (AB_DBG) write(0,*) " | free dtset ", idtset
             call dtsetFree(token%dtsets(idtset))
             if (AB_DBG) write(0,*) " | free OK"
          end do
          deallocate(token%dtsets)
          nullify(token%dtsets)
          if (AB_DBG) write(0,*) " | general free OK"
!!$       ! We remove token from the list.
!!$       call _dtset_get_from_id(prev, dtsetsId - 1)
!!$       if (associated(prev)) then
!!$          prev%next => token%next
!!$       else
!!$          my_dtsets => token%next
!!$       end if
!!$       deallocate(token)
!!$       nb_dtsets = max(nb_dtsets - 1, 0)
       end if
    end if
    if (AB_DBG) write(0,*) "AB module: Free done"
  end subroutine ab_dtset_free

  subroutine ab_dtset_get_ndtset(dtsetsId, value, errno)
    integer, intent(in) :: dtsetsId
    integer, intent(out) :: value
    integer, intent(out) :: errno

    type(dtsets_list), pointer :: token

    call ab_dtset_get_from_id(token, dtsetsId)
    if (associated(token)) then
       value = size(token%dtsets) - 1
       errno = AB_NO_ERROR
    else
       errno = AB_ERROR_DTSET_OBJ
    end if
  end subroutine ab_dtset_get_ndtset

  include "ab_dtset_f90_get.f90"
end module ab_dtset
