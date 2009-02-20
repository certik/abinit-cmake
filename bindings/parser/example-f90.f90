program toto

  use ab_dtset

  integer :: id, errno

  integer :: natom, ndtset, ixc, idtset, brvltt
  double precision, allocatable :: coord(:)
  integer :: dims(7), ndims
  character(len = 256) :: filename

  call getarg(1, filename)

  call ab_dtset_new(id, filename, len(filename))

  call ab_dtset_get_ndtset(id, ndtset, errno)
  if (errno /= AB_NO_ERROR) goto 1000
  write(*,*) ndtset
  do idtset = 1, ndtset, 1
     write(*,*) "Dataset", idtset
     call ab_dtset_get_integer(id, natom, ab_dtset_natom, ndtset, errno)
     if (errno /= AB_NO_ERROR) goto 1000
     write(*,*) "natom", natom
     call ab_dtset_get_shape(id, dims, ndims, ab_dtset_xred_orig, ndtset, errno)
     if (errno /= AB_NO_ERROR) goto 1000
     write(*,*) "dims of xred_orig", ndims, " - ", dims(1:ndims)
     allocate(coord(product(dims(1:ndims))))
     call ab_dtset_get_real_array(id, coord, size(coord), ab_dtset_xred_orig, ndtset, errno)
     if (errno /= AB_NO_ERROR) goto 1000
     write(*,*) "xred_orig", coord
     call ab_dtset_get_real(id, coord(1), ab_dtset_wvl_hgrid, ndtset, errno)
     if (errno /= AB_NO_ERROR) goto 1000
     write(*,*) "wvl_hgrid", coord(1)
     deallocate(coord)
     call ab_dtset_get_integer(id, ixc, ab_dtset_ixc, ndtset, errno)
     if (errno /= AB_NO_ERROR) goto 1000
     write(*,*) "ixc", ixc
     call ab_dtset_get_integer(id, brvltt, ab_dtset_brvltt, ndtset, errno)
     if (errno /= AB_NO_ERROR) goto 1000
     write(*,*) "brvltt", brvltt
  end do

  1000 continue
  call ab_dtset_free(id)
  if (allocated(coord)) deallocate(coord)
  if (errno /= AB_NO_ERROR) then
     write(0, *) "Error!"
  end if

end program toto
