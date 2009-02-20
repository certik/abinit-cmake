
      DOUBLE PRECISION FUNCTION VOLCEL( C )

c Copyright (C) 1999-2003 (P. Ordejon, J. Junquera)
c This file is distributed under the terms of the
c GNU General Public License, see ~abinit/COPYING
c or http://www.gnu.org/copyleft/gpl.txt .

C  CALCULATES THE VOLUME OF THE UNIT CELL

      IMPLICIT NONE

      DOUBLE PRECISION C(3,3)
      VOLCEL = ( C(2,1)*C(3,2) - C(3,1)*C(2,2) ) * C(1,3) +
     .         ( C(3,1)*C(1,2) - C(1,1)*C(3,2) ) * C(2,3) +
     .         ( C(1,1)*C(2,2) - C(2,1)*C(1,2) ) * C(3,3)
      VOLCEL = ABS( VOLCEL )
      END FUNCTION VOLCEL
