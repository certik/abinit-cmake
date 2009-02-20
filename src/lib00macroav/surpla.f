
c Copyright (C) 1999-2003 (P. Ordejon, J. Junquera)
c This file is distributed under the terms of the
c GNU General Public License, see ~abinit/COPYING
c or http://www.gnu.org/copyleft/gpl.txt .

      DOUBLE PRECISION FUNCTION SURPLA( C )

C  CALCULATES THE SRFACE OF THE UNIT CELL NORMAL TO THE INTERFACE

      IMPLICIT NONE

      DOUBLE PRECISION C(3,3)
      SURPLA = ( C(2,1)*C(3,2) - C(3,1)*C(2,2) ) **2 +
     .         ( C(3,1)*C(1,2) - C(1,1)*C(3,2) ) **2 +
     .         ( C(1,1)*C(2,2) - C(2,1)*C(1,2) ) **2
      SURPLA = SQRT( ABS( SURPLA ) )
      END FUNCTION SURPLA
