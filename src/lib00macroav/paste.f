! Copyright (C) 1999-2003 (P. Ordejon, J. Junquera)
! Modified by TD (2008)
! This file is distributed under the terms of the
! GNU General Public License, see ~abinit/COPYING
! or http://www.gnu.org/copyleft/gpl.txt .

      CHARACTER(LEN=26) FUNCTION PASTE( STR1, STR2 )

      IMPLICIT NONE

C CONCATENATES THE STRINGS STR1 AND STR2 REMOVING BLANKS IN BETWEEN

      CHARACTER(LEN=*) :: STR1, STR2
      INTEGER :: L
      DO L = LEN( STR1 ), 1, -1
         IF (STR1(L:L) .NE. ' ') EXIT
      END DO
      PASTE = STR1(1:L)//STR2
      END FUNCTION PASTE


      CHARACTER(LEN=26) FUNCTION PASTEB( STR1, STR2 )

! CONCATENATES THE STRINGS STR1 AND STR2 LEAVING ONLY ONE BLANK IN BETWEEN

      IMPLICIT NONE

      CHARACTER(LEN=*) :: STR1, STR2 
      CHARACTER(LEN=1), PARAMETER :: BLANK=' '
      INTEGER :: L
      DO L = LEN( STR1 ), 1, -1
         IF (STR1(L:L) .NE. ' ') EXIT
      END DO
      PASTEB = STR1(1:L)//BLANK
      PASTEB = PASTEB(1:L+1)//STR2
      END FUNCTION PASTEB

