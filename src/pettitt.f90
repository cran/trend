SUBROUTINE pettitt (n, x, pval, tau, Kt)
!    Copyright (C) 2015, 2016  Thorsten Pohlert
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!    This subroutine computes the Ui statistics, the Kt = max(abs(Ui),
!    the approximate p-value and the location of the Pettitt-Test of 
!    change-point detection
!
  IMPLICIT NONE
  INTEGER, INTENT(in)          ::   n
  INTEGER                      ::   i, j, k, Kt, tau
  INTEGER, DIMENSION(n)        ::   U, Ustar
  REAL                         ::   pval
  REAL,DIMENSION(n), INTENT(in)::   x
  INTEGER                      :: sgn

  EXTERNAL sgn

  ! Initalize
  DO i = 1, n
     U(i) = 0
  END DO

  ! Get U
  DO i = 1, n-1
     DO k = 1, i
        DO j = i+1,n
           U(i) = sgn(x(k) - x(j)) + U(i)
        END DO
     END DO
  END DO

  Ustar = ABS(U)
  ! Get statistics
  Kt = MAXVAL(Ustar)
  ! Get tau
  tau = MAXLOC(Ustar,1)
  ! Get probability
  pval = 2.0 * EXP(( -6.0 * REAL(Kt)**2.0) /  &
             (REAL(n)**3.0 + REAL(n)**2.0))

  RETURN
END SUBROUTINE pettitt
 
INTEGER FUNCTION sgn(a)
! This function is the sgn function and returns +1 for positive 
! reals, -1 for negative reals and else 0.
  REAL, INTENT(in):: a 
  IF (a > 0.0) THEN
     sgn = 1
  ELSE IF (a < 0.0) THEN
     sgn = -1
  ELSE
     sgn = 0
  END IF
  RETURN
END FUNCTION sgn
