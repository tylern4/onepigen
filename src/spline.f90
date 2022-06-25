!------------------------------------------------------------------------
!EV The subroutines SPLINE and SPLINT can be used for interpolating a
!   function of one independent variable.
!------------------------------------------------------------------------

SUBROUTINE spline(x, y, n, yp1, ypn, y2)
    !
    ! Given arrays y(1:n) and x(1:n) containing a tabulated function and its
    ! independent variable, respectively, with x(1) < x(2) < ... < x(n), and given
    ! values yp1 and ypn for the 1st derivative of the interpolating function at
    ! the points 1 and n, respectively, this routine returns an array y2(1:n) of
    ! length n which contains the 2nd derivative of the interpolating function at
    ! the tabulated points x(i).  If yp1 and/or ypn are equal to 10^30 or larger,
    ! the routine is signaled to set the corresponding boundary condition for a
    ! natural spline, with zero 2nd derivative on that boundary.  From "Numerical
    ! Recipes", Ch. 3.3, p. 109.
    IMPLICIT NONE
    ! Passed variables:
    INTEGER n
    REAL x(n), y(n), yp1, ypn, y2(n)
    ! Local variables:
    INTEGER i, j, NMAX
    PARAMETER (NMAX = 500)        ! Maximum expected value of n
    REAL p, qn, sig, un, u(NMAX)
    !
    if (yp1.gt..99e30) then
        ! Natural lower boundary condition:
        y2(1) = 0.e0
        u(1) = 0.e0
    else
        ! Specified 1st derivative yp1:
        y2(1) = -.5e0
        u(1) = (3.e0 / (x(2) - x(1))) * ((y(2) - y(1)) / (x(2) - x(1)) - yp1)
    endif
    ! We need at least N+1 points to construct an interpolating polynomial
    ! of degree N:
    if (n.lt.4) stop ' Too few points for cubic spline'
    !
    ! Decomposition loop of the tridiagonal algorithm.
    ! y2,u: temporary storage of the decomposed factors.
    do i = 2, n - 1
        sig = (x(i) - x(i - 1)) / (x(i + 1) - x(i - 1))
        p = sig * y2(i - 1) + 2.e0
        y2(i) = (sig - 1.e0) / p
        u(i) = (6.e0 * ((y(i + 1) - y(i)) / (x(i + 1) - x(i)) - (y(i) - y(i - 1))&
                / (x(i) - x(i - 1))) / (x(i + 1) - x(i - 1)) - sig * u(i - 1)) / p
    enddo
    if (ypn.gt..99e30) then
        ! Natural upper boundary condition:
        qn = 0.e0
        un = 0.e0
    else
        ! Specified 1st derivative ypn:
        qn = .5e0
        un = (3.e0 / (x(n) - x(n - 1))) * (ypn - (y(n) - y(n - 1)) / (x(n) - x(n - 1)))
    endif
    y2(n) = (un - qn * u(n - 1)) / (qn * y2(n - 1) + 1.e0)
    !
    ! Backsubstitution loop of the tridiagonal algorithm.
    do j = n - 1, 1, -1
        y2(j) = y2(j) * y2(j + 1) + u(j)
    enddo
    !
    return
END