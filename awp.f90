!************************************************************************
!*                                                                      *
!* Solve an ordinary system of first order differential equations using *
!* -------------------------------------------------------------------- *
!* automatic step size control                                          *
!* ----------------------------                                         *
!*                                                                      *
!* Programming language: ANSI C                                         *
!* Author:               Klaus Niederdrenk (FORTRAN)                    *
!* Adaptation:           Juergen Dietel, Computer Center, RWTH Aachen   *
!* Source:               existing C, Pascal, QuickBASIC and FORTRAN     *
!*                       codes                                          *
!* Date:                 6.2.1992, 10.2.1995                            *
!*                                                                      *
!*                       F90 Release By J-P Moreau, Paris.              *
!*                       (www.jpmoreau.fr)                              *
!* -------------------------------------------------------------------- *
!* REF.: "Numerical Algorithms with C, By Gisela Engeln-Muellges        *
!*        and Frank Uhlig, Springer-Verlag, 1996" [BIBLI 11].           *
!************************************************************************

!debug only
Subroutine PrintVec(name, n, V)
character*(*) name
integer n,  i
real*8 V(0:n-1)
  print *, name
  write(*,10) (V(i), i=0,n-1)
  print *,' '
  return
10 format(10E12.6)
End


! Maximum norm of a difference vector ........
real*8 Function dist_max(vector1,vector2, n)
       real*8  vector1(0:n-1),  &
               vector2(0:n-1)
       integer n
!************************************************************************
!* Compute the maximum norm of the difference of two [0..n-1] vectors   *
!*                                                                      *
!* global name used:                                                    *
!* ================                                                     *
!*   None                                                               *
!************************************************************************
  real*8  abstand,   &     !reference value for computation of distance
          hilf             !distance of two vector elements
  integer i

  abstand = 0.d0
  do i=n-1, 0, -1
    hilf = DABS(vector1(i)-vector2(i))
    if (hilf > abstand)  abstand = hilf
  end do
  dist_max = abstand
  return
End


! Runge-Kutta embedding formulas of 2nd, 3rd degree ....................
Subroutine ruku23(x,y,bspn,n,h,y2,y3,yhilf,k1,k2,k3)
   real*8 x, y(0:n-1)
   integer bspn, n
   real*8 h, y2(0:n-1), y3(0:n-1)
! auxiliary vectors
   real*8 yhilf(0:n-1),k1(0:n-1),k2(0:n-1),k3(0:n-1)
!************************************************************************
!* Compute 2nd and 3rd order approximates y2, y3 at x + h starting with *
!* a solution y at x by using Runge-Kutta embedding formulas on the     *
!* first order system of n differential equations   y' = f(x,y) , as    *
!* supplied by  dgl().                                                  *
!*                                                                      *
!* Input parameters:                                                    *
!* =================                                                    *
!* x    x-value of left end point                                       *
!* y    y-values at x                                                   *
!* bspn # example                                                       *
!* n    number of differential equations                                *
!* h    step size                                                       *
!*                                                                      *
!* yhilf,k1,k2,k3: auxiliary vectors defined in awp.                    *
!*                                                                      *
!* Output parameters:                                                   *
!* ==================                                                   *
!* y2   2nd order approximation for y at x + h                          *
!* y3   3rd order approximation for y at x + h                          *
!*                                                                      *
!* External function:                                                   *
!* =================                                                    *
!* dgl  function that evaluates the right hand side of the system       *
!*      y' = f(x,y)                                                     *
!*      (See file t_dgls.f90).                                          * 
!************************************************************************
  integer i               !loop variable

  call dgl(bspn,n,x,y,k1)
  do i=0, n-1
    yhilf(i) = y(i) + h * k1(i)
  end do
  call dgl(bspn,n,x+h,yhilf,k2)
  do i=0, n-1
    yhilf(i) = y(i) + 0.25D0 * h * (k1(i) + k2(i))
  end do
  call dgl(bspn,n,x+0.5*h,yhilf,k3)
  do i=0, n-1
    y2(i) = y(i) + 0.5D0 * h * (k1(i) + k2(i))
    y3(i) = y(i) + h / 6.D0 * (k1(i) + k2(i) + 4.D0 * k3(i))
  end do
  return
End


! England's Einbettungs formulas of 4th and 5th degree ................
Subroutine engl45(x,y,bspn,n,h,y4,y5,yhilf,k1,k2,k3,k4,k5,k6)
  real*8    x                   !starting point of integration ........
  real*8    y(0:n-1)            !initial value at x ...................
  integer   bspn,       &       !# example
            n                   !number of differential equations .....
  real*8    h                   !step size ............................
  real*8   y4(0:n-1),   &       !4th order approximation for y at x + h
           y5(0:n-1)            !5th order approximation for y at x + h
  !auxiliary vectors
  real*8 yhilf(0:n-1),k1(0:n-1),k2(0:n-1),k3(0:n-1),k4(0:n-1),k5(0:n-1),k6(0:n-1)
!************************************************************************
!* Compute 4th and 5th order approximates y4, y5 at x + h starting with *
!* a solution y at x by using the England embedding formulas on the     *
!* first order system of n differential equations   y' = f(x,y) , as    *
!* supplied by  dgl().                                                  *
!*                                                                      *
!* Input parameters:                                                    *
!* =================                                                    *
!* x    initial x-value                                                 *
!* y    y-values at x, type pVEC                                        *
!* n    number of differential equations                                *
!* h    step size                                                       *
!*                                                                      *
!* yhilf, k1..K6: auxiliary vectors                                     *
!*                                                                      *
!* Output parameters:                                                   *
!* ==================                                                   *
!* y4   4th order approximation for y at x + h (pVEC)                   *
!* y5   5th order approximation for y at x + h (pVEC)                   *
!*                                                                      *
!* External function:                                                   *
!* =================                                                    *
!* dgl  function that evaluates the right hand side of the system       *
!*      y' = f(x,y)                                                     *
!*      (See file t_dgls.f90).                                          * 
!************************************************************************
  integer i               !loop variable
  call dgl(bspn,n,x,y,k1)
  do i=0, n-1
    yhilf(i) = y(i) + 0.5d0 * h * k1(i)
  end do
  call dgl(bspn,n,x+0.5d0*h,yhilf,k2)
  do i=0, n-1
    yhilf(i) = y(i) + (0.25D0 * h * (k1(i) + k2(i)))
  end do
  call dgl(bspn,n,x+0.5d0*h,yhilf,k3)
  do i=0, n-1
    yhilf(i) = y(i) + h * (-k2(i) + 2.D0 * k3(i))
  end do
  call dgl(bspn,n,x+h,yhilf,k4)
  do i=0, n-1
    yhilf(i) = y(i) + h / 27.D0 * (7.D0 * k1(i) + 10.D0 * k2(i) + k4(i))
  end do
  call dgl(bspn,n,x+(2.d0/3.d0)*h,yhilf,k5)
  do i=0, n-1
    yhilf(i) = y(i) + h / 625.d0 * (28.d0 * k1(i) - 125.d0 * k2(i) + &
               546.d0 * k3(i) + 54.d0 * k4(i) - 378.d0 * k5(i))
  end do
  call dgl(bspn,n,x+h/5.d0,yhilf,k6)
  do i=0, n-1
    y4(i) = y(i) + h / 6.d0 * (k1(i) + 4.d0 * k3(i) + k4(i))
    y5(i) = y(i) + h / 336.d0 * (14.d0 * k1(i) + 35.d0 * k4(i) +  &
            162.d0 * k5(i) + 125.d0 * k6(i))
  end do
  return
End


! embedding formulas of Prince-Dormand of 4./5. order .........................
Subroutine prdo45(x,y,bspn,n,h,y4,y5,steif1,steifanz,steif2, &
                  yhilf,k1,k2,k3,k4,k5,k6,k7, g6,g7)         
  real*8    x,           &          !starting point of integration .....
            y(0:n-1)                !initial value at x ................
  integer   bspn,        &          !# example .........................
            n                       !number of DEs .....................
  real*8    h,           &          !step size .........................
            y4(0:n-1),   &          !solution of 4th order at x+h ......
            y5(0:n-1)               !solution of 5th order at x+h ......
  !auxiliary flags
  integer steif1,steifanz,steif2
  !auxiliary vectors
  real*8 yhilf(0:n-1)
  real*8 k1(0:n-1),k2(0:n-1),k3(0:n-1),k4(0:n-1),k5(0:n-1),k6(0:n-1),k7(0:n-1)
  real*8 g6(0:n-1),g7(0:n-1)
!************************************************************************
!* Compute 4th and 5th order approximates y4, y5 at x + h starting with *
!* a solution y at x by using the Prince-Dormand embedding formulas on  *
!* the first order system of n differential equations y' = f(x,y) , as  *
!* supplied by  dgl().                                                  *
!* Simultaneously we perform two tests for stiffness whose results are  *
!* stored in steif1 and steif2.                                         *
!*                                                                      *
!* Input parameters:                                                    *
!* =================                                                    *
!* x    initial x-value                                                 *
!* y    y-values at x (pVEC)                                            *
!* n    number of differential equations                                *
!* h    step size                                                       *
!*                                                                      *
!* yhilf, k1..k7,g6,g7: auxiliary vectors.                              *
!*                                                                      *
!* Output parameters:                                                   *
!* ==================                                                   *
!* y4   4th order approximation for y at x + h                          *
!* y5   5th order approximation for y at x + h                          *
!*                                                                      *
!* External function:                                                   *
!* =================                                                    *
!* dgl  function that evaluates the right hand side of the system       *
!*      y' = f(x,y)                                                     *
!*      (See file t_dgls.f90).                                          * 
!***********************************************************************}
  integer i            !loop variable
  integer steifa       !Flag which is set if the second test for stiffness
                       !Shampine und Hiebert) is positive; otherwise the
                       !flag is erased.
  real*8 dist_max

  call dgl(bspn,n,x,y,k1)
  do i=0, n-1
    yhilf(i) = y(i) + 0.2D0 * h * k1(i)
  end do
  call dgl(bspn,n,x+0.2D0*h,yhilf,k2)
  do i=0, n-1
    yhilf(i) = y(i) + 0.075D0 * h * (k1(i) + 3.D0 * k2(i))
  end do
  call dgl(bspn,n,x+0.3D0*h,yhilf,k3)
  do i=0, n-1
    yhilf(i) = y(i) + h / 45.D0 * (44.D0 * k1(i) - 168.D0 * k2(i) + 160.D0 * k3(i))
  end do
  call dgl(bspn,n,x+0.8D0*h,yhilf,k4)
  do i=0, n-1
    yhilf(i) = y(i) + h / 6561.D0 * (19372.D0 * k1(i) - 76080.D0 * k2(i)  &
               + 64448.D0 * k3(i) - 1908.D0 * k4(i))
  end do
  call dgl(bspn,n,x+(8.D0/9.D0)*h,yhilf,k5)
  do i=0, n-1
    g6(i) = y(i) + h / 167904.D0 * (477901.D0 * k1(i) - 1806240.D0 * k2(i)  &
            + 1495424.D0 * k3(i) + 46746.D0 * k4(i) - 45927.D0 * k5(i))
  end do
  call dgl(bspn,n,x+h,g6,k6)
  do i=0, n-1
    g7(i) = y(i) + h / 142464.D0 * (12985.D0 * k1(i) + 64000.D0 * k3(i)  &
            + 92750.D0 * k4(i) - 45927.D0 * k5(i) + 18656.D0 * k6(i))
  end do
  call dgl(bspn,n,x+h,g7,k7)
  do i=0, n-1
    y5(i) = g7(i)
    y4(i) = y(i) + h / 21369600.D0 * (1921409.D0 * k1(i) + 9690880.D0 * k3(i)  &
            + 13122270.D0 * k4(i)- 5802111.D0 * k5(i) + 1902912.D0 * k6(i)     &
            + 534240.D0 * k7(i))
  end do

! Test for stiffness via dominant eigenvalue

  if (dist_max(k7, k6, n) > 3.3D0 * dist_max(g7, g6, n))  steif1 = 1

! one step in steffness test of Shampine & Hiebert

  do i=0, n-1
    g6(i) = h * (2.2D0 * k2(i) + 0.13D0 * k4(i) + 0.144D0 * k5(i))
    g7(i) = h * (2.134D0 * k1(i) + 0.24D0 * k3(i) + 0.1D0 * k6(i))
  end do

  if (dist_max(g6, g7, n) < dist_max(y4, y5, n))  then
    steifa = 1
  else  
    steifa = 0
  end if

  if (steifa > 0) then
    steifanz = steifanz + 1
    if (steifanz >= 3) steif2 = 1
  else
    steifanz = 0
  end if
  return

End


! Find the maximum norm of a REAL vector ............................
real*8 Function norm_max(vektor, n)    
       real*8  vektor(0:n-1)               ! vector .................
       integer n                           ! length of vector .......
! ************************************************************************
! *          Return the maximum norm of a [0..n-1] vector v              *
! *                                                                      *
! ************************************************************************
    real*8 norm,  &                                         ! local max.
    betrag                                    ! magnitude of a component
    integer i
  
    norm=0.D0
    do i=0, n-1
      betrag = dabs(vektor(i))
      if (betrag > norm)  norm = betrag
    end do   
    norm_max = norm
	return
End

! 1st order DESs with automatic step size control ...........................
Subroutine awp(x,xend,bspn,n,y,epsabs,epsrel,h,methode,fmax,aufrufe,fehler, &
               yhilf, k1, k2, k3, k4, k5, k6)
   real*8 x,             &        !initial/final x value ..............
          xend                    !desired end point ..................
  integer bspn,          &        !# example 
          n                       !number of DEs ......................
   real*8 y(0:n-1),      &        !initial/final y value ..............
          epsabs,        &        !absolute error bound ...............
          epsrel,        &        !relative error bound ...............
          h                       !initial/final step size ............
  integer methode,       &        !desired method (3, 6, 7) ...........
          fmax,          &        !maximal # of calls of  dgl() .......
          aufrufe,       &        !actual # of calls of  dgl() ........
          fehler                  !error code .........................
  !auxiliary vectors
  real*8 yhilf(0:9),k1(0:9),k2(0:9),k3(0:9),k4(0:9),k5(0:9),k6(0:9)
!************************************************************************
!* Compute the solution y of a system of first order ordinary           *
!* differential equations       y' = f(x,y)   at xend from the given    *
!* initial data (x0, y0).                                               *
!* We use automatic step size control internally so that the error of   *
!* y (absolutely or relatively) lies within the given error bounds      *
!* epsabs and epsrel.                                                   *
!*                                                                      *
!* Input parameters:                                                    *
!* =================                                                    *
!* x        initial value for x                                         *
!* y        initial values for y(0:n-1)                                 *
!* bspn     # example                                                   *
!* n        number of differential equations                            *
!* dgl      function that evaluates the right hand side of the system   *
!*          y' = f(x,y)  (see t_dgls) (Removed from list of parameters) *
!* xend     end of integration; xend > x0                               *
!* h        initial step size                                           *
!* epsabs   absolute error bound; >= 0; if = 0 we only check the        *
!*          relative error.                                             *
!* epsrel   relative error bound; >= 0; if = 0 we check only the        *
!*          absolute eror.                                              *
!* fmax     max number of evaluations of right hand side in dgl()       *
!* methode  chooses the method                                          *
!*          = 3: Runge-Kutta method of 2nd/3rd order                    *
!*          = 6: England formula of 4th/5th order                       *
!*          = 7: Formula of Prince-Dormand of 4th/5th order             *
!*                                                                      *
!* Output parameters:                                                   *
!* ==================                                                   *
!* x        final x-value for iteration. If fehler = 0  we usually have *
!*            x = xend.                                                 *
!* y        final y-values for the solution at x                        *
!* h        final step size used; leave for subsequent calls            *
!* aufrufe  actual number of calls of dgl()                             *
!*                                                                      *
!* Return value (fehler):                                               *
!* =====================                                                *
!* = 0: all ok                                                          *
!* = 1: both error bounds chosen too small for the given mach. constant *
!* = 2: xend <= x0                                                      *
!* = 3: h <= 0                                                          *
!* = 4: n <= 0                                                          *
!* = 5: more right hand side calls than allowed: aufrufe > fmax,        *
!*      x and h contain the current values when stop occured.           *
!* = 6: improper input for embedding formula                            *
!* = 7: lack of available memory (not used here)                        *
!* = 8: Computations completed, but the Prince Dormand formula stiff-   *
!*      ness test indicates possible stiffness.                         *
!* = 9: Computations completed, but both Prince Dormand formula stiff-  *
!*      ness tests indicate possible stiffness. Use method for stiff    *
!*      systems instead !                                               *
!* =10: aufrufe > fmax, see error code 5; AND the Prince Dormand formula*
!*      indicates stiffness; retry using a stiff DE solver !            *
!*                                                                      *
!************************************************************************
  real*8, parameter :: MACH_EPS=1.2d-16
  real*8, parameter :: MACH_2 = 100.d0 * MACH_EPS
  real*8   &    
  xend_h,  &       !|xend| - MACH_2, carrying same sign as xend
  ymax,    &       !Maximum norm of newest approximation of max
                   !order.
  hhilf,   &       !aux storage for the latest value of h
                   !produced by step size control. It is saved
                   !here in order to avoid to return a `h' that
                   !resulted from an arbitrary reduction at the
                   !end of the interval.
  diff,    &       !distance of the two approximations from the
                   !embedding formula.
  s                !indicates acceptance level for results from
                   !embeding formula.

  !approximate solution of low order
  real*8, pointer :: y_bad(:)
  !ditto of high order
  real*8, pointer ::y_good(:)
   
  real*8 mach_1    !machine constant dependent variable which
                   !avoids using too little steps near xend.
  integer i,  &    !Loop variable
  amEnde,     &    !flag that shows if the end of the interval
                   !can be reached with the actual step size.
  fertig           !flag indicating end of iterations.

  integer    &
  steif1,    &     !Flag, that is set in prdo45() if its
                   !stiffness test (dominant eigenvalue)
                   !indicates so. Otherwise no changes.
  steifanz,  &     !counter for number of successive successes
                   !of stiffness test of Shampine and Hiebert in
                   !prdo45().
  steif2,    &     !Flag, set in prdo45(), when the stiffness
                   !test of  Shampine and Hiebert wa successful
                   !three times in a row; otherwise no changes.
  ialloc           !allocation status (must be zero).
   
  !additional auxiliary vectors for prdo45 method (not used here)
  real*8 k7(0:n-1), g6(0:n-1), g7(0:n-1)

  real*8 dist_max, norm_max, Min, Max   !external functions
 
  !allocate auxiliary vectors
  allocate(y_bad(0:n-1), stat=ialloc)
  allocate(y_good(0:n-1), stat=ialloc)

  if (ialloc.ne.0) then
    print *,' awp: memory full !'!
    return
  end if

  !initialize some variables
  fehler   = 0                           
  mach_1   = MACH_EPS**0.75d0
  amEnde   = 0
  fertig   = 0
  steif1   = 0
  steif2   = 0
  steifanz = 0
  aufrufe  = 1
  ymax     = norm_max(y, n)

  if (xend >= 0.d0) then 
    xend_h = xend * (1.d0 - MACH_2)
  else 
    xend_h = xend * (1.d0 + MACH_2)
  end if

! ----------------------- check inputs ----------------------
  if (epsabs <= MACH_2 * ymax.and.epsrel <= MACH_2) then
    fehler=1
    return
  end if
  if (xend_h < x) then
    fehler= 2
    return
  end if
  if (h < MACH_2 * DABS(x)) then
    fehler=3
    return
  end if
  if (n <= 0) then
    fehler=4
    return
  end if
  if (methode.ne.3.and.methode.ne.6.and.methode.ne.7) then
    fehler=6
    return
  end if
 
! **********************************************************************
! *                                                                    *
! *                       I t e r a t i o n s                          *
! *                                                                    *
! **********************************************************************
  if (x + h > xend_h) then            
                                    !almost at end point ?
    hhilf  = h                      !A shortened step might be
    h      = xend - x               !enough.
    amEnde = 1
  end if

  do while(fertig.eq.0)             !solve DE system by integrating from
                                    !x0 to xend by suitable steps.
    Select Case(methode)
      Case(3) 
	    call ruku23(x, y, bspn,n, h, y_bad, y_good, yhilf,k1,k2,k3)
      Case(6)
	    call engl45(x, y, bspn, n, h, y_bad, y_good, yhilf,k1,k2,k3,k4,k5,k6)
      Case(7)
	    call prdo45(x, y, bspn, n, h, y_bad, y_good, steif1,steifanz,steif2, &
		            yhilf,k1,k2,k3,k4,k5,k6,k7, g6,g7)
    End Select

    aufrufe = aufrufe + methode

	diff = dist_max(y_bad, y_good, n)

    if (diff < MACH_2) then         !compute s
      s = 2.d0
    else
      ymax = norm_max(y_good, n)
      s    = DSQRT(h * (epsabs + epsrel * ymax) / diff)
      if (methode.ne.3)   s = DSQRT(s)
    end if

    if (s > 1.d0) then              !integration acceptable ?
      do i=0, n-1                   !accept highest order solution
        y(i) = y_good(i)            !move x
      end do
      
	  x = x + h

      if (amEnde.ne.0) then         !at end of interval ?
        fertig = 1                  !stop iteration
        if (methode.eq.7) then
          if (steif1>0.or.steif2>0)  fehler = 8
          if (steif1>0.and.steif2>0) fehler = 9
        end if
      else if (aufrufe > fmax) then !too many calls of   dgl() ?
        hhilf  = h                  !save actual step size
        fehler = 5                  !report error and stop
        fertig = 1
        if (methode.eq.7.and.(steif1>0.or.steif2>0))  fehler = 10
      else                          !Integration was successful
                                    !not at the interval end ?
        h = h * Min(2.d0, 0.98d0*s) !increase step size for next
                                    !step properly, at most by
                                    !factor two. Value `0.98*s' is
                                    !used in order to avoid that
                                    !the theoretical value s is
                                    !exceeded by accidental
                                    !rounding errors.
        if (x + h > xend_h) then    !nearly reached xend ?
          hhilf  = h                !=> One further step with
          h      = xend - x         !reduced step size might be
          amEnde = 1                !enough.
          if (h < mach_1 * DABS(xend)) then  !very close to xend ?
            fertig = 1                       !finish iteration.
          end if
        end if
      end if
    else                            !step unsuccessful ?
                                    !before repeating this step
      h = h * Max(0.5d0, 0.98d0*s)  !reduce step size properly, at
                                    !most by factor 1/2 (for factor
      amEnde = 0                    !0.98: see above).
    end if

  end do
 
  h = hhilf           !return the latest step size computed by step
                      !size control and  error code to the caller.

  return
End

! -------------------------- END  uawp.f90 ----------------------------
