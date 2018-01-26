!************************************************************************
!*                                                                      *
!* Solve a first order system of DEs using the implicit Gear method     *
!* of order 4.                                                          *
!*                                                                      *
!* Programming language: ANSI C                                         *
!* Author:               Klaus Niederdrenk, 1.22.1996 (FORTRAN 77)      *
!* Adaptation:           Juergen Dietel, Computing Center, RWTH Aachen  *
!* Source:               FORTRAN 77 source code                         *
!* Date:                 2.26.1996                                      *
!*                                                                      *
!*                       F90 Release By J-P Moreau, Paris.              *
!*                       (www.jpmoreau.fr)                              *
!* -------------------------------------------------------------------- *
!* REF.: "Numerical Algorithms with C, By Gisela Engeln-Muellges        *
!*        and Frank Uhlig, Springer-Verlag, 1996" [BIBLI 11].           *
!************************************************************************

!print a REAL square matrix with name (debug only)
Subroutine PrintMat(Name, n, mat)
  character*(*) Name
  integer n
  real*8 mat(0:n-1,0:n-1)
  integer i,j
  print *, Name
  print *,((mat(i,j),j=0,n-1),i=0,n-1)
  return
end

!Gear  method of 4th order for DESs of 1st order
Subroutine gear4(x,xend,bspn,n,y,epsabs,epsrel,h,fmax,aufrufe,fehler)
  real*8 x,     &               !starting or end point ................
  xend                          !desired end point (> x) ..............
  integer bspn, &               !# of example .........................
  n                             !number of DEs ........................
  real*8 y(0:n-1),  &           !initial value or solution ............
  epsabs,       &               !absolute error bound .................
  epsrel,       &               !relative error bound .................
  h                             !starting or final step size ..........
  integer fmax, &               !maximal number of calls of dgl() .....
  aufrufe,      &               !actual number of calls of dgl() ......
  fehler                        !error code ...........................
!************************************************************************
!* Compute the value of the solution at xend, starting with the IVP.    *
!* We use the implicit method of Gear of 4th order with step size       *
!* control which is especially suited for stiff DEs.                    *
!* The local step size control insures that the two error bounds are met*
!* The number of function calls of the right hand side is limited by    *
!* fmax. This function can be used inside a loop to find solutions at   *
!* a specified point to arbitrary accuracy.                             *
!*                                                                      *
!* Input parameters:                                                    *
!* =================                                                    *
!* x        initial value x0                                            *
!* xend     final value for the integration (xend > x0)                 *
!* bspn     # example                                                   *
!* n        number of DEs                                               *
!* dgl      Function to compute the right hand side f(x0,y0) for the    *
!*          system of DEs                                               *
!* y        [0..n-1] solution vector y0 of the system of DEs at x0      *
!* epsabs   absolute error bound (>= 0); if zero, we only check for the *
!*          relative error.                                             *
!* epsrel   relative error bound (>= 0); if zero, we only check for the *
!*          absolute error.                                             *
!* h        given step size for first step                              *
!* fmax     maximal number of calls of the right hand side of the system*
!*                                                                      *
!* Output parameters:                                                   *
!* =================                                                    *
!* x        final x-value of the integration; normally equal to xend    *
!* h        final step size; keep for further calls                     *
!* y        [0..n-1]-vector, the solution of the system of DEs at x     *
!* aufrufe  counter of calls of dgl()                                   *
!*                                                                      *
!* erroe code (fehler):                                                 *
!* ====================                                                 *
!*   =   0:  all ok                                                     *
!*   =   1:  Both error bounds too small                                *
!*   =   2:  xend <= x0                                                 *
!*   =   3:  Step size h <= 0                                           *
!*   =   4:  n <= 0                                                     *
!*   =   5:  # calls > fmax; we have not reached the desired xend;      *
!*           repeat function call with actual values of x, y, h.        *
!*   =   6:  Jacobi matrix is singular; x, y, h are the last values     *
!*           that could be computed with accuracy                       *
!*   =   7:  ran out of memory (not used here)                          *
!*   =   8:  error when calling gauss() for the second time;            *
!*           should not occur.                                          *
!*                                                                      *
!************************************************************************
  real*8, parameter :: ONE = 1.d0

  real*8 eps1,  &         !MACH_EPS**0.75; used instead of MACH_EPS
                          !to check whether we have reached xend in
                          !order to avoid minuscule step sizes
  eps2,         &         !100 * MACH_EPS; for comparison with zero
  hs                      !optimal step size for Jacobi matrix
                          !approximation

  real*8, pointer :: hilf(:)            !pointer to a (0..n-1)-vector
  real*8, pointer :: zj(:,:), zjp1(:,:) !pointer to a (0:4,0:n-1)-matrix 
  real*8, pointer :: f(:),ykp1(:)       !pointers to (0..n-1)-vectors
  real*8, pointer :: fs(:,:) ,fsg(:,:)  !pointer to (0..n-1,0..n-1)-matrices
  real*8, pointer :: con(:)             !pointers to a (0..n-1)-vector
  integer,pointer :: perm(:)            !pointer to a (0..n-1) permutation
                                        !vector for Gauss elimination

  real*8 sg,       &                    !sign of xend
  xe                                    !|xend| - eps2, carrying the sign of xend                             }
  integer amEnde                        !Flag, indicating that we shall reach
                                        !xend with the current step
  real*8 ymax,     &                    !Maximum norm of the current
                                        !approximate value of y
  dummy,           &                    !aux storage for h
  xka,             &
  xke, hka, hk1, diff, eps, q, halt, quot1, quot2, quot3, quot4

  integer done, nochmal                 !Boolean
  integer aufrufe_awp,  &
  signdet,              &               !sign of determinant in Gauss
  iter, i, k

  real*8, pointer ::dum(:), dum1(:)     !dummy vectors for gauss
  real*8  MACH_EPS                      !machine smallest real number
  real*8  dist_max, norm_max, Max, Min

  !auxiliary vectors
  real*8 yhilf(0:n-1),k1(0:n-1),k2(0:n-1),k3(0:n-1),k4(0:n-1),k5(0:n-1),k6(0:n-1)

! ------------------------- Initialize ------------------------------
  dummy =0.D0; done=0
  MACH_EPS = 1.2D-16                    !for IBM PC
  eps1  = MACH_EPS**0.75
  eps2  = 100.D0 * MACH_EPS
  hs    = 10.D0 * DSQRT(MACH_EPS)

! ----------------- Allocate matrices & vectors -------------------
  allocate(hilf(0:n-1),stat=ialloc)
  allocate(zj(0:4,0:n-1),stat=ialloc)
  allocate(zjp1(0:4,0:n-1),stat=ialloc)
  allocate(f(0:n-1),stat=ialloc)
  allocate(ykp1(0:n-1),stat=ialloc)
  allocate(fs(0:n-1,0:n-1),stat=ialloc)
  allocate(fsg(0:n-1,0:n-1),stat=ialloc)
  allocate(con(0:n-1),stat=ialloc)
  allocate(perm(0:n-1),stat=ialloc)
  allocate(dum(0:n-1),stat=ialloc)
  allocate(dum1(0:n-1),stat=ialloc)

  if (xend >= 0.D0) then
    sg = 1.D0
  else 
    sg = -1.D0
  end if
  xe       = (1.d0 - sg * eps2) * xend
  fehler   = 0
  aufrufe  = 1
  amEnde   = 0

! --------- check input parameters -------------------
  ymax = norm_max(y, n)
  if (epsabs <= eps2 * ymax.and.epsrel <= eps2) then
    fehler = 1
  else if (xe < x) then
    fehler = 2
  else if (h < eps2 * DABS(x)) then
    fehler = 3
  else if (n <= 0) then
    fehler = 4
  end if
  if (fehler>0) return 

! ------------ first integration step ---------------
  if (x + h > xe) then
    h      = xend - x
    dummy  = h
    amEnde = 1
  end if
  do i=0, n-1
    hilf(i)=y(i)
  end do
  xka = x;
  xke = xka;
  hka = 0.25d0 * h
  hk1 = hka

  do k = 1, 4
    xke = xke + hka

    call awp(xka, xke, bspn, n, hilf, epsabs, epsrel, hk1,6, fmax - aufrufe, aufrufe_awp, &
	fehler, yhilf, k1, k2, k3, k4, k5, k6)    
	
	aufrufe = aufrufe + aufrufe_awp

    if (fehler.ne.0) return
    do i=0, n-1 
	  zjp1(k,i) = hilf(i)
    end do
  end do
! Dbug M.O
  call dgl(bspn,n,x,y,f)
  aufrufe = aufrufe + 1

! ---------- Compute first Gear-Nordsieck approximation -------------------
  do i = 0, n-1
    zj(0,i) = y(i)
    zj(1,i) = h * f(i)
    zj(2,i) = ONE / 24.d0 * (35.d0 * y(i) - 104.d0 * zjp1(1,i) +  &
              114.d0 * zjp1(2,i) - 56.d0 * zjp1(3,i) + 11.d0 * zjp1(4,i))
    zj(3,i) = ONE / 12.d0 * (-5.d0 * y(i) + 18.d0  * zjp1(1,i) -  &
              24.d0 * zjp1(2,i) + 14.d0 * zjp1(3,i) - 3.d0 * zjp1(4,i))
    zj(4,i) = ONE / 24.d0 * (y(i) - 4.d0 * zjp1(1,i) + 6.d0 * zjp1(2,i) -  &
              4.d0 * zjp1(3,i) + zjp1(4,i))
  end do

! ------------------------ adjust step size --------------------------
  do while (done.eq.0)
  
    ! --- use Newton method for an implicit approximation ---

    do i=0, n-1
      ykp1(i) = zj(0,i) + zj(1,i) + zj(2,i) + zj(3,i) + zj(4,i)
    end do
! Dbug M.O
    call dgl(bspn,n,x+h,ykp1,f)

    do k=0, n-1
      !copy vector ykp1 in hilf
      do i=0, n-1
	    hilf(i)=ykp1(i)
      end do
      hilf(k) = hilf(k) - hs
	! Dbug M.O
      call dgl(bspn,n,x+h,hilf,dum);
      do i=0, n-1 
	    fs(k,i)=dum(i)
      end do 
      do i=0, n-1
        fs(k,i) = -h * 0.48D0 * (f(i) - fs(k,i)) / hs
      end do
      fs(k,k) = fs(k,k) + ONE
    end do

    !update number of calls to dgl
    aufrufe = aufrufe + n + 1

    do i=0, n-1
      con(i) = ykp1(i) - 0.48D0 * (zj(1,i) + 2.D0 * zj(2,i) +  &
               3.D0 * zj(3,i) + 4.D0 * zj(4,i))
      do k=0, n-1
        fsg(k,i) = fs(i,k)
      end do
    end do

    call gauss(1, n, fsg, fsg, perm, dum, dum1, signdet,fehler)

    if (fehler>0) then  !error in gauss ?
      fehler = 6
      return
    end if

    do iter=1, 3
      do i=0, n-1
        hilf(i) = - ykp1(i)
        do k=0, n-1
          hilf(i) = hilf(i) + fs(k,i) * ykp1(k)
        end do
        hilf(i) = h * 0.48d0 * f(i) + hilf(i) + con(i)
      end do

      call gauss(2, n, fsg, fsg, perm, hilf, ykp1, signdet,fehler)

      if (fehler > 0) then
        fehler=8
        return
      end if
! Dbug M.O
      call dgl(bspn,n,x+h,ykp1,f)

    end do
    !update number of calls to dgl
    aufrufe = aufrufe + 3

!   ---- Compute corresponding Gear-Nordsieck approximation ----

    do i=0, n-1
      hilf(i) = h * f(i) - zj(1,i) - 2.d0 * zj(2,i) -  &
                3.d0 * zj(3,i) - 4.d0 * zj(4,i)
    end do

    do i=0, n-1
      zjp1(0,i) = ykp1(i)
      zjp1(1,i) = h * f(i)
      zjp1(2,i) = zj(2,i) + 3.d0 * zj(3,i) + 6.d0 * zj(4,i) + 0.7d0 * hilf(i)
      zjp1(3,i) = zj(3,i) + 4.d0 * zj(4,i) + 0.2d0 * hilf(i)
      zjp1(4,i) = zj(4,i) + 0.02d0 * hilf(i)
    end do

!   --- decide whether to accept last step ---

!   copy vector zjp1(4) in hilf and zj(4) in con 
    do i=0, n-1
      hilf(i) = zjp1(4,i)
      con(i)  = zj(4,i)
    end do

    diff = dist_max(hilf, con, n)
    ymax = norm_max(ykp1, n)
    eps  = (epsabs + epsrel * ymax) / 6.d0
    q    = DSQRT(DSQRT(eps / diff)) / 1.2d0

    if (diff < eps) then

      ! --- accept last step; prepare for next one ---

      x = x + h
      !copy vector ykp1 in y
      do i=0, n-1 
	    y(i) = ykp1(i)
      end do

!     stop integration, if interval end xend has been reached
!     or if there has been too many function dgl calls.

      nochmal = 0
      do while (nochmal.eq.0)
        if (amEnde.ne.0) then
          h = dummy
          return
        else if (aufrufe > fmax) then
          fehler = 5
          return
        end if

        ! --- adjust step size for next step ---
        halt = h
        h    = MIN(q, 2.d0) * h

        if (x + h >= xe) then
          dummy  = h
          h      = xend - x
          amEnde = 1

          ! --- close enough to xend => stop integration ---
          if (h < eps1 * DABS(xend))  nochmal = 1
        end if
		if (nochmal==0) goto 10
      end do

!     ------ compute Gear-Nordsieck approximation -----
!     ------ for the next step                    -----
10    quot1 = h / halt
      quot2 = quot1**2
      quot3 = quot2 * quot1
      quot4 = quot3 * quot1
      do i=0, n-1
        zj(0,i) = zjp1(0,i)
        zj(1,i) = quot1 * zjp1(1,i)
        zj(2,i) = quot2 * zjp1(2,i)
        zj(3,i) = quot3 * zjp1(3,i)
        zj(4,i) = quot4 * zjp1(4,i)
	  end do
    else
!     ------ repeat last step with smaller step size;  -----
!     -------- adjust Gear-Nordsieck approximation ---------
      halt  = h
      h     = Max(0.5d0, q) * h
      quot1 = h / halt
      quot2 = quot1**2
      quot3 = quot2 * quot1
      quot4 = quot3 * quot1
      do i=0, n-1
        zj(1,i) = quot1 * zj(1,i)
        zj(2,i) = quot2 * zj(2,i)
        zj(3,i) = quot3 * zj(3,i)
        zj(4,i) = quot4 * zj(4,i)
      end do
      amEnde = 0
    end if
  end do !while done=0
return
End !gear4

! -------------------------- END  gear.f90 -------------------------
