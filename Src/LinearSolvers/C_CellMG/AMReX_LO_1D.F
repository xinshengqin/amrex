#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include "AMReX_LO_BCTYPES.H"
#include "AMReX_LO_F.H"
#include "AMReX_ArrayLim.H"

c-----------------------------------------------------------------------
      subroutine FORT_HARMONIC_AVERAGEEC (
     $     c, DIMS(c),
     $     f, DIMS(f),
     $     lo, hi, nc,
     $     cdir
     $     )
c
      integer nc
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer cdir
      integer DIMDEC(f)
      REAL_T f(DIMV(f),nc)
      integer DIMDEC(c)
      REAL_T c(DIMV(c),nc)
c
      integer i,n
c
      do n = 1, nc
         do i = lo(1), hi(1)
            c(i,n) =  f(2*i,n)
         end do
      end do
c
      end
c-----------------------------------------------------------------------
      subroutine FORT_AVERAGEEC (
     $     c, DIMS(c),
     $     f, DIMS(f),
     $     lo, hi, nc,
     $     cdir
     $     )
c
      integer nc
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer cdir
      integer DIMDEC(f)
      REAL_T f(DIMV(f),nc)
      integer DIMDEC(c)
      REAL_T c(DIMV(c),nc)
c
      integer i,n
c
      do n = 1, nc
         do i = lo(1), hi(1)
            c(i,n) = f(2*i,n)
         end do
      end do
c
      end
c-----------------------------------------------------------------------
      subroutine FORT_AVERAGECC (
     $     c, DIMS(c),
     $     f, DIMS(f),
     $     lo, hi, nc
     $     )
c
      integer nc
      integer DIMDEC(f)
      integer DIMDEC(c)
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      REAL_T f(DIMV(f),nc)
      REAL_T c(DIMV(c),nc)
c
      integer i,n
c
      do n = 1, nc
         do i = lo(1), hi(1)
            c(i,n) =  (f(2*i+1,n) + f(2*i,n))*half
         end do
      end do
c
      end
c-----------------------------------------------------------------------
      subroutine FORT_APPLYBC (
     $     flagden, flagbc, maxorder,
     $     phi,   DIMS(phi),
     $     cdir, bct, bcl,
     $     bcval, DIMS(bcval),
     $     mask,  DIMS(mask),
     $     den,   DIMS(den),
     $     lo, hi, nc,
     $     h
     $     )
c
c     If the boundary is of Neumann type, set the ghost cell value to
c     that of the outermost point in the valid data (2nd order accurate)
c     and then fill the "den" array with the value "1"
c     
c     
c     If flagbc==1:
c     
c     If the boundary is of Dirichlet type, construct a polynomial
c     interpolation through the boundary location and internal points
c     (at locations x(-1:len-2) that generates the ghost cell value (at
c     location xInt).  Then fill the ghost cell with the interpolated value.
c     If flagden==1, load the "den" array with the interpolation
c     coefficient corresponding to outermost point in the valid region
c     ( the coef(0) corresponding to the location x(0) )
c      
c     Note: 
c     The bc type = LO_REFLECT_ODD is a special type of boundary condition.
c      
      integer maxorder
      integer nc, cdir, flagden, flagbc
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer DIMDEC(phi)
      REAL_T phi(DIMV(phi),nc)
      integer DIMDEC(den)
      REAL_T den(DIMV(den))
      integer DIMDEC(bcval)
      REAL_T bcval(DIMV(bcval),nc)
      integer DIMDEC(mask)
      integer mask(DIMV(mask))
      integer bct
      REAL_T bcl
      REAL_T h(BL_SPACEDIM)
c
      integer i,n
      logical is_dirichlet
      logical is_neumann
c
      integer lenx
      integer m
c
      integer Lmaxorder
      integer maxmaxorder
      parameter(maxmaxorder=4)
      REAL_T x(-1:maxmaxorder-2)
      REAL_T coef(-1:maxmaxorder-2)
      REAL_T xInt
c
      is_dirichlet(i) = ( i .eq. LO_DIRICHLET )
      is_neumann(i)   = ( i .eq. LO_NEUMANN )
c
      if ( maxorder .eq. -1 ) then
         Lmaxorder = maxmaxorder
      else
         Lmaxorder = MIN(maxorder,maxmaxorder)
      end if
      lenx = MIN(hi(1)-lo(1), Lmaxorder-2)
c
c     TODO:
c     In order for this to work with growing multigrid, must
c     sort xa[] because it is possible for the xb value to lay
c     within this range.
c     
c     The Left face of the grid
c
      if(cdir .eq. 0) then
         if (is_neumann(bct)) then
            do n = 1, nc
               phi(lo(1)-1,n) = merge(
     $              phi(lo(1),n),
     $              phi(lo(1)-1,n),
     $              mask(lo(1)-1) .gt. 0)
            end do
            if ( flagden .eq. 1) then
               den(lo(1)) = 1.0D0
            end if
         else if (is_dirichlet(bct)) then
            do m=0,lenx
               x(m) = m + 0.5D0
            end do
            x(-1) = - bcl/h(1)
            xInt = - 0.5D0
            call polyInterpCoeff(xInt, x, lenx+2, coef)
            do n = 1, nc
               if ( flagbc .eq. 1 ) then
                  phi(lo(1)-1,n) = merge(
     $                 bcval(lo(1)-1,n)*coef(-1),
     $                 phi(lo(1)-1,n),
     $                 mask(lo(1)-1) .gt. 0)
               else
                  phi(lo(1)-1,n) = merge(
     $                 0.0D0,
     $                 phi(lo(1)-1,n),
     $                 mask(lo(1)-1) .gt. 0)
               end if
               do m = 0, lenx
                  phi(lo(1)-1,n) = merge(
     $                 phi(lo(1)-1,n)
     $                 + phi(lo(1)+m,n)*coef(m),
     $                 phi(lo(1)-1,n),
     $                 mask(lo(1)-1) .gt. 0)
               end do
            end do
            if ( flagden .eq. 1 ) then
               den(lo(1)) = merge(coef(0), 0.0D0,
     $              mask(lo(1)-1) .gt. 0)
            end if

         else if ( bct .eq. LO_REFLECT_ODD ) then

            do n = 1, nc
               phi(lo(1)-1,n) = merge(
     $             -phi(lo(1),n),
     $              phi(lo(1)-1,n),
     $              mask(lo(1)-1) .gt. 0)
            end do
            if ( flagden .eq. 1 ) then
               den(lo(1)) = merge(-1.0D0, 0.0D0,
     $              mask(lo(1)-1) .gt. 0)
            end if

         else 
            print *,'UNKNOWN BC ON LEFT FACE IN APPLYBC'
            call bl_error("stop")
         end if
      end if
c     
c     The Right face of the grid
c
      if(cdir .eq. 2) then
         if(is_neumann(bct)) then
            do n = 1, nc
               phi(hi(1)+1,n) = merge(
     $              phi(hi(1),n),
     $              phi(hi(1)+1,n),
     $              mask(hi(1)+1) .gt. 0)
            end do
	    if ( flagden .eq. 1 ) then
                den(hi(1)) = 1.0D0
	    end if
         else if (is_dirichlet(bct)) then
            do m=0,lenx
               x(m) = m + 0.5D0
            end do
            x(-1) = - bcl/h(1)
            xInt = - 0.5D0
            call polyInterpCoeff(xInt, x, lenx+2, coef)
            do n = 1, nc 
               if ( flagbc .eq. 1 ) then
                  phi(hi(1)+1,n) = merge(
     $                 bcval(hi(1)+1,n)*coef(-1),
     $                 phi(hi(1)+1,n),
     $                 mask(hi(1)+1) .gt. 0)
               else
                  phi(hi(1)+1,n) = merge(
     $                 0.0D0,
     $                 phi(hi(1)+1,n),
     $                 mask(hi(1)+1) .gt. 0)
               end if
               do m = 0, lenx
                  phi(hi(1)+1,n) = merge(
     $                 phi(hi(1)+1,n)
     $                 + phi(hi(1)-m,n)*coef(m),
     $                 phi(hi(1)+1,n),
     $                 mask(hi(1)+1) .gt. 0)
               end do
            end do
            if ( flagden .eq. 1 ) then
               den(hi(1))   = merge(coef(0), 0.0D0,
     $              mask(hi(1)+1) .gt. 0)
            end if

         else if ( bct .eq. LO_REFLECT_ODD ) then

            do n = 1, nc
               phi(hi(1)+1,n) = merge(
     $             -phi(hi(1),n),
     $              phi(hi(1)+1,n),
     $              mask(hi(1)+1) .gt. 0)
            end do
            if ( flagden .eq. 1 ) then
               den(hi(1)) = merge(-1.0D0, 0.0D0,
     $              mask(hi(1)+1) .gt. 0)
            end if

         else
            print *,'UNKNOWN BC ON RIGHT FACE IN APPLYBC'
            call bl_error("stop")
         end if
      end if
c
      end
