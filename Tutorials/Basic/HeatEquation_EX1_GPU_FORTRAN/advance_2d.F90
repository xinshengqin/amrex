subroutine compute_flux (lo, hi, phi, philo, phihi, &
     fluxx, fxlo, fxhi, fluxy, fylo, fyhi, dx) bind(C, name="compute_flux")

  use amrex_fort_module, only : amrex_real
  implicit none

  integer lo(2), hi(2), philo(2), phihi(2), fxlo(2), fxhi(2), fylo(2), fyhi(2)
  real(amrex_real), intent(in)    :: phi  (philo(1):phihi(1),philo(2):phihi(2))
  real(amrex_real), intent(inout) :: fluxx( fxlo(1): fxhi(1), fxlo(2): fxhi(2))
  real(amrex_real), intent(inout) :: fluxy( fylo(1): fyhi(1), fylo(2): fyhi(2))
  real(amrex_real), intent(in)    :: dx(2)
  integer :: idir

  ! x-fluxes
  idir = 1
  call compute_flux_doit(lo, [hi(1)+1, hi(2)], phi, philo, phihi, fluxx, fxlo, fxhi, dx, idir)
  ! y-fluxes
  idir = 2
  call compute_flux_doit(lo, [hi(1), hi(2)+1], phi, philo, phihi, fluxy, fylo, fyhi, dx, idir)

  ! do    j = lo(2), hi(2)
  !    do i = lo(1), hi(1)+1
  !       fluxx(i,j) = ( phi(i,j) - phi(i-1,j) ) / dx(1)
  !    end do
  ! end do

  ! ! y-fluxes
  ! do    j = lo(2), hi(2)+1
  !    do i = lo(1), hi(1)
  !       fluxy(i,j) = ( phi(i,j) - phi(i,j-1) ) / dx(2)
  !    end do
  ! end do

end subroutine compute_flux


subroutine update_phi (lo, hi, phiold, polo, pohi, phinew, pnlo, pnhi, &
     fluxx, fxlo, fxhi, fluxy, fylo, fyhi, dx, dt) bind(C, name="update_phi")

  use amrex_fort_module, only : amrex_real
  implicit none

  integer lo(2), hi(2), polo(2), pohi(2), pnlo(2), pnhi(2), fxlo(2), fxhi(2), fylo(2), fyhi(2)
  real(amrex_real), intent(in)    :: phiold(polo(1):pohi(1),polo(2):pohi(2))
  real(amrex_real), intent(inout) :: phinew(pnlo(1):pnhi(1),pnlo(2):pnhi(2))
  real(amrex_real), intent(in   ) :: fluxx (fxlo(1):fxhi(1),fxlo(2):fxhi(2))
  real(amrex_real), intent(in   ) :: fluxy (fylo(1):fyhi(1),fylo(2):fyhi(2))
  real(amrex_real), intent(in)    :: dx(2)
  real(amrex_real), value         :: dt

  ! local variables
  integer i,j
  real(amrex_real) :: dtdx(2)

  dtdx = dt/dx

  do    j = lo(2), hi(2)
     do i = lo(1), hi(1)

        phinew(i,j) = phiold(i,j) &
             + dtdx(1) * (fluxx(i+1,j  ) - fluxx(i,j)) &
             + dtdx(2) * (fluxy(i  ,j+1) - fluxy(i,j))

     end do
  end do

end subroutine update_phi

#ifdef CUDA
  attributes(device) &
#endif
subroutine compute_flux_doit (lo, hi, phi, p_lo, p_hi, flx, f_lo, f_hi, dx, idir)

  use amrex_fort_module, only: rt => amrex_real, get_loop_bounds_2d

  implicit none

  integer,  intent(in   ) :: lo(2), hi(2), p_lo(2), p_hi(2), f_lo(2), f_hi(2)
  real(rt), intent(in   ) :: phi(p_lo(1):p_hi(1),p_lo(2):p_hi(2))
  real(rt), intent(inout) :: flx(f_lo(1):f_hi(1),f_lo(2):f_hi(2))
  real(rt), intent(in   ) :: dx(2)
  integer,  intent(in   ) :: idir

  ! local variables
  integer :: i, j
  integer :: blo(2),bhi(2)

  call get_loop_bounds_2d(blo, bhi, lo, hi)

  do    j = blo(2), bhi(2)
     do i = blo(1), bhi(1)
        if (idir .eq. 1) then
           ! x-flux
           flx(i,j) = ( phi(i,j) - phi(i-1,j) ) / dx(1)
        else if (idir .eq. 2) then
           ! y-flux
           flx(i,j) = ( phi(i,j) - phi(i,j-1) ) / dx(2)
	else
	   ! this should not happen
	   stop 1
        endif
     end do
  end do
end subroutine compute_flux_doit

