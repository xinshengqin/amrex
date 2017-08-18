module advance_2d
    implicit none
    private
    public :: compute_flux_auto, update_phi
    contains

subroutine compute_flux (lo, hi, phi, philo, phihi, &
     fluxx, fxlo, fxhi, fluxy, fylo, fyhi, dx & 
#ifdef CUDA
     , idx, device_id &
#endif
     ) bind(C, name="compute_flux")

  use amrex_fort_module, only : amrex_real
#ifdef CUDA
  use cuda_module, only: threads_and_blocks, stream_from_index
  use cuda_module, only: numThreads, numBlocks, cuda_streams 
#endif
  implicit none

  integer lo(2), hi(2), philo(2), phihi(2), fxlo(2), fxhi(2), fylo(2), fyhi(2)
  real(amrex_real), intent(in)    :: phi  (philo(1):phihi(1),philo(2):phihi(2))
  real(amrex_real), intent(inout) :: fluxx( fxlo(1): fxhi(1), fxlo(2): fxhi(2))
  real(amrex_real), intent(inout) :: fluxy( fylo(1): fyhi(1), fylo(2): fyhi(2))
  real(amrex_real), intent(in)    :: dx(2)
  integer :: idir
#ifdef CUDA
  attributes(device) :: phi, fluxx, fluxy
  integer, value :: idx, device_id
#endif

#ifdef CUDA
    ! call threads_and_blocks(lo, hi+1, numBlocks, numThreads)
    ! call compute_flux_kernel<<<numBlocks, numThreads, 0, cuda_streams(stream_from_index(idx),device_id)>>> &
    !         (lo(1), lo(2), hi(1), hi(2), &
    !          phi, philo(1), philo(2), phihi(1), phihi(2), &
    !          fluxx, fxlo(1), fxlo(2), fxhi(1), fxhi(2), &
    !          fluxy, fylo(1), fylo(2), fyhi(1), fyhi(2), &
    !          dx(1), dx(2))
    call compute_flux_auto(lo, hi, phi, philo, phihi, &
         fluxx, fxlo, fxhi, fluxy, fylo, fyhi, dx & 
         , idx, device_id) 
#else
    ! x-fluxes
    idir = 1
    call compute_flux_doit(lo, [hi(1)+1, hi(2)], phi, philo, phihi, fluxx, fxlo, fxhi, dx, idir)
    ! y-fluxes
    idir = 2
    call compute_flux_doit(lo, [hi(1), hi(2)+1], phi, philo, phihi, fluxy, fylo, fyhi, dx, idir)
#endif


end subroutine compute_flux

#ifdef CUDA
attributes(global) &
subroutine compute_flux_kernel(lox, loy, hix, hiy, &
    phi, philox, philoy, phihix, phihiy, &
    fx, fxlox, fxloy, fxhix, fxhiy, &
    fy, fylox, fyloy, fyhix, fyhiy, &
    dx, dy)

  use amrex_fort_module, only : amrex_real, map_cuda_threads_to_space
  implicit none
  integer, value, intent(in) :: lox, loy, hix, hiy
  integer, value, intent(in) :: philox, philoy, phihix, phihiy
  integer, value, intent(in) :: fxlox, fxloy, fxhix, fxhiy
  integer, value, intent(in) :: fylox, fyloy, fyhix, fyhiy
  real(amrex_real), value, intent(in) :: dx, dy
  real(amrex_real), intent(in)    :: phi(philox:phihix,philoy:phihiy)
  real(amrex_real), intent(inout) :: fx(fxlox: fxhix, fxloy: fxhiy)
  real(amrex_real), intent(inout) :: fy(fylox: fyhix, fyloy: fyhiy)

  ! local variable
  integer :: idir
  integer :: i,j
  logical :: has_work
  ! x-fluxes
  idir = 1
  call map_cuda_threads_to_space([lox,loy], [hix+1, hiy], i, j, has_work)
  if (has_work) then
      call compute_flux_gpu(i, j, &
      phi, [philox, philoy] , [phihix, phihiy] , & 
      fx, [fxlox, fxloy] , [fxhix, fxhiy] , & 
      [dx, dy] , idir)
  endif
  ! y-fluxes
  idir = 2
  call map_cuda_threads_to_space([lox,loy], [hix, hiy+1], i, j, has_work)
  if (has_work) then
  call compute_flux_gpu(i, j, &
      phi, [philox, philoy] , [phihix, phihiy] , & 
      fy, [fylox, fyloy] , [fyhix, fyhiy] , & 
      [dx, dy] , idir)
  endif

end subroutine compute_flux_kernel
#endif



subroutine update_phi (lo, hi, phiold, polo, pohi, phinew, pnlo, pnhi, &
     fluxx, fxlo, fxhi, fluxy, fylo, fyhi, dx, dt &
#ifdef CUDA
     , idx, device_id &
#endif
    ) bind(C, name="update_phi")

  use amrex_fort_module, only : amrex_real
#ifdef CUDA
  use cuda_module, only: threads_and_blocks, stream_from_index
  use cuda_module, only: numThreads, numBlocks, cuda_streams 
#endif
  implicit none

  integer lo(2), hi(2), polo(2), pohi(2), pnlo(2), pnhi(2), fxlo(2), fxhi(2), fylo(2), fyhi(2)
  real(amrex_real), intent(in)    :: phiold(polo(1):pohi(1),polo(2):pohi(2))
  real(amrex_real), intent(inout) :: phinew(pnlo(1):pnhi(1),pnlo(2):pnhi(2))
  real(amrex_real), intent(in   ) :: fluxx (fxlo(1):fxhi(1),fxlo(2):fxhi(2))
  real(amrex_real), intent(in   ) :: fluxy (fylo(1):fyhi(1),fylo(2):fyhi(2))
  real(amrex_real), intent(in)    :: dx(2)
  real(amrex_real), value, intent(in)    :: dt

#ifdef CUDA
  attributes(device) :: phiold, phinew, fluxx, fluxy
  integer, value :: idx, device_id
#endif


#ifdef CUDA
    call threads_and_blocks(lo, hi, numBlocks, numThreads)
    call update_phi_kernel<<<numBlocks, numThreads, 0, cuda_streams(stream_from_index(idx),device_id)>>> &
            (lo(1), lo(2), hi(1), hi(2), &
             phiold, polo(1), polo(2), pohi(1), pohi(2), &
             phinew, pnlo(1), pnlo(2), pnhi(1), pnhi(2), &
             fluxx, fxlo(1), fxlo(2), fxhi(1), fxhi(2), &
             fluxy, fylo(1), fylo(2), fyhi(1), fyhi(2), &
             dx(1), dx(2), dt)
#else
    call update_phi_doit(lo, hi, phiold, polo, pohi, phinew, pnlo, pnhi, &
     fluxx, fxlo, fxhi, fluxy, fylo, fyhi, dx, dt)
#endif


end subroutine update_phi


#ifdef CUDA
attributes(global) &
subroutine update_phi_kernel(lox, loy, hix, hiy, &
    phiold, polox, poloy, pohix, pohiy, &
    phinew, pnlox, pnloy, pnhix, pnhiy, &
    fx, fxlox, fxloy, fxhix, fxhiy, &
    fy, fylox, fyloy, fyhix, fyhiy, &
    dx, dy, dt)

  use amrex_fort_module, only : amrex_real
  implicit none
  integer, value, intent(in) :: lox, loy, hix, hiy
  integer, value, intent(in) :: polox, poloy, pohix, pohiy
  integer, value, intent(in) :: pnlox, pnloy, pnhix, pnhiy
  integer, value, intent(in) :: fxlox, fxloy, fxhix, fxhiy
  integer, value, intent(in) :: fylox, fyloy, fyhix, fyhiy
  real(amrex_real), value, intent(in) :: dx, dy, dt
  real(amrex_real), intent(in)    :: phiold(polox:pohix,poloy:pohiy)
  real(amrex_real), intent(inout) :: phinew(pnlox:pnhix,pnloy:pnhiy)
  real(amrex_real), intent(in)    :: fx(fxlox: fxhix, fxloy: fxhiy)
  real(amrex_real), intent(in)    :: fy(fylox: fyhix, fyloy: fyhiy)

  call update_phi_doit([lox,loy], [hix,hiy], &
      phiold, [polox, poloy] , [pohix, pohiy] , & 
      phinew, [pnlox, pnloy] , [pnhix, pnhiy] , & 
      fx, [fxlox, fxloy] , [fxhix, fxhiy] , & 
      fy, [fylox, fyloy] , [fyhix, fyhiy] , & 
      [dx, dy] , dt)
end subroutine update_phi_kernel
#endif


#ifdef CUDA
  attributes(device) &
#endif
subroutine compute_flux_doit (lo, hi, phi, p_lo, p_hi, flx, f_lo, f_hi, dx, idir)

  use amrex_fort_module, only: rt => amrex_real, get_loop_bounds

  implicit none

  integer,  intent(in   ) :: lo(2), hi(2), p_lo(2), p_hi(2), f_lo(2), f_hi(2)
  real(rt), intent(in   ) :: phi(p_lo(1):p_hi(1),p_lo(2):p_hi(2))
  real(rt), intent(inout) :: flx(f_lo(1):f_hi(1),f_lo(2):f_hi(2))
  real(rt), intent(in   ) :: dx(2)
  integer,  intent(in   ) :: idir

  ! local variables
  integer :: i, j
  integer :: blo(2),bhi(2)

  call get_loop_bounds(blo, bhi, lo, hi)

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

#ifdef CUDA
  attributes(device) &
#endif
subroutine update_phi_doit (lo, hi, phiold, polo, pohi, phinew, pnlo, pnhi, &
     fluxx, fxlo, fxhi, fluxy, fylo, fyhi, dx, dt) 
  use amrex_fort_module, only : amrex_real, get_loop_bounds
  implicit none

  integer :: lo(2), hi(2), polo(2), pohi(2), pnlo(2), pnhi(2), fxlo(2), fxhi(2), fylo(2), fyhi(2)
  real(amrex_real), intent(in)    :: phiold(polo(1):pohi(1),polo(2):pohi(2))
  real(amrex_real), intent(inout) :: phinew(pnlo(1):pnhi(1),pnlo(2):pnhi(2))
  real(amrex_real), intent(in   ) :: fluxx (fxlo(1):fxhi(1),fxlo(2):fxhi(2))
  real(amrex_real), intent(in   ) :: fluxy (fylo(1):fyhi(1),fylo(2):fyhi(2))
  real(amrex_real), intent(in)    :: dx(2)
  real(amrex_real), value, intent(in) :: dt
  ! local variables
  integer :: i, j
  integer :: blo(2),bhi(2)
  real(amrex_real) :: dtdx(2)
  dtdx = dt/dx

  call get_loop_bounds(blo, bhi, lo, hi)
  do    j = blo(2), bhi(2)
     do i = blo(1), bhi(1)

        phinew(i,j) = phiold(i,j) &
             + dtdx(1) * (fluxx(i+1,j  ) - fluxx(i,j)) &
             + dtdx(2) * (fluxy(i  ,j+1) - fluxy(i,j))

     end do
  end do
end subroutine update_phi_doit

#ifdef CUDA
  attributes(device) &
subroutine compute_flux_gpu (i, j, phi, p_lo, p_hi, flx, f_lo, f_hi, dx, idir)

  use amrex_fort_module, only: rt => amrex_real, get_loop_bounds

  implicit none

  integer,  intent(in   ) :: p_lo(2), p_hi(2), f_lo(2), f_hi(2)
  real(rt), intent(in   ) :: phi(p_lo(1):p_hi(1),p_lo(2):p_hi(2))
  real(rt), intent(inout) :: flx(f_lo(1):f_hi(1),f_lo(2):f_hi(2))
  real(rt), intent(in   ) :: dx(2)
  integer,  intent(in   ) :: idir
  integer,  intent(in   ) :: i, j

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
end subroutine compute_flux_gpu
#endif

! test pgfortran kernel loop directives
subroutine compute_flux_auto(lo, hi, phi, philo, phihi, &
     fluxx, fxlo, fxhi, fluxy, fylo, fyhi, dx & 
     , idx, device_id &
     )

  use amrex_fort_module, only : amrex_real
  use cuda_module, only: threads_and_blocks, stream_from_index
  use cuda_module, only: numThreads, numBlocks, cuda_streams 

  implicit none

  integer lo(2), hi(2), philo(2), phihi(2), fxlo(2), fxhi(2), fylo(2), fyhi(2)
  real(amrex_real), intent(in)    :: phi  (philo(1):phihi(1),philo(2):phihi(2))
  real(amrex_real), intent(inout) :: fluxx( fxlo(1): fxhi(1), fxlo(2): fxhi(2))
  real(amrex_real), intent(inout) :: fluxy( fylo(1): fyhi(1), fylo(2): fyhi(2))
  real(amrex_real), intent(in)    :: dx(2)
  integer :: idir
  attributes(device) :: phi, fluxx, fluxy
  integer, value :: idx, device_id


  ! locals
  real(amrex_real) :: dx_x, dx_y
  integer :: i,j


  dx_x = dx(1)
  dx_y = dx(2)

  ! x-fluxes
  !$cuf kernel do(2) <<<*, (16,16), 0, cuda_streams(stream_from_index(idx),device_id)>>> 
  do    j = lo(2), hi(2)
     do i = lo(1), hi(1)+1
         fluxx(i,j) = ( phi(i,j) - phi(i-1,j) ) / dx_x
     end do
  end do

  ! y-fluxes
  !$cuf kernel do(2) <<<*, (16,16), 0, cuda_streams(stream_from_index(idx),device_id)>>> 
  do    j = lo(2), hi(2)+1
     do i = lo(1), hi(1)
        fluxy(i,j) = ( phi(i,j) - phi(i,j-1) ) / dx_y
     end do
  end do

end subroutine 

end module advance_2d
