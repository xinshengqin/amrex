module advect_module
#ifdef CUDA
    use cuda_module, only: threads_and_blocks, stream_from_index
    use cuda_module, only: numThreads, numBlocks, cuda_streams 
    use iso_c_binding
#endif
    implicit none
    contains
subroutine advect(time, lo, hi, &
     &            uin , ui_lo, ui_hi, &
     &            uout, uo_lo, uo_hi, &
     &            vx  , vx_lo, vx_hi, &
     &            vy  , vy_lo, vy_hi, &
     &            vz  , vz_lo, vz_hi, &
     &            flxx, fx_lo, fx_hi, &
     &            flxy, fy_lo, fy_hi, &
     &            flxz, fz_lo, fz_hi, &
     &            dx,dt &
#ifdef CUDA
                  , idx, device_id, tag &
#endif
                ) bind(C, name="advect")
  
  use advect_doit_module, only: advect_doit

  implicit none

  integer, intent(in) :: lo(3), hi(3)
  double precision, intent(in) :: dx(3), dt, time
  integer, intent(in) :: ui_lo(3), ui_hi(3)
  integer, intent(in) :: uo_lo(3), uo_hi(3)
  integer, intent(in) :: vx_lo(3), vx_hi(3)
  integer, intent(in) :: vy_lo(3), vy_hi(3)
  integer, intent(in) :: vz_lo(3), vz_hi(3)
  integer, intent(in) :: fx_lo(3), fx_hi(3)
  integer, intent(in) :: fy_lo(3), fy_hi(3)
  integer, intent(in) :: fz_lo(3), fz_hi(3)
  double precision, intent(in   ) :: uin (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3))
  double precision, intent(inout) :: uout(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2),uo_lo(3):uo_hi(3))
  double precision, intent(in   ) :: vx  (vx_lo(1):vx_hi(1),vx_lo(2):vx_hi(2),vx_lo(3):vx_hi(3))
  double precision, intent(in   ) :: vy  (vy_lo(1):vy_hi(1),vy_lo(2):vy_hi(2),vy_lo(3):vy_hi(3))
  double precision, intent(in   ) :: vz  (vz_lo(1):vz_hi(1),vz_lo(2):vz_hi(2),vz_lo(3):vz_hi(3))
  double precision, intent(  out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3))
  double precision, intent(  out) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3))
  double precision, intent(  out) :: flxz(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3))

  integer :: i, j, k
  integer :: glo(3), ghi(3)
  double precision :: dtdx(3), umax, vmax, wmax

#ifdef CUDA
  integer, value, intent(in) :: idx, device_id
  attributes(device) :: uin, uout, vx, vy, vz, flxx, flxy, flxz
  integer(kind=c_intptr_t), value, intent(in) :: tag
#endif 

  call advect_doit(time, lo, hi, &
     &             uin , ui_lo, ui_hi, &
     &             uout, uo_lo, uo_hi, &
     &             vx  , vx_lo, vx_hi, &
     &             vy  , vy_lo, vy_hi, &
     &             vz  , vz_lo, vz_hi, &
     &             flxx, fx_lo, fx_hi, &
     &             flxy, fy_lo, fy_hi, &
     &             flxz, fz_lo, fz_hi, &
     &             dx,dt &
#ifdef CUDA
                   , idx, device_id, tag &
#endif
              )

end subroutine advect
end module advect_module
