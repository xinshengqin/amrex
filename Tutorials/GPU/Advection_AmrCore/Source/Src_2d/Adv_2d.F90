module advect_module
    use cudafor, only: cudaMemcpyAsync, cudaMemcpyHostToDevice, & 
        cudaMemcpyDeviceToHost, cudaDeviceSynchronize, cudaMemcpy
    use cuda_module, only: threads_and_blocks, stream_from_index
    use cuda_module, only: numThreads, numBlocks, cuda_streams 
    use amrex_fort_module, only: amrex_real
    implicit none
    contains
subroutine advect(time, lo, hi, &
     &            uin , ui_lo, ui_hi, &
     &            uout, uo_lo, uo_hi, &
     &            vx  , vx_lo, vx_hi, &
     &            vy  , vy_lo, vy_hi, &
     &            flxx, fx_lo, fx_hi, &
     &            flxy, fy_lo, fy_hi, &
     &            dx,dt &
#ifdef CUDA
                  , idx, device_id &
#endif
                 ) bind(C, name="advect")
  
  use advect_cuda_module, only: advect_cuda
  implicit none

  integer, intent(in) :: lo(2), hi(2)
  double precision, intent(in) :: dx(2), dt, time
  integer, intent(in) :: ui_lo(2), ui_hi(2)
  integer, intent(in) :: uo_lo(2), uo_hi(2)
  integer, intent(in) :: vx_lo(2), vx_hi(2)
  integer, intent(in) :: vy_lo(2), vy_hi(2)
  integer, intent(in) :: fx_lo(2), fx_hi(2)
  integer, intent(in) :: fy_lo(2), fy_hi(2)
  double precision, intent(in   ) :: uin (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2))
  double precision, intent(inout) :: uout(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2))
  double precision, intent(in   ) :: vx  (vx_lo(1):vx_hi(1),vx_lo(2):vx_hi(2))
  double precision, intent(in   ) :: vy  (vy_lo(1):vy_hi(1),vy_lo(2):vy_hi(2))
  double precision, intent(  out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2))
  double precision, intent(  out) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2))


#ifdef CUDA
  integer, value, intent(in) :: idx, device_id
  attributes(device) :: uin, uout, vx, vy, flxx, flxy
  ! integer :: uin_size, uout_size, vx_size, vy_size, flxx_size, flxy_size, cudaResult
#endif 
  call advect_cuda(time, lo, hi, &
     &            uin , ui_lo, ui_hi, &
     &            uout, uo_lo, uo_hi, &
     &            vx  , vx_lo, vx_hi, &
     &            vy  , vy_lo, vy_hi, &
     &            flxx, fx_lo, fx_hi, &
     &            flxy, fy_lo, fy_hi, &
     &            dx,dt &
#ifdef CUDA
                  , idx, device_id &
#endif
                  )
end subroutine advect
end module advect_module
