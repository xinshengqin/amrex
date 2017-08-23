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
  ! attributes(device) :: uin, uout, vx, vy, flxx, flxy
  integer :: uin_size, uout_size, vx_size, vy_size, flxx_size, flxy_size, cudaResult
#endif 
  double precision, device, allocatable :: uin_d(:,:)
  double precision, device, allocatable :: uout_d(:,:)
  double precision, device, allocatable :: vx_d(:,:)
  double precision, device, allocatable :: vy_d(:,:)
  double precision, device, allocatable :: flxx_d(:,:)
  double precision, device, allocatable :: flxy_d(:,:)
  allocate(uin_d (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2)))
  allocate(uout_d(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2)))
  allocate(vx_d  (vx_lo(1):vx_hi(1),vx_lo(2):vx_hi(2)))
  allocate(vy_d  (vy_lo(1):vy_hi(1),vy_lo(2):vy_hi(2)))
  allocate(flxx_d(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2)))
  allocate(flxy_d(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2)))

  uin_size = (ui_hi(1)-ui_lo(1)+1)*(ui_hi(2)-ui_lo(2)+1)
  uout_size = (uo_hi(1)-uo_lo(1)+1)*(uo_hi(2)-uo_lo(2)+1)
  vx_size = (vx_hi(1)-vx_lo(1)+1)*(vx_hi(2)-vx_lo(2)+1)
  vy_size = (vy_hi(1)-vy_lo(1)+1)*(vy_hi(2)-vy_lo(2)+1)
  flxx_size = (fx_hi(1)-fx_lo(1)+1)*(fx_hi(2)-fx_lo(2)+1)
  flxy_size = (fy_hi(1)-fy_lo(1)+1)*(fy_hi(2)-fy_lo(2)+1)


  cudaResult = cudaMemcpyAsync(uin_d, uin, uin_size, cudaMemcpyHostToDevice, cuda_streams(stream_from_index(idx),device_id))
  cudaResult = cudaMemcpyAsync(vx_d, vx, vx_size, cudaMemcpyHostToDevice, cuda_streams(stream_from_index(idx),device_id))
  cudaResult = cudaMemcpyAsync(vy_d, vy, vy_size, cudaMemcpyHostToDevice, cuda_streams(stream_from_index(idx),device_id))
!   call advect_cuda(time, lo, hi, &
!      &            uin , ui_lo, ui_hi, &
!      &            uout, uo_lo, uo_hi, &
!      &            vx  , vx_lo, vx_hi, &
!      &            vy  , vy_lo, vy_hi, &
!      &            flxx, fx_lo, fx_hi, &
!      &            flxy, fy_lo, fy_hi, &
!      &            dx,dt &
! #ifdef CUDA
!                   , idx, device_id &
! #endif
!                   )
  call advect_cuda(time, lo, hi, &
     &            uin_d , ui_lo, ui_hi, &
     &            uout_d, uo_lo, uo_hi, &
     &            vx_d  , vx_lo, vx_hi, &
     &            vy_d  , vy_lo, vy_hi, &
     &            flxx_d, fx_lo, fx_hi, &
     &            flxy_d, fy_lo, fy_hi, &
     &            dx,dt &
#ifdef CUDA
                  , idx, device_id &
#endif
                  )

  cudaResult = cudaMemcpyAsync(uout, uout_d, uout_size, cudaMemcpyDeviceToHost, cuda_streams(stream_from_index(idx),device_id))
  cudaResult = cudaMemcpyAsync(flxx, flxx_d, flxx_size, cudaMemcpyDeviceToHost, cuda_streams(stream_from_index(idx),device_id))
  cudaResult = cudaMemcpyAsync(flxy, flxy_d, flxy_size, cudaMemcpyDeviceToHost, cuda_streams(stream_from_index(idx),device_id))
  deallocate(uin_d)
  deallocate(uout_d)
  deallocate(vx_d)
  deallocate(vy_d)
  deallocate(flxx_d)
  deallocate(flxy_d)
end subroutine advect
end module advect_module
