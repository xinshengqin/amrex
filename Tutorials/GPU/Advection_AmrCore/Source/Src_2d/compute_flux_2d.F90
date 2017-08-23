module compute_flux_module

  implicit none

  private

  public :: compute_flux_2d

contains

  subroutine compute_flux_2d(lo, hi, dt, dx, &
                             phi,ph_lo,ph_hi, &
                             umac,  u_lo,  u_hi, &
                             vmac,  v_lo,  v_hi, &
                             flxx, fx_lo, fx_hi, &
                             flxy, fy_lo, fy_hi, &
                             phix_1d, phiy_1d, phix, phiy, slope, glo, ghi &
#ifdef CUDA
                             , idx, device_id &
#endif
                             )

#ifdef CUDA
    use cudafor, only: cudaMemcpyAsync, cudaMemcpyHostToDevice, & 
        cudaMemcpyDeviceToHost, cudaDeviceSynchronize, cudaMemcpy
    use cuda_module, only: threads_and_blocks, stream_from_index
    use cuda_module, only: numThreads, numBlocks, cuda_streams 
    ! use slope_module_cuda, only: slopex
    use slope_module_cuda, only: slopex_cuf, slopey_cuf
#else
    use slope_module, only: slopex, slopey
#endif

    integer, intent(in) :: lo(2), hi(2), glo(2), ghi(2)
    double precision, intent(in) :: dt, dx(2)
    integer, intent(in) :: ph_lo(2), ph_hi(2)
    integer, intent(in) ::  u_lo(2),  u_hi(2)
    integer, intent(in) ::  v_lo(2),  v_hi(2)
    integer, intent(in) :: fx_lo(2), fx_hi(2)
    integer, intent(in) :: fy_lo(2), fy_hi(2)
    double precision, intent(in   ) :: phi (ph_lo(1):ph_hi(1),ph_lo(2):ph_hi(2))
    double precision, intent(in   ) :: umac( u_lo(1): u_hi(1), u_lo(2): u_hi(2))
    double precision, intent(in   ) :: vmac( v_lo(1): v_hi(1), v_lo(2): v_hi(2))
    double precision, intent(  out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2))
    double precision, intent(  out) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2))
    double precision, dimension(glo(1):ghi(1),glo(2):ghi(2)) :: &
         phix_1d, phiy_1d, phix, phiy, slope
         
    integer :: i, j, k
    double precision :: hdtdx(2)

#ifdef CUDA
    attributes(device) :: phi, umac, vmac, flxx, flxy, phix_1d, phiy_1d, phix, phiy, slope
    ! device 
    integer, intent(in) :: idx, device_id
    double precision :: hdtdx_x, hdtdx_y
    integer :: cudaResult
    ! integer :: phi_size, umac_size, vmac_size, slope_size, phix_1d_size, phiy_1d_size, phix_size, phiy_size, &
    !     flxx_size, flxy_size
    ! double precision, allocatable, device :: phi_d(:,:), slope_d(:,:), umac_d(:,:), vmac_d(:,:)
    ! double precision, allocatable, device :: phix_1d_d(:,:), phiy_1d_d(:,:)
    ! double precision, allocatable, device :: phix_d(:,:), phiy_d(:,:)
    ! double precision, allocatable, device :: flxx_d(:,:), flxy_d(:,:)
    ! TODO: can replace this with our own device memory allocator
    ! allocate(phix_1d_d(glo(1):ghi(1),glo(2):ghi(2)))
    ! allocate(phiy_1d_d(glo(1):ghi(1),glo(2):ghi(2)))
    ! allocate(phix_d(glo(1):ghi(1),glo(2):ghi(2)))
    ! allocate(phiy_d(glo(1):ghi(1),glo(2):ghi(2)))
    ! allocate(slope_d(glo(1):ghi(1),glo(2):ghi(2)))
    ! allocate(phi_d(ph_lo(1):ph_hi(1),ph_lo(2):ph_hi(2)))
    ! allocate(umac_d( u_lo(1): u_hi(1), u_lo(2): u_hi(2)))
    ! allocate(vmac_d( v_lo(1): v_hi(1), v_lo(2): v_hi(2)))
    ! allocate(flxx_d( fx_lo(1): fx_hi(1), fx_lo(2): fx_hi(2)))
    ! allocate(flxy_d( fy_lo(1): fy_hi(1), fy_lo(2): fy_hi(2)))
#endif


    ! phi_size = (ph_hi(1)-ph_lo(1)+1) * (ph_hi(2)-ph_lo(2)+1)
    ! slope_size = (ghi(1)-glo(1)+1) * (ghi(2)-glo(2)+1)
    ! umac_size = (u_hi(1)-u_lo(1)+1) * (u_hi(2)-u_lo(2)+1)
    ! vmac_size = (v_hi(1)-v_lo(1)+1) * (v_hi(2)-v_lo(2)+1)
    ! phix_1d_size = (ghi(1)-glo(1)+1) * (ghi(2)-glo(2)+1)
    ! phiy_1d_size = (ghi(1)-glo(1)+1) * (ghi(2)-glo(2)+1)
    ! phix_size = (ghi(1)-glo(1)+1) * (ghi(2)-glo(2)+1)
    ! phiy_size = (ghi(1)-glo(1)+1) * (ghi(2)-glo(2)+1)
    ! flxx_size = (fx_hi(1)-fx_lo(1)+1) * (fx_hi(2)-fx_lo(2)+1)
    ! flxy_size = (fy_hi(1)-fy_lo(1)+1) * (fy_hi(2)-fy_lo(2)+1)

    ! cudaResult = cudaMemcpyAsync(phi_d, phi, phi_size, cudaMemcpyHostToDevice, cuda_streams(stream_from_index(idx),device_id))
    ! cudaResult = cudaMemcpyAsync(umac_d, umac, umac_size, cudaMemcpyHostToDevice, cuda_streams(stream_from_index(idx),device_id))
    ! cudaResult = cudaMemcpyAsync(vmac_d, vmac, vmac_size, cudaMemcpyHostToDevice, cuda_streams(stream_from_index(idx),device_id))

    hdtdx = 0.5*(dt/dx)

    ! call slopex(glo, ghi, &
    call slopex_cuf(glo, ghi, &
                phi, ph_lo, ph_hi, &
                slope, glo, ghi &
#ifdef CUDA
                ,idx,device_id &
#endif
                )

    ! compute phi on x faces using umac to upwind; ignore transverse terms
    hdtdx_x = hdtdx(1)
    hdtdx_y = hdtdx(2)

    !$cuf kernel do(2) <<<(*,*), (16,16), 0, cuda_streams(stream_from_index(idx),device_id)>>> 
    do    j = lo(2)-1, hi(2)+1
       do i = lo(1)  , hi(1)+1
          if (umac(i,j) .lt. 0.d0) then
             phix_1d(i,j) = phi(i  ,j) - (0.5d0 + hdtdx_x*umac(i,j))*slope(i  ,j)
          else
             phix_1d(i,j) = phi(i-1,j) + (0.5d0 - hdtdx_x*umac(i,j))*slope(i-1,j)
          end if
       end do
    end do


    call slopey_cuf(glo, ghi, &
                phi, ph_lo, ph_hi, &
                slope, glo, ghi, idx, device_id)


    ! compute phi on y faces using umac to upwind; ignore transverse terms
    !$cuf kernel do(2) <<<(*,*), (16,16), 0, cuda_streams(stream_from_index(idx),device_id)>>> 
    do    j = lo(2)  , hi(2)+1
       do i = lo(1)-1, hi(1)+1

          if (vmac(i,j) .lt. 0.d0) then
             phiy_1d(i,j) = phi(i,j  ) - (0.5d0 + hdtdx_y*vmac(i,j))*slope(i,j  )
          else
             phiy_1d(i,j) = phi(i,j-1) + (0.5d0 - hdtdx_y*vmac(i,j))*slope(i,j-1)
          end if

       end do
    end do

    ! update phi on x faces by adding in y-transverse terms
    !$cuf kernel do(2) <<<(*,*), (16,16), 0, cuda_streams(stream_from_index(idx),device_id)>>> 
    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)+1

          if (umac(i,j) .lt. 0.d0) then
             phix(i,j) = phix_1d(i,j) &
                  - hdtdx_y*( 0.5d0*(vmac(i  ,j+1)+vmac(i  ,j)) * (phiy_1d(i  ,j+1)-phiy_1d(i  ,j)) )
          else
             phix(i,j) = phix_1d(i,j) &
                  - hdtdx_y*( 0.5d0*(vmac(i-1,j+1)+vmac(i-1,j)) * (phiy_1d(i-1,j+1)-phiy_1d(i-1,j)) )
          end if

          ! compute final x-fluxes
          flxx(i,j) = phix(i,j)*umac(i,j)

       end do
    end do


    ! update phi on y faces by adding in x-transverse terms
    !$cuf kernel do(2) <<<(*,*), (16,16), 0, cuda_streams(stream_from_index(idx),device_id)>>> 
    do    j = lo(2), hi(2)+1
       do i = lo(1), hi(1)

          if (vmac(i,j) .lt. 0.d0) then
             phiy(i,j) = phiy_1d(i,j) &
                  - hdtdx_x*( 0.5d0*(umac(i+1,j  )+umac(i,j  )) * (phix_1d(i+1,j  )-phix_1d(i,j  )) )
          else
             phiy(i,j) = phiy_1d(i,j) &
                  - hdtdx_x*( 0.5d0*(umac(i+1,j-1)+umac(i,j-1)) * (phix_1d(i+1,j-1)-phix_1d(i,j-1)) )
          end if

          ! compute final y-fluxes
          flxy(i,j) = phiy(i,j)*vmac(i,j)

       end do
    end do

    ! cudaResult = cudaMemcpyAsync(slope, slope_d, slope_size, cudaMemcpyDeviceToHost, cuda_streams(stream_from_index(idx),device_id))
    ! cudaResult = cudaMemcpyAsync(phix_1d, phix_1d_d, phix_1d_size, cudaMemcpyDeviceToHost, cuda_streams(stream_from_index(idx),device_id))
    ! cudaResult = cudaMemcpyAsync(phiy_1d, phiy_1d_d, phiy_1d_size, cudaMemcpyDeviceToHost, cuda_streams(stream_from_index(idx),device_id))
    ! cudaResult = cudaMemcpyAsync(phix, phix_d, phix_size, cudaMemcpyDeviceToHost, cuda_streams(stream_from_index(idx),device_id))
    ! cudaResult = cudaMemcpyAsync(phiy, phiy_d, phiy_size, cudaMemcpyDeviceToHost, cuda_streams(stream_from_index(idx),device_id))
    ! cudaResult = cudaMemcpyAsync(flxx, flxx_d, flxx_size, cudaMemcpyDeviceToHost, cuda_streams(stream_from_index(idx),device_id))
    ! cudaResult = cudaMemcpyAsync(flxy, flxy_d, flxy_size, cudaMemcpyDeviceToHost, cuda_streams(stream_from_index(idx),device_id))


    ! deallocate(phix_1d_d)
    ! deallocate(phiy_1d_d)
    ! deallocate(phix_d)
    ! deallocate(phiy_d)
    ! deallocate(umac_d)
    ! deallocate(vmac_d)
    ! deallocate(slope_d)
    ! deallocate(phi_d)
    ! deallocate(flxx_d)
    ! deallocate(flxy_d)

  end subroutine compute_flux_2d

end module compute_flux_module
