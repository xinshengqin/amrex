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
                             phix_1d, phiy_1d, phix, phiy, slope, glo, ghi)

    use slope_module, only: slopey
    use slope_module_cuda, only: slopex_d
    use cuda_module, only: threads_and_blocks, stream_from_index
    use cuda_module, only: numThreads, numBlocks, cuda_streams 
    ! use slope_module, only: slopex

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

    double precision :: hdtdx_x, hdtdx_y

    ! device 
    double precision, allocatable, device :: phix_1d_d(:,:), phi_d(:,:), slope_d(:,:), umac_d(:,:)
    allocate(phix_1d_d(glo(1):ghi(1),glo(2):ghi(2)))
    allocate(slope_d(glo(1):ghi(1),glo(2):ghi(2)))
    allocate(phi_d(ph_lo(1):ph_hi(1),ph_lo(2):ph_hi(2)))
    allocate(umac_d( u_lo(1): u_hi(1), u_lo(2): u_hi(2)))


    phi_d = phi
    umac_d = umac

    hdtdx = 0.5*(dt/dx)

    ! call slopex(glo, ghi, &
    !             phi, ph_lo, ph_hi, &
    !             slope, glo, ghi)
    call slopex_d(glo, ghi, &
                phi_d, ph_lo, ph_hi, &
                slope_d, glo, ghi,0,0)
    slope = slope_d

    ! compute phi on x faces using umac to upwind; ignore transverse terms
    ! do    j = lo(2)-1, hi(2)+1
    !    do i = lo(1)  , hi(1)+1

    !       if (umac(i,j) .lt. 0.d0) then
    !          phix_1d(i,j) = phi(i  ,j) - (0.5d0 + hdtdx(1)*umac(i,j))*slope(i  ,j)
    !       else
    !          phix_1d(i,j) = phi(i-1,j) + (0.5d0 - hdtdx(1)*umac(i,j))*slope(i-1,j)
    !       end if

    !    end do
    ! end do

    hdtdx_x = hdtdx(1)
    hdtdx_y = hdtdx(2)
    !$cuf kernel do(2) <<< (*,*), (16,16) >>>
    do    j = lo(2)-1, hi(2)+1
       do i = lo(1)  , hi(1)+1
          if (umac_d(i,j) .lt. 0.d0) then
             phix_1d_d(i,j) = phi_d(i  ,j) - (0.5d0 + hdtdx_x*umac_d(i,j))*slope_d(i  ,j)
          else
             phix_1d_d(i,j) = phi_d(i-1,j) + (0.5d0 - hdtdx_x*umac_d(i,j))*slope_d(i-1,j)
          end if
       end do
    end do
    phix_1d = phix_1d_d


    ! kernel
    ! call threads_and_blocks([lo(1), lo(2)-1], [hi(1)+1, hi(2)+1], numBlocks, numThreads)
    ! call sum_up_kernel<<<numBlocks, numThreads, 0 >>> &
    !     (lo(1), lo(2)-1, hi(1)+1, hi(2)+1, &
    !     phi_d,ph_lo(1),ph_hi(1),ph_lo(2),ph_hi(2), &
    !     umac_d,u_lo(1), u_hi(1), u_lo(2), u_hi(2), &
    !     slope_d,glo(1),ghi(1),glo(2),ghi(2), &
    !     phix_1d_d,glo(1),ghi(1),glo(2),ghi(2), &
    !     hdtdx_x, hdtdx_y)

    ! phix_1d = phix_1d_d

    deallocate(phix_1d_d)
    deallocate(slope_d)
    deallocate(phi_d)
    deallocate(umac_d)

    call slopey(glo, ghi, &
                phi, ph_lo, ph_hi, &
                slope, glo, ghi)

    ! compute phi on y faces using umac to upwind; ignore transverse terms
    do    j = lo(2)  , hi(2)+1
       do i = lo(1)-1, hi(1)+1

          if (vmac(i,j) .lt. 0.d0) then
             phiy_1d(i,j) = phi(i,j  ) - (0.5d0 + hdtdx(2)*vmac(i,j))*slope(i,j  )
          else
             phiy_1d(i,j) = phi(i,j-1) + (0.5d0 - hdtdx(2)*vmac(i,j))*slope(i,j-1)
          end if

       end do
    end do

    ! update phi on x faces by adding in y-transverse terms
    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)+1

          if (umac(i,j) .lt. 0.d0) then
             phix(i,j) = phix_1d(i,j) &
                  - hdtdx(2)*( 0.5d0*(vmac(i  ,j+1)+vmac(i  ,j)) * (phiy_1d(i  ,j+1)-phiy_1d(i  ,j)) )
          else
             phix(i,j) = phix_1d(i,j) &
                  - hdtdx(2)*( 0.5d0*(vmac(i-1,j+1)+vmac(i-1,j)) * (phiy_1d(i-1,j+1)-phiy_1d(i-1,j)) )
          end if

          ! compute final x-fluxes
          flxx(i,j) = phix(i,j)*umac(i,j)

       end do
    end do

    ! update phi on y faces by adding in x-transverse terms
    do    j = lo(2), hi(2)+1
       do i = lo(1), hi(1)

          if (vmac(i,j) .lt. 0.d0) then
             phiy(i,j) = phiy_1d(i,j) &
                  - hdtdx(1)*( 0.5d0*(umac(i+1,j  )+umac(i,j  )) * (phix_1d(i+1,j  )-phix_1d(i,j  )) )
          else
             phiy(i,j) = phiy_1d(i,j) &
                  - hdtdx(1)*( 0.5d0*(umac(i+1,j-1)+umac(i,j-1)) * (phix_1d(i+1,j-1)-phix_1d(i,j-1)) )
          end if

          ! compute final y-fluxes
          flxy(i,j) = phiy(i,j)*vmac(i,j)

       end do
    end do

  end subroutine compute_flux_2d

    attributes(global) &
    subroutine sum_up_kernel(lo_x, lo_y, hi_x, hi_y, &
                    phi, phi_lo_x, phi_lo_y, phi_hi_x, phi_hi_y, &
                    umac, umac_lo_x, umac_lo_y, umac_hi_x, umac_hi_y, &
                    slope, slope_lo_x, slope_lo_y, slope_hi_x, slope_hi_y, &
                    phix_1d, phix_1d_lo_x, phix_1d_lo_y, phix_1d_hi_x, phix_1d_hi_y, &
                    hdtdx_x, hdtdx_y) 


        use amrex_fort_module, only: get_loop_bounds
        implicit none

        integer, value, intent(in) :: lo_x, lo_y, hi_x, hi_y
        integer, value, intent(in) :: phi_lo_x, phi_lo_y, phi_hi_x, phi_hi_y
        integer, value, intent(in) :: umac_lo_x, umac_lo_y, umac_hi_x, umac_hi_y
        integer, value, intent(in) :: slope_lo_x, slope_lo_y, slope_hi_x, slope_hi_y
        integer, value, intent(in) :: phix_1d_lo_x, phix_1d_lo_y, phix_1d_hi_x, phix_1d_hi_y

        double precision, intent(in ) :: phi(phi_lo_x:phi_hi_x, phi_lo_y:phi_hi_y)
        double precision, intent(in ) :: umac(umac_lo_x:umac_hi_x, umac_lo_y:umac_hi_y)
        double precision, intent(in ) :: slope(slope_lo_x:slope_hi_x, slope_lo_y:slope_hi_y)
        double precision, intent(out) :: phix_1d(phix_1d_lo_x:phix_1d_hi_x, phix_1d_lo_y:phix_1d_hi_y)
        double precision, value, intent(in)  :: hdtdx_x, hdtdx_y

        integer :: blo(2), bhi(2)
        integer :: i, j

        call get_loop_bounds(blo, bhi, [lo_x, lo_y], [hi_x, hi_y])
        ! compute phi on x faces using umac to upwind; ignore transverse terms
        do    j = blo(2), bhi(2)
           do i = blo(1), bhi(1)

              if (umac(i,j) .lt. 0.d0) then
                 phix_1d(i,j) = phi(i  ,j) - (0.5d0 + hdtdx_x*umac(i,j))*slope(i  ,j)
              else
                 phix_1d(i,j) = phi(i-1,j) + (0.5d0 - hdtdx_x*umac(i,j))*slope(i-1,j)
              end if

           end do
        end do
    end subroutine sum_up_kernel

end module compute_flux_module
