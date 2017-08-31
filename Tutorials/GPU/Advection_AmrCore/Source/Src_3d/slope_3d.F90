module slope_module
#ifdef CUDA
  use cuda_module, only: threads_and_blocks, stream_from_index
  use cuda_module, only: numThreads, numBlocks, cuda_streams 
  use iso_c_binding
#endif
 
  implicit none

  double precision, parameter:: four3rd=4.d0/3.d0, sixth=1.d0/6.d0
  
  private
 
  public :: slopex, slopey, slopez
 
contains
 
  subroutine slopex(lo, hi, &
                    q, qlo, qhi, &
                    dq, dqlo, dqhi &
#ifdef CUDA
                    , idx, device_id, tag &
#endif
                    )
#ifdef CUDA
    use mempool_module, only : gpu_allocate_hold, gpu_deallocate
#else
    use mempool_module, only : bl_allocate, bl_deallocate
#endif
    implicit none

    integer, intent(in) :: lo(3), hi(3), qlo(3), qhi(3), dqlo(3), dqhi(3)
    double precision, intent(in ) ::  q( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3))
    double precision, intent(out) :: dq(dqlo(1):dqhi(1),dqlo(2):dqhi(2),dqlo(3):dqhi(3))

    integer :: i, j, k
    ! double precision, dimension(lo(1)-1:hi(1)+1,lo(2):hi(2),lo(3):hi(3)) :: dsgn, dlim, df, dcen
    double precision :: dlft, drgt, dq1

#ifdef CUDA
    attributes(device) :: q, dq
    integer, intent(in) :: idx, device_id
    integer(kind=c_intptr_t), value, intent(in) :: tag
#endif

#ifdef CUDA
    double precision, dimension(:,:,:), pointer, contiguous, device :: dsgn, dlim, df, dcen
#else
    ! Some compiler may not support 'contiguous'.  Remove it in that case.
    double precision, dimension(:,:,:), pointer, contiguous :: dsgn, dlim, df, dcen
#endif

#ifdef CUDA
    call gpu_allocate_hold(dsgn  , tag, device_id, lo(1)-1, hi(1)+1, lo(2), hi(2), lo(3), hi(3)) 
    call gpu_allocate_hold(dlim  , tag, device_id, lo(1)-1, hi(1)+1, lo(2), hi(2), lo(3), hi(3)) 
    call gpu_allocate_hold(df    , tag, device_id, lo(1)-1, hi(1)+1, lo(2), hi(2), lo(3), hi(3)) 
    call gpu_allocate_hold(dcen  , tag, device_id, lo(1)-1, hi(1)+1, lo(2), hi(2), lo(3), hi(3)) 
#else
    call bl_allocate(dsgn, lo(1)-1, hi(1)+1, lo(2), hi(2), lo(3), hi(3))
    call bl_allocate(dlim, lo(1)-1, hi(1)+1, lo(2), hi(2), lo(3), hi(3))
    call bl_allocate(df  , lo(1)-1, hi(1)+1, lo(2), hi(2), lo(3), hi(3))
    call bl_allocate(dcen, lo(1)-1, hi(1)+1, lo(2), hi(2), lo(3), hi(3))
#endif

    ! TODO: move below to slopex_doit to remove pointerness of dsgn, dlim, df, dcen
    !$cuf kernel do(3) <<<*, (8,8,4), 0, cuda_streams(stream_from_index(idx),device_id)>>> 
    do    k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1)-1, hi(1)+1
             dlft = q(i  ,j,k) - q(i-1,j,k)
             drgt = q(i+1,j,k) - q(i  ,j,k)
             dcen(i,j,k) = .5d0 * (dlft+drgt)
             dsgn(i,j,k) = sign(1.d0, dcen(i,j,k))
             if (dlft*drgt .ge. 0.d0) then
                dlim(i,j,k) = 2.d0 * min( abs(dlft), abs(drgt) )
             else
                dlim(i,j,k) = 0.d0
             endif
             df(i,j,k) = dsgn(i,j,k)*min( dlim(i,j,k), abs(dcen(i,j,k)) )
          end do
       end do
    end do
          
    !$cuf kernel do(3) <<<*, (8,8,4), 0, cuda_streams(stream_from_index(idx),device_id)>>> 
    do    k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             dq1 = four3rd*dcen(i,j,k) - sixth*(df(i+1,j,k) + df(i-1,j,k))
             dq(i,j,k) = dsgn(i,j,k)*min(dlim(i,j,k),abs(dq1))
          end do
       end do
    end do

#ifdef CUDA
    ! call gpu_deallocate(dsgn, device_id)
    ! call gpu_deallocate(dlim, device_id)
    ! call gpu_deallocate(df  , device_id)
    ! call gpu_deallocate(dcen, device_id)
#else
    call bl_deallocate(dsgn)
    call bl_deallocate(dlim)
    call bl_deallocate(df)
    call bl_deallocate(dcen)
#endif

  end subroutine slopex


  subroutine slopey(lo, hi, &
                    q, qlo, qhi, &
                    dq, dqlo, dqhi&
#ifdef CUDA
                    , idx, device_id, tag &
#endif
                    )

#ifdef CUDA
    use mempool_module, only : gpu_allocate_hold, gpu_deallocate
#else
    use mempool_module, only : bl_allocate, bl_deallocate
#endif

    integer, intent(in) :: lo(3), hi(3), qlo(3), qhi(3), dqlo(3), dqhi(3)
    double precision, intent(in ) ::  q( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3))
    double precision, intent(out) :: dq(dqlo(1):dqhi(1),dqlo(2):dqhi(2),dqlo(3):dqhi(3))

#ifdef CUDA
    integer, intent(in) :: idx, device_id
    attributes(device) :: q, dq
    integer(kind=c_intptr_t), value, intent(in) :: tag
#endif

#ifdef CUDA
    double precision, dimension(:,:,:), pointer, contiguous, device :: dsgn, dlim, df, dcen
#else
    ! Some compiler may not support 'contiguous'.  Remove it in that case.
    double precision, dimension(:,:,:), pointer, contiguous :: dsgn, dlim, df, dcen
#endif

#ifdef CUDA
    call gpu_allocate_hold(dsgn, tag, device_id, lo(1), hi(1), lo(2)-1, hi(2)+1, lo(3), hi(3))
    call gpu_allocate_hold(dlim, tag, device_id, lo(1), hi(1), lo(2)-1, hi(2)+1, lo(3), hi(3))
    call gpu_allocate_hold(df  , tag, device_id, lo(1), hi(1), lo(2)-1, hi(2)+1, lo(3), hi(3))
    call gpu_allocate_hold(dcen, tag, device_id, lo(1), hi(1), lo(2)-1, hi(2)+1, lo(3), hi(3))

#else
    call bl_allocate(dsgn, lo(1), hi(1), lo(2)-1, hi(2)+1, lo(3), hi(3))
    call bl_allocate(dlim, lo(1), hi(1), lo(2)-1, hi(2)+1, lo(3), hi(3))
    call bl_allocate(df  , lo(1), hi(1), lo(2)-1, hi(2)+1, lo(3), hi(3))
    call bl_allocate(dcen, lo(1), hi(1), lo(2)-1, hi(2)+1, lo(3), hi(3))
#endif

    call slopey_doit(lo, hi, &
                     q, qlo, qhi, &
                     dq, dqlo, dqhi, &
                     dsgn, dlim, df, dcen, (/lo(1),lo(2)-1,lo(3)/), (/hi(1),hi(2)+1,hi(3)/) &
#ifdef CUDA
                    , idx, device_id &
#endif
                     )

#ifdef CUDA
    ! call gpu_deallocate(dsgn, device_id)
    ! call gpu_deallocate(dlim, device_id)
    ! call gpu_deallocate(df  , device_id)
    ! call gpu_deallocate(dcen, device_id)
#else
    call bl_deallocate(dsgn)
    call bl_deallocate(dlim)
    call bl_deallocate(df)
    call bl_deallocate(dcen)
#endif

  end subroutine slopey

  subroutine slopey_doit(lo, hi, &
                         q, qlo, qhi, &
                         dq, dqlo, dqhi, &
                         dsgn, dlim, df, dcen, ddlo, ddhi &
#ifdef CUDA
                    , idx, device_id &
#endif
                         )

    integer, intent(in) :: lo(3), hi(3), qlo(3), qhi(3), dqlo(3), dqhi(3), &
         ddlo(3), ddhi(3)
    double precision, intent(in ) ::  q  ( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3))
    double precision, intent(out) :: dq  (dqlo(1):dqhi(1),dqlo(2):dqhi(2),dqlo(3):dqhi(3))
    double precision              :: dsgn(ddlo(1):ddhi(1),ddlo(2):ddhi(2),ddlo(3):ddhi(3))
    double precision              :: dlim(ddlo(1):ddhi(1),ddlo(2):ddhi(2),ddlo(3):ddhi(3))
    double precision              :: df  (ddlo(1):ddhi(1),ddlo(2):ddhi(2),ddlo(3):ddhi(3))
    double precision              :: dcen(ddlo(1):ddhi(1),ddlo(2):ddhi(2),ddlo(3):ddhi(3))

    integer :: i, j, k
    double precision :: dlft, drgt, dq1
#ifdef CUDA
    integer, intent(in) :: idx, device_id
    attributes(device) :: q, dq, dsgn, dlim, df, dcen
#endif

    !$cuf kernel do(3) <<<*, (8,8,4), 0, cuda_streams(stream_from_index(idx),device_id)>>> 
    do k = lo(3), hi(3)
       do j    = lo(2)-1, hi(2)+1
          do i = lo(1)  , hi(1)
             dlft = q(i,j  ,k) - q(i,j-1,k)
             drgt = q(i,j+1,k) - q(i,j  ,k)
             dcen(i,j,k) = .5d0 * (dlft+drgt)
             dsgn(i,j,k) = sign( 1.d0, dcen(i,j,k) )
             if (dlft*drgt .ge. 0.d0) then
                dlim(i,j,k) = 2.d0 * min( abs(dlft), abs(drgt) )
             else
                dlim(i,j,k) = 0.d0
             endif
             df(i,j,k) = dsgn(i,j,k)*min( dlim(i,j,k),abs(dcen(i,j,k)) )
          end do
       end do
    end do
       
    !$cuf kernel do(3) <<<*, (8,8,4), 0, cuda_streams(stream_from_index(idx),device_id)>>> 
    do k = lo(3), hi(3)
       do j    = lo(2), hi(2)
          do i = lo(1), hi(1)
             dq1 = four3rd*dcen(i,j,k) - sixth*( df(i,j+1,k) + df(i,j-1,k) )
             dq(i,j,k) = dsgn(i,j,k)*min(dlim(i,j,k),abs(dq1))
          end do
       end do
    end do

  end subroutine slopey_doit


  subroutine slopez(lo, hi, &
                    q, qlo, qhi, &
                    dq, dqlo, dqhi &
#ifdef CUDA
                    , idx, device_id, tag &
#endif
                    )

#ifdef CUDA
    use mempool_module, only : gpu_allocate_hold, gpu_deallocate
#else
    use mempool_module, only : bl_allocate, bl_deallocate
#endif

    integer, intent(in) :: lo(3), hi(3), qlo(3), qhi(3), dqlo(3), dqhi(3)
    double precision, intent(in ) ::  q( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3))
    double precision, intent(out) :: dq(dqlo(1):dqhi(1),dqlo(2):dqhi(2),dqlo(3):dqhi(3))

#ifdef CUDA
    integer, intent(in) :: idx, device_id
    attributes(device) :: q, dq
    integer(kind=c_intptr_t), value, intent(in) :: tag
#endif

#ifdef CUDA
    double precision, dimension(:,:,:), pointer, contiguous, device :: dsgn, dlim, df, dcen
#else

    ! Some compiler may not support 'contiguous'.  Remove it in that case.
    double precision, dimension(:,:,:), pointer, contiguous :: dsgn, dlim, df, dcen
#endif

#ifdef CUDA
    call gpu_allocate_hold(dsgn, tag, device_id, lo(1), hi(1), lo(2), hi(2), lo(3)-1, hi(3)+1)
    call gpu_allocate_hold(dlim, tag, device_id, lo(1), hi(1), lo(2), hi(2), lo(3)-1, hi(3)+1)
    call gpu_allocate_hold(df  , tag, device_id, lo(1), hi(1), lo(2), hi(2), lo(3)-1, hi(3)+1)
    call gpu_allocate_hold(dcen, tag, device_id, lo(1), hi(1), lo(2), hi(2), lo(3)-1, hi(3)+1)
#else
    call bl_allocate(dsgn, lo(1), hi(1), lo(2), hi(2), lo(3)-1, hi(3)+1)
    call bl_allocate(dlim, lo(1), hi(1), lo(2), hi(2), lo(3)-1, hi(3)+1)
    call bl_allocate(df  , lo(1), hi(1), lo(2), hi(2), lo(3)-1, hi(3)+1)
    call bl_allocate(dcen, lo(1), hi(1), lo(2), hi(2), lo(3)-1, hi(3)+1)
#endif

    call slopez_doit(lo, hi, &
                     q, qlo, qhi, &
                     dq, dqlo, dqhi, &
                     dsgn, dlim, df, dcen, &
                     (/lo(1),lo(2),lo(3)-1/), (/hi(1),hi(2),hi(3)+1/) &
#ifdef CUDA
                    , idx, device_id &
#endif
                     )

#ifdef CUDA
    ! call gpu_deallocate(dsgn, device_id)
    ! call gpu_deallocate(dlim, device_id)
    ! call gpu_deallocate(df  , device_id)
    ! call gpu_deallocate(dcen, device_id)
#else
    call bl_deallocate(dsgn)
    call bl_deallocate(dlim)
    call bl_deallocate(df)
    call bl_deallocate(dcen)
#endif

  end subroutine slopez

  subroutine slopez_doit(lo, hi, &
                         q, qlo, qhi, &
                         dq, dqlo, dqhi, &
                         dsgn, dlim, df, dcen, ddlo, ddhi &
#ifdef CUDA
                    , idx, device_id &
#endif
                         )

    integer, intent(in) :: lo(3), hi(3), qlo(3), qhi(3), dqlo(3), dqhi(3), &
         ddlo(3), ddhi(3)
    double precision, intent(in ) ::  q  ( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3))
    double precision, intent(out) :: dq  (dqlo(1):dqhi(1),dqlo(2):dqhi(2),dqlo(3):dqhi(3))
    double precision              :: dsgn(ddlo(1):ddhi(1),ddlo(2):ddhi(2),ddlo(3):ddhi(3))
    double precision              :: dlim(ddlo(1):ddhi(1),ddlo(2):ddhi(2),ddlo(3):ddhi(3))
    double precision              :: df  (ddlo(1):ddhi(1),ddlo(2):ddhi(2),ddlo(3):ddhi(3))
    double precision              :: dcen(ddlo(1):ddhi(1),ddlo(2):ddhi(2),ddlo(3):ddhi(3))

    integer :: i, j, k
    double precision :: dlft, drgt, dq1
#ifdef CUDA
    integer, intent(in) :: idx, device_id
    attributes(device) :: q, dq, dsgn, dlim, df, dcen
#endif


    ! first compute Fromm slopes
    !$cuf kernel do(3) <<<*, (8,8,4), 0, cuda_streams(stream_from_index(idx),device_id)>>> 
    do k       = lo(3)-1, hi(3)+1
       do j    = lo(2)  , hi(2)
          do i = lo(1)  , hi(1)
             dlft = q(i,j,k  ) - q(i,j,k-1)
             drgt = q(i,j,k+1) - q(i,j,k  )
             dcen(i,j,k) = .5d0 * (dlft+drgt)
             dsgn(i,j,k) = sign( 1.d0, dcen(i,j,k) )
             if (dlft*drgt .ge. 0.d0) then
                dlim(i,j,k) = 2.d0 * min( abs(dlft), abs(drgt) )
             else
                dlim(i,j,k) = 0.d0
             endif
             df(i,j,k) = dsgn(i,j,k)*min( dlim(i,j,k),abs(dcen(i,j,k)) )
          end do
       end do
    end do
       
    !$cuf kernel do(3) <<<*, (8,8,4), 0, cuda_streams(stream_from_index(idx),device_id)>>> 
    ! Now compute limited fourth order slopes
    do k       = lo(3), hi(3)
       do j    = lo(2), hi(2)
          do i = lo(1), hi(1)
             dq1 = four3rd*dcen(i,j,k) - sixth*( df(i,j,k+1) + df(i,j,k-1) )
             dq(i,j,k) = dsgn(i,j,k)*min(dlim(i,j,k),abs(dq1))
          end do
       end do
    end do

  end subroutine slopez_doit

end module slope_module 
