module slope_module_cuda
  use cuda_module, only: threads_and_blocks, stream_from_index
  use cuda_module, only: numThreads, numBlocks, cuda_streams 
  use cudafor, only: cudaMemset
 
  implicit none

  double precision, parameter:: four3rd=4.d0/3.d0, sixth=1.d0/6.d0
  
  private
 
  public :: slopex, slopey_cuf, slopex_cuf
 
contains
  subroutine slopex_cuf(lo, hi, &
                    q, qlo, qhi, &
                    dq, dqlo, dqhi &
                    , idx, device_id)

    implicit none

    integer, intent(in) :: lo(2), hi(2), qlo(2), qhi(2), dqlo(2), dqhi(2)
    ! double precision, intent(in ) ::  q( qlo(1): qhi(1), qlo(2): qhi(2))
    ! double precision, intent(out) :: dq(dqlo(1):dqhi(1),dqlo(2):dqhi(2))
    double precision, device, intent(in ) ::  q( qlo(1): qhi(1), qlo(2): qhi(2))
    double precision, device, intent(out) :: dq(dqlo(1):dqhi(1),dqlo(2):dqhi(2))

    integer :: i, j, cudaResult
    double precision, device, allocatable :: dsgn(:), dlim(:), df(:), dcen(:)
    double precision  :: dlft, drgt, dq1
    ! TODO: if I add device attribute to these scalars, the result is incorrect
    ! double precision, device  :: dlft, drgt, dq1
    integer, intent(in) :: idx, device_id

    ! TODO: can replace this with customized memory allocator
    allocate(dsgn(lo(1)-1:hi(1)+1)) 
    allocate(dlim(lo(1)-1:hi(1)+1))
    allocate(df(lo(1)-1:hi(1)+1))
    allocate(dcen(lo(1)-1:hi(1)+1)) 

    do j = lo(2), hi(2)
       ! first compute Fromm slopes
       !$cuf kernel do <<<*, 128, 0, cuda_streams(stream_from_index(idx),device_id)>>> 
       do i = lo(1)-1, hi(1)+1
          dlft = q(i  ,j) - q(i-1,j)
          drgt = q(i+1,j) - q(i  ,j)
          dcen(i) = .5d0 * (dlft+drgt)
          dsgn(i) = sign(1.d0, dcen(i))
          if (dlft*drgt .ge. 0.d0) then
             dlim(i) = 2.d0 * min( abs(dlft), abs(drgt) )
          else
             dlim(i) = 0.d0
          endif
          df(i) = dsgn(i)*min( dlim(i), abs(dcen(i)) )
       end do

       ! Now limited fourth order slopes
       !$cuf kernel do <<<*, 128, 0, cuda_streams(stream_from_index(idx),device_id)>>> 
       do i = lo(1), hi(1)
          dq1 = four3rd*dcen(i) - sixth*(df(i+1) + df(i-1))
          dq(i,j) = dsgn(i)*min(dlim(i),abs(dq1))
       end do
    enddo

    deallocate(dsgn)
    deallocate(dlim)
    deallocate(df)
    deallocate(dcen)

  end subroutine slopex_cuf
 
  subroutine slopex(lo, hi, &
                    q_d, qlo, qhi, &
                    dq_d, dqlo, dqhi, idx, device_id)



    implicit none

    integer, intent(in) :: lo(2), hi(2), qlo(2), qhi(2), dqlo(2), dqhi(2)
    double precision, device, intent(in ) ::  q_d( qlo(1): qhi(1), qlo(2): qhi(2))
    double precision, device, intent(out) :: dq_d(dqlo(1):dqhi(1),dqlo(2):dqhi(2))


    integer, intent(in) :: idx, device_id


    call threads_and_blocks(lo, hi, numBlocks, numThreads)
    call slopex_kernel<<<numBlocks, numThreads, 0, cuda_streams(stream_from_index(idx),device_id)>>> &
        (lo(1), lo(2), hi(1), hi(2), &
         q_d, qlo(1), qlo(2), qhi(1), qhi(2), &
         dq_d, dqlo(1), dqlo(2), dqhi(1), dqhi(2))

  end subroutine slopex


  subroutine slopey_cuf(lo, hi, &
                    q, qlo, qhi, &
                    dq, dqlo, dqhi, idx, device_id)

    use mempool_module, only : bl_allocate, bl_deallocate

    integer, intent(in) :: lo(2), hi(2), qlo(2), qhi(2), dqlo(2), dqhi(2)
    double precision, device, intent(in ) ::  q( qlo(1): qhi(1), qlo(2): qhi(2))
    double precision, device, intent(out) :: dq(dqlo(1):dqhi(1),dqlo(2):dqhi(2))

    double precision, device, allocatable :: dsgn(:,:), dlim(:,:), df(:,:), dcen(:,:)
    
    integer, intent(in) :: idx, device_id

#ifdef CUDA
    attributes(device) :: q, dq
#endif

    ! ! Some compiler may not support 'contiguous'.  Remove it in that case.
    ! double precision, dimension(:,:), pointer, contiguous :: dsgn, dlim, df, dcen

    ! call bl_allocate(dsgn, lo(1), hi(1), lo(2)-1, hi(2)+1)
    ! call bl_allocate(dlim, lo(1), hi(1), lo(2)-1, hi(2)+1)
    ! call bl_allocate(df  , lo(1), hi(1), lo(2)-1, hi(2)+1)
    ! call bl_allocate(dcen, lo(1), hi(1), lo(2)-1, hi(2)+1)

    allocate(dsgn(lo(1):hi(1), lo(2)-1:hi(2)+1))
    allocate(dlim(lo(1):hi(1), lo(2)-1:hi(2)+1))
    allocate(df  (lo(1):hi(1), lo(2)-1:hi(2)+1))
    allocate(dcen(lo(1):hi(1), lo(2)-1:hi(2)+1))

    call slopey_doit(lo, hi, &
                     q, qlo, qhi, &
                     dq, dqlo, dqhi, &
                     dsgn, dlim, df, dcen, (/lo(1),lo(2)-1/), (/hi(1),hi(2)+1/), idx, device_id)

    ! call bl_deallocate(dsgn)
    ! call bl_deallocate(dlim)
    ! call bl_deallocate(df)
    ! call bl_deallocate(dcen)
    ! deallocate(dsgn)
    ! deallocate(dlim)
    ! deallocate(df)
    ! deallocate(dcen)

  end subroutine slopey_cuf

  subroutine slopey_doit(lo, hi, &
                         q, qlo, qhi, &
                         dq, dqlo, dqhi, &
                         dsgn, dlim, df, dcen, ddlo, ddhi &
                         , idx, device_id &
                         )

    integer, intent(in) :: lo(2), hi(2), qlo(2), qhi(2), dqlo(2), dqhi(2), &
         ddlo(2), ddhi(2)
    integer, intent(in) :: idx, device_id
    double precision, device, intent(in ) ::  q  ( qlo(1): qhi(1), qlo(2): qhi(2))
    double precision, device, intent(out) :: dq  (dqlo(1):dqhi(1),dqlo(2):dqhi(2))
    double precision, device              :: dsgn(ddlo(1):ddhi(1),ddlo(2):ddhi(2))
    double precision, device              :: dlim(ddlo(1):ddhi(1),ddlo(2):ddhi(2))
    double precision, device              :: df  (ddlo(1):ddhi(1),ddlo(2):ddhi(2))
    double precision, device              :: dcen(ddlo(1):ddhi(1),ddlo(2):ddhi(2))

    integer :: i, j 
    double precision :: dlft, drgt, dq1


    ! first compute Fromm slopes
    !$cuf kernel do(2) <<<(*,*), (16,16), 0, cuda_streams(stream_from_index(idx),device_id)>>> 
    do j    = lo(2)-1, hi(2)+1
       do i = lo(1)  , hi(1)
          dlft = q(i,j  ) - q(i,j-1)
          drgt = q(i,j+1) - q(i,j  )
          dcen(i,j) = .5d0 * (dlft+drgt)
          dsgn(i,j) = sign( 1.d0, dcen(i,j) )
          if (dlft*drgt .ge. 0.d0) then
             dlim(i,j) = 2.d0 * min( abs(dlft), abs(drgt) )
          else
             dlim(i,j) = 0.d0
          endif
          df(i,j) = dsgn(i,j)*min( dlim(i,j),abs(dcen(i,j)) )
       end do
    end do

    ! Now compute limited fourth order slopes
    !$cuf kernel do(2) <<<(*,*), (16,16), 0, cuda_streams(stream_from_index(idx),device_id)>>> 
    do j    = lo(2), hi(2)
       do i = lo(1), hi(1)
          dq1 = four3rd*dcen(i,j) - sixth*( df(i,j+1) + df(i,j-1) )
          dq(i,j) = dsgn(i,j)*min(dlim(i,j),abs(dq1))
       end do
    end do

  end subroutine slopey_doit

    attributes(global) &
    subroutine slopex_kernel(lo_x, lo_y, hi_x, hi_y, &
                    q, qlo_x, qlo_y, qhi_x, qhi_y, &
                    dq, dqlo_x, dqlo_y, dqhi_x, dqhi_y)

    use amrex_fort_module, only: get_loop_bounds
    implicit none

    integer, value, intent(in) :: lo_x, lo_y, hi_x, hi_y
    integer, value, intent(in) :: qlo_x, qlo_y, qhi_x, qhi_y
    integer, value, intent(in) :: dqlo_x, dqlo_y, dqhi_x, dqhi_y
    double precision, intent(in ) ::  q( qlo_x: qhi_x, qlo_y: qhi_y)
    double precision, intent(out) :: dq(dqlo_x:dqhi_x,dqlo_y:dqhi_y)

    integer :: blo(2), bhi(2)

    call get_loop_bounds(blo, bhi, [lo_x, lo_y], [hi_x, hi_y])
    call slopex_doit(blo, bhi, &
                     q,  [qlo_x, qlo_y] , [qhi_x, qhi_y], &
                     dq, [dqlo_x, dqlo_y] , [dqhi_x, dqhi_y])
    end subroutine slopex_kernel

    ! this is the work for each thread 
    ! lo == hi == space index of the cell this thread should work on
    attributes(device) &
    subroutine slopex_doit(lo, hi, &
                    q, qlo, qhi, &
                    dq, dqlo, dqhi)

    implicit none

    integer, intent(in) :: lo(2), hi(2), qlo(2), qhi(2), dqlo(2), dqhi(2)
    double precision, intent(in ) ::  q( qlo(1): qhi(1), qlo(2): qhi(2))
    double precision, intent(out) :: dq(dqlo(1):dqhi(1),dqlo(2):dqhi(2))

    integer :: i, j
    double precision, dimension(lo(1)-1:hi(1)+1) :: dsgn, dlim, df, dcen
    double precision :: dlft, drgt, dq1

    do j = lo(2), hi(2)

       ! first compute Fromm slopes
       do i = lo(1)-1, hi(1)+1
          dlft = q(i  ,j) - q(i-1,j)
          drgt = q(i+1,j) - q(i  ,j)
          dcen(i) = .5d0 * (dlft+drgt)
          dsgn(i) = sign(1.d0, dcen(i))
          if (dlft*drgt .ge. 0.d0) then
             dlim(i) = 2.d0 * min( abs(dlft), abs(drgt) )
          else
             dlim(i) = 0.d0
          endif
          df(i) = dsgn(i)*min( dlim(i), abs(dcen(i)) )
       end do

       ! Now limited fourth order slopes
       do i = lo(1), hi(1)
          dq1 = four3rd*dcen(i) - sixth*(df(i+1) + df(i-1))
          dq(i,j) = dsgn(i)*min(dlim(i),abs(dq1))
       end do
    enddo
    end subroutine slopex_doit

end module slope_module_cuda
