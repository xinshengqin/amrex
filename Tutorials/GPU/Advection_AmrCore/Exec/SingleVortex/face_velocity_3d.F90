module get_face_velocity_module
#ifdef CUDA
  use iso_c_binding
  use cuda_module, only: threads_and_blocks, stream_from_index
  use cuda_module, only: numThreads, numBlocks, cuda_streams 
#endif
  implicit none
  contains

  subroutine get_face_velocity(level, time, &
       vx, vx_l1, vx_l2, vx_l3, vx_h1, vx_h2, vx_h3, &
       vy, vy_l1, vy_l2, vy_l3, vy_h1, vy_h2, vy_h3, &
       vz, vz_l1, vz_l2, vz_l3, vz_h1, vz_h2, vz_h3, &
       dx, prob_lo &
  #ifdef CUDA
       , idx, device_id, tag &
  #endif
       ) bind(C, name="get_face_velocity")
  
    implicit none
  
    integer, intent(in) :: level
    double precision, intent(in) :: time
    integer, intent(in) :: vx_l1, vx_l2, vx_l3, vx_h1, vx_h2, vx_h3
    integer, intent(in) :: vy_l1, vy_l2, vy_l3, vy_h1, vy_h2, vy_h3
    integer, intent(in) :: vz_l1, vz_l2, vz_l3, vz_h1, vz_h2, vz_h3
    double precision, intent(out) :: vx(vx_l1:vx_h1,vx_l2:vx_h2,vx_l3:vx_h3)
    double precision, intent(out) :: vy(vy_l1:vy_h1,vy_l2:vy_h2,vy_l3:vy_h3)
    double precision, intent(out) :: vz(vz_l1:vz_h1,vz_l2:vz_h2,vz_l3:vz_h3)
    double precision, intent(in) :: dx(3), prob_lo(3)
  
    integer :: i, j, k, plo(2), phi(2)
    double precision :: x, y, z
  #ifdef CUDA
      attributes(device) :: vx, vy, vz
      integer, value, intent(in) :: idx, device_id
      integer(kind=c_intptr_t), value, intent(in) :: tag
  #endif
    call get_face_velocity_doit(level, time, &
       vx, vx_l1, vx_l2, vx_l3, vx_h1, vx_h2, vx_h3, &
       vy, vy_l1, vy_l2, vy_l3, vy_h1, vy_h2, vy_h3, &
       vz, vz_l1, vz_l2, vz_l3, vz_h1, vz_h2, vz_h3, &
       dx, prob_lo &
  #ifdef CUDA
       , idx, device_id, tag &
  #endif
       )
  
  end subroutine get_face_velocity

  subroutine get_face_velocity_doit(level, time, &
       vx, vx_l1, vx_l2, vx_l3, vx_h1, vx_h2, vx_h3, &
       vy, vy_l1, vy_l2, vy_l3, vy_h1, vy_h2, vy_h3, &
       vz, vz_l1, vz_l2, vz_l3, vz_h1, vz_h2, vz_h3, &
       dx, prob_lo &
  #ifdef CUDA
       , idx, device_id, tag &
  #endif
       ) 
  
  #ifdef CUDA
      use mempool_module, only : gpu_allocate_hold
  #else
      use mempool_module, only : bl_allocate, bl_deallocate
  #endif
  
    implicit none
  
    integer, intent(in) :: level
    double precision, intent(in) :: time
    integer, intent(in) :: vx_l1, vx_l2, vx_l3, vx_h1, vx_h2, vx_h3
    integer, intent(in) :: vy_l1, vy_l2, vy_l3, vy_h1, vy_h2, vy_h3
    integer, intent(in) :: vz_l1, vz_l2, vz_l3, vz_h1, vz_h2, vz_h3
    double precision, intent(out) :: vx(vx_l1:vx_h1,vx_l2:vx_h2,vx_l3:vx_h3)
    double precision, intent(out) :: vy(vy_l1:vy_h1,vy_l2:vy_h2,vy_l3:vy_h3)
    double precision, intent(out) :: vz(vz_l1:vz_h1,vz_l2:vz_h2,vz_l3:vz_h3)
    double precision, intent(in) :: dx(3), prob_lo(3)
  
    integer :: i, j, k, plo(2), phi(2)
    double precision :: x, y, z
    double precision :: prob_lo_x, prob_lo_y, prob_lo_z, dx_x, dx_y, dx_z
  #ifdef CUDA
    double precision, dimension(:,:), pointer, contiguous, device :: psi
  #else
    double precision, pointer, contiguous :: psi(:,:)
  #endif
    double precision, parameter :: M_PI = 3.141592653589793238462643383279502884197d0
  
  #ifdef CUDA
      attributes(device) :: vx, vy, vz
      integer, value, intent(in) :: idx, device_id
      integer(kind=c_intptr_t), value, intent(in) :: tag
  #endif
  
  
    plo(1) = min(vx_l1-1, vy_l1-1)
    plo(2) = min(vx_l2-1, vy_l2-1)
    phi(1) = max(vx_h1  , vy_h1+1)
    phi(2) = max(vx_h2+1, vy_h2  )

    prob_lo_x = prob_lo(1)
    prob_lo_y = prob_lo(2)
    prob_lo_z = prob_lo(3)
    dx_x = dx(1)
    dx_y = dx(2)
    
  #ifdef CUDA
    call gpu_allocate_hold(psi , tag, device_id, plo(1), phi(1), plo(2), phi(2)) 
  #else
    call bl_allocate(psi, plo(1), phi(1), plo(2), phi(2))
  #endif
  
    ! streamfunction psi
    !$cuf kernel do(2) <<<*, (16,16), 0, cuda_streams(stream_from_index(idx),device_id)>>> 
    do j = plo(2), phi(2)
       do i = plo(1), phi(1)
          y = (dble(j)+0.5d0)*dx_y + prob_lo_y
          x = (dble(i)+0.5d0)*dx_x + prob_lo_x
          psi(i,j) =  sin(M_PI*x)**2 * sin(M_PI*y)**2 * cos (M_PI*time/2.d0) * (1.d0 / M_PI)
       end do
    end do
    
    ! x velocity
    !$cuf kernel do(3) <<<*, (8,8,4), 0, cuda_streams(stream_from_index(idx),device_id)>>> 
    do k = vx_l3, vx_h3
      do j = vx_l2, vx_h2
        do i = vx_l1, vx_h1
          y = (dble(j)+0.5d0) * dx_y + prob_lo_y
          x = dble(i) * dx_x + prob_lo_x
          vx(i,j,k) =  -( (psi(i,j+1)+psi(i-1,j+1)) - (psi(i,j-1)+psi(i-1,j-1)) ) * (0.25d0/dx_y)
        end do
      end do
    end do
  
    ! y velocity
    !$cuf kernel do(3) <<<*, (8,8,4), 0, cuda_streams(stream_from_index(idx),device_id)>>> 
    do k = vy_l3, vy_h3
      do j = vy_l2, vy_h2
        do i = vy_l1, vy_h1
           y = dble(j) * dx_y + prob_lo_y
           x = (dble(i)+0.5d0) * dx_x + prob_lo_x
           vy(i,j,k) = ( (psi(i+1,j)+psi(i+1,j-1)) - (psi(i-1,j)+psi(i-1,j-1)) ) * (0.25d0/dx_x)
        end do
      end do
    end do
  
  #ifdef CUDA
    !$cuf kernel do(3) <<<*, (8,8,4), 0, cuda_streams(stream_from_index(idx),device_id)>>> 
    do k = vz_l3, vz_h3
      do j = vz_l2, vz_h2
        do i = vz_l1, vz_h1
          vz(i,j,k) = 1.d0
        end do
      end do
    end do
  #else
    vz = 1.d0
  #endif
  
  #ifdef CUDA
    ! will deallocate later
  #else
    call bl_deallocate(psi)
  #endif
  
  end subroutine get_face_velocity_doit

end module get_face_velocity_module
