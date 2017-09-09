module get_face_velocity_wrapper_module
#ifdef NCUDA
  use iso_c_binding
#endif
  use get_face_velocity_module, only: get_face_velocity_doit
  implicit none
  contains

  subroutine get_face_velocity(level, time, &
       vx, vx_l1, vx_l2, vx_l3, vx_h1, vx_h2, vx_h3, &
       vy, vy_l1, vy_l2, vy_l3, vy_h1, vy_h2, vy_h3, &
       vz, vz_l1, vz_l2, vz_l3, vz_h1, vz_h2, vz_h3, &
       dx, prob_lo &
  #ifdef NCUDA
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
  #ifdef NCUDA
      attributes(device) :: vx, vy, vz
      integer, value, intent(in) :: idx, device_id
      integer(kind=c_intptr_t), value, intent(in) :: tag
  #endif
    call get_face_velocity_doit(level, time, &
       vx, vx_l1, vx_l2, vx_l3, vx_h1, vx_h2, vx_h3, &
       vy, vy_l1, vy_l2, vy_l3, vy_h1, vy_h2, vy_h3, &
       vz, vz_l1, vz_l2, vz_l3, vz_h1, vz_h2, vz_h3, &
       dx, prob_lo &
  #ifdef NCUDA
       , idx, device_id, tag &
  #endif
       )
  
  end subroutine get_face_velocity

end module get_face_velocity_wrapper_module
