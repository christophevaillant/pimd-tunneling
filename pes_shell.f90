module pes_shell
  use pes,wp=>pes_wp
  implicit none
contains
  !==================================================
  ! init pes
  !==================================================
  subroutine pes_init()
    character (len=*), parameter :: dname='./coefs/'
    call pes0_init (dname)
    call pes1_init (pes_h4c3o2s_sysall)

    return
  end subroutine pes_init
  !==================================================
  ! fit in cartesian
  !==================================================
  function f(x) result(pot)
    real(kind=wp),dimension(:),intent(in)::x
    real(kind=wp)::pot
    ! ::::::::::::::::::::
    real(kind=wp),dimension(3,size(x)/3)::xn

    xn=reshape(x,(/3,size(x)/3/))

    pot=pes_h4c3o2s_pot(xn)
    
    !pot=pot-fmin
    !if ((pot-fmin)<-1.e-5) then
    !   print *, "hole!!!"
    !end if 
   
    return
  end function f
end module pes_shell
