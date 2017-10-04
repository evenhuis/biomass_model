!--------------------------------------------------------------------------------
module py_interface_mod
!--------------------------------------------------------------------------------
use model_mod
use ext_driver_mod


public  sim_de, py_setup_drivers, py_set_spline
contains

   

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function sim_de( np, theta, nt, time  )result(f)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
integer, parameter :: nstate=11, naux=10
integer,         intent(in) :: np,        nt
double precision,intent(in) :: theta(np) ,time(nt)
double precision            :: f(nstate+naux,nt)

integer                     :: n

f = sim_PR_de( theta, time )
return
end function sim_de

!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!function py_get_rate( np, theta, nx, x ) result( rate )
!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!integer,         intent(in) :: nx,        np
!double precision,intent(in) :: x(nx),     theta(np) 
!double precision            :: rate(nx)
!
!integer :: n
!
!rate = plot_rate( theta, x )
!return
!end function py_get_rate

!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!function py_get_Pmax( np, theta ) result( Pmax)
!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!integer,         intent(in) :: np
!double precision,intent(in) :: theta(np)
!double precision            :: Pmax
!
!Pmax = Pmax_hinge(theta)
!return
!end function py_get_Pmax

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
subroutine py_setup_drivers( fname )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
character*(*),intent(in) :: fname
call initialise()
call balance_reaction()
call setup_drivers(fname)
return
end subroutine py_setup_drivers

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
subroutine py_set_spline( varname, ncp, Cp )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! sets the spline coefficents up for the variables in the model
character*(*),intent(in)    :: varname
integer, intent(in)         :: ncp
double precision,intent(in) :: Cp(ncp)
call set_spline_var(varname, Cp )
return
end subroutine  py_set_spline



end module py_interface_mod


