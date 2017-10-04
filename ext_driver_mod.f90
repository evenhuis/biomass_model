!--------------------------------------------------------------------------------
module ext_driver_mod
!--------------------------------------------------------------------------------
use read_data_mod
use spline_mod
! this module contains the functions and their setup for the external drivers

implicit none

! A series of models for algal growth
! the state vector is                  Ranges
! - - - - - - - - - - - - - - - - - - - - - - - - 
! external
!  1 O2  conc    (   umol/L)      ~  200 - 400
!  2 DIC conc    (   umol/L)      ~  100 - 3000
!  3 Alk conc    (   umol/L)      ~ 1000 - 3000
!  4 N   conc    (   umol/L)      ~    0 - 800
!  5 P   conc    (   umol/L)      ~    0 -  60
! - - - - - - - - - - - - - - - - - - - -  - - -  -
! internal
!  6 N   store   (   umol/L)      ~    0 - 0.05
!  7 P   store   (   umol/L)      ~    0 - 0.05
!  8 Protien     (   umol/L)        
!  9 Carb        (   umol/L)        
! 10 Lipid       (   umol/L)               
                                             !   C     H     O     N     P
double precision,dimension(5) ::  co2_stoich=(/  1.,   0.,   2.0 , 0.  , 0.   /),  &  ! DIC source
                                    n_stoich=(/  0.,   0.,   3.0 , 1.  , 0.   /),  &  ! N   source
                                    p_stoich=(/  0.,   0.,   4.0 , 0.  , 1.   /),  &  ! P   source 
                                 prot_stoich=(/  6.,  12.,   1.5 , 0.7 , 0.020/),  &  ! 50, 10, 10,10 ,3 by weight
                                 carb_stoich=(/  6.,  12.,   6.  , 0.  , 0.   /),  &  ! glucose type molecule
                                 lipd_stoich=(/ 36.,  80.,   6.  , 0.  , 0.   /),  &  ! trigly with chain length 12
                                 atomic_mass=(/ 12.,   1.,  16.  , 14. , 31.  /)


double precision :: react_mat(12,5)  

double precision, allocatable,dimension(:,:) :: &
   ext_light_arr, &
   ext_gas_arr,   &     ! CO2 file
   ext_air_arr,   &     ! Air file
   ext_temp_arr,  &  
   ext_sal_arr,   &
   ext_K1f_arr,   &
   ext_K2f_arr


logical ::  light_ext = .false., &
            air_ext   = .false., &
            temp_ext  = .false., &
            sal_ext   = .false., &
            K1f_ext   = .false., &
            K2f_ext   = .false.


contains

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
subroutine balance_reaction(show)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
logical, optional :: show

! this setups up the matrices to balance to the reactions
double precision :: nh2o
logical          :: show_eq=.false.
integer :: i

if( present(show) ) show_eq = show

react_mat = 0.

! Reaction 1: Carbohydrate synthesis
!              CO2 +  H2O    ->  carbs + O2
react_mat( 9,2) =   1.0                      ! moles of carb

react_mat( 2,2) =    -carb_stoich(1)         ! number of moles of DIC
nh2o            =-0.5*carb_stoich(2)         ! number of moles of H2O
react_mat( 1,2) =-0.5*(carb_stoich(3) + react_mat(2,2)*co2_stoich(3) + nh2o)    ! number of moles O2
if( show_eq )then
   write(6,'(a7  ,1(" + ",a7  )," -> ",a7  ,1(" + ",a7  ))') &
         "CO2",     "H2O",       "Carbs",   "O2"
   write(6,'(f7.3,1(" + ",f7.3)," -> ",f7.3,1(" + ",f7.3))') &
         -react_mat(2,2),-nh2o,  react_mat(9,2),react_mat(1,2) 
endif


! Reaction 2: Protien synthesis
!       Carbs  + Ni + Pi + H2O ->    Protien  + O2
react_mat(  8,1) = 1.0

react_mat( 6,1)  =-prot_stoich(4)                   ! N
react_mat( 7,1)  =-prot_stoich(5)                   ! P
react_mat( 9,1)  =-prot_stoich(1)/carb_stoich(1)    ! carbs given by carbon
nh2o             =-0.5*( prot_stoich(2) + carb_stoich(2)*react_mat(9,1) )
react_mat( 1,1)  =-0.5*( prot_stoich(3) + carb_stoich(3)*react_mat(9,1) + nh2o &
                                        +    n_stoich(3)*react_mat(6,1)        &
                                        +    p_stoich(3)*react_mat(7,1)  )

if( show_eq )then
   write(6,'(a7  ,3(" + ",a7  )," -> ",a7  ,1(" + ",a7  ))') &
              "Carbs","NO3","PO4","H2O","Protein","O2"
   write(6,'(f7.3,3(" + ",f7.3)," -> ",f7.3,1(" + ",f7.3))') &
         -react_mat( 9,1),-react_mat( 6,1),-react_mat( 7,1), -nh2o,react_mat(8,1), react_mat(1,1)
endif

! Reaction 3: Lipid synthesis
!     Carbs + H2O   -> Lipids + O2
react_mat(10,3)  =  1.d0
react_mat( 9,3)  = -lipd_stoich(1)/carb_stoich(1)  ! number of moles of car consumed
nh2o             = -0.5*( lipd_stoich(2) + carb_stoich(2)*react_mat(10,3) )
react_mat( 2,3)  = -0.5*( lipd_stoich(3) + carb_stoich(3)*react_mat(10,3) + nh2o)
if( show_eq )then
   write(6, '(a7  ,1(" + ",a7  ),"  -> ",a7  ,1(" + ",a7  ))') &
            "Carbs" ,"H2O","Lipid","O2"
   write(6,'(f7.3,1(" + ",f7.3),"  -> ",f7.3,1(" + ",f7.3))') &
         -react_mat( 9,3),-nh2o, react_mat(10,3),react_mat(2,3)
endif

! Reaction 3: N assim
!     Ne            ->  Ni
react_mat( 3,4) = +1    ! alkalinty increase
react_mat( 4,4) = -1    ! NO3- decrease
react_mat( 6,4) = +1    ! NO2 incrase

! Reaction 3: N assim
!     Pe            ->  Pi
react_mat( 5,5) = -1
react_mat( 7,5) = +1

if( show_eq )then
   2 format(12(1x,f6.3))
   do i = 1,5
      write(6,2) react_mat(:,i)
   enddo
endif


return
end subroutine balance_reaction

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! External drivers
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function ext_light( t ) result ( I )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double precision,intent(in) :: t
double precision            :: I

double precision :: tmod, ti

tmod = (t - floor( t ))*24

if( light_ext )then
   I =  lint_2D( t, ext_light_arr, 2 )
else
   I=0.
   if( 9.<tmod .and. tmod<21.) I=1.
endif

return
end function ext_light

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function ext_gas( t ) result( spg )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double precision,intent(in) :: t
double precision            :: spg(2)

double precision :: tmod, ti

ti = t*24.d0*60d0                  ! time in minutes
ti = ti-30.
tmod = ti-floor(ti/120.)*120.

if( .not. air_ext )then
   spg(1) = 1.d0
   if( tmod < 10.d0 )  spg(1) = 0.d0
else
   spg(1) = lint_2D( t, ext_air_arr,2)
endif
spg(2) = lint_2D( t, ext_gas_arr,2)

return
end function ext_gas

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function ext_temp( t ) result( temp )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double precision,intent(in) :: t
double precision            :: temp

temp = lint_2D( t, ext_temp_arr,2) 

return
end function ext_temp

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function ext_sal ( t ) result( sal )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double precision,intent(in) :: t
double precision            :: sal

if( sal_ext )then
   sal  = lint_2D( t, ext_sal_arr,2) 
else
   sal = 33.
endif

return
end function ext_sal

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function ext_K1f( t ) result( K1f )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double precision,intent(in) :: t
double precision            :: K1f

if( K1f_ext )then
   K1f  = lint_2D( t, ext_K1f_arr,0)
else
   K1f = 1.
endif

return
end function ext_K1f

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function ext_K2f( t ) result( K2f )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double precision,intent(in) :: t
double precision            :: K2f

if( K2f_ext )then
   K2f  = lint_2D( t, ext_K2f_arr,0)
else
   K2f = 1.
endif

return
end function ext_K2f



! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! util functions

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
subroutine setup_drivers( filename )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
character*(*), intent(in) :: filename

character*(200) :: line, key,val,  var,dataf
integer         :: ius           ! underscore location
logical         :: lex
integer :: i

call init_arrays()

open ( unit=1, file=filename )
do
   read(1,'(a200)',end=5)    line
   if( len(trim(line))== 0 ) cycle

   read(line,*,err=10,end=5) key,val   

   ! split the key on _
   ius = scan( key,"_")
   if( ius > 0 )then
      var      = key(1    :ius-1) 
      dataf    = key(ius+1:    )
      
      select case(trim(dataf))
      case("file")
         inquire(file=trim(val), exist=lex)
         if( .not. lex )then
            print *,'setup_drivers : file not found'
            print *,    trim(val)
            stop 1
         endif

         open(unit=2,file=val )
         select case(trim(var))
         case("light") 
            call read_data_nf(ext_light_arr,2)
            light_ext = .true.
         case("gas"  ) ; call read_data_nf(ext_gas_arr  ,2)
         case("air"  ) 
            call read_data_nf(ext_air_arr  ,2)   
            air_ext = .true.
         case("temp" ) ; call read_data_nf(ext_temp_arr ,2)   
         end select
         close(unit=2)

      case("val")
         select case(trim(var))
         case("light") ; read(val,*) ext_light_arr(2,1) ; light_ext = .true.
         case("gas"  ) ; read(val,*) ext_gas_arr  (2,1) 
         case("temp" ) ; read(val,*) ext_temp_arr (2,1) ; temp_ext = .true.
         case("sal"  ) ; read(val,*) ext_sal_arr (2,1)  ; sal_ext = .true.
         end select 
      end select
   endif

   10 continue
enddo
5 continue
close( unit=1 )

return
end subroutine setup_drivers

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
subroutine init_arrays()
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if( allocated( ext_light_arr ) )deallocate(ext_light_arr)
if( allocated( ext_gas_arr   ) )deallocate(ext_gas_arr  )
if( allocated( ext_air_arr   ) )deallocate(ext_air_arr  )
if( allocated( ext_temp_arr  ) )deallocate(ext_temp_arr )
if( allocated( ext_sal_arr  )  )deallocate(ext_sal_arr )

allocate( ext_light_arr(2,1) ) ; ext_light_arr = 0
allocate( ext_gas_arr  (2,1) ) ; ext_gas_arr   = 0
allocate( ext_air_arr  (2,1) ) ; ext_air_arr   = 0
allocate( ext_temp_arr (2,1) ) ; ext_temp_arr  = 0
allocate( ext_sal_arr  (2,1) ) ; ext_sal_arr   = 0
return
end subroutine init_arrays

end module ext_driver_mod
