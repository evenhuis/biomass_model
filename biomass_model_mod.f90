!--------------------------------------------------------------------------------
module model_mod
!--------------------------------------------------------------------------------
use ODE_mod
use carbonate_chemistry_mod
use spline_mod
use ext_driver_mod
use resize_mod

implicit none
integer,public  , parameter :: nstate =11, naux=10  ! size of state vector + aux variables  
double precision, parameter :: pi=3.141592654, eps=5d-9

integer ::  nmu0,  nmu1,    &
            np0   ,np1,     &
            nr0   ,nr1,     &
            nkl10 ,nkl11,   &
            nkl20 ,nkl21,   &
            nKm10 ,nKm11,   &
            nKm20 ,nKm21,   &
            nccm0, nccm1,   &
            nta0 , nta1,    &
            ntauR0, ntauR1, &
            ntauP0, ntauP1, &
            nfler0, nfler1, &
            nPQd0, nPQd1,   &
            nPQn0, nPQn1,   &
            nNCd0, nNCd1,   &
            nPRd0, nPRd1,   &
            nPRn0, nPRn1,   &
            nCRd0, nCRd1,   &
            nCRn0, nCRn1,   &
            nLRd0, nLRd1,   &
            nLRn0, nLRn1

double precision,allocatable,dimension(:) :: mut, muc,     &
                                             Pt,Pc,        & !
                                             Rt,Rc,        &
                                             kla1t, kla1c, &     ! rate of O2 production
                                             kla2t, kla2c, &
                                             kM1t,  kM1c,  &
                                             kM2t,  kM2c, &
                                             fccmt, fccmc, &
                                             tat,   tac,   &
                                             tauRt, tauRc, &
                                             tauPt, tauPc, &
                                             flert, flerc, &
                                             PQdt,  PQdc,  &  ! photosynthetic quotient day
                                             PQnt,  PQnc,  &  !                         night
                                             NCdt,  NCdc,  &  ! Nitrogen-Carbon ratio
                                             PRdt,  PRdc,  &  ! protein synthesis rate [protien uM / protein uM ] h
                                             PRnt,  PRnc,  & 
                                             CRdt,  CRdc,  &  ! Carb    synthesis rate [ carb   uM / protein uM ] h 
                                             CRnt,  CRnc,  &
                                             LRdt,  LRdc,  &  ! Lipid   synthesis rate [ lipid  uM / protein uM ] h
                                             LRnt,  LRnc      
   
                                             


! carbonate chemistry parameters
double precision :: S =33.,                     &  !Salinity
                    TK=273.15+22.,              &  ! Temp in kelvin
                    TA =2300., &                    ! alkinity in uM
                    KM=500.,  &
                    xO2 (2) = (/0.2094,0.0000/),&   ! partial pressure of O2
                    xCO2(2) = (/400d-6,1.00d0/), &  ! partial pressure of CO2
                    K1f, K2f 


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

contains

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function sim_PR_de( theta, time_steps ) result( y )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double precision, intent(in) :: theta(:)
double precision, intent(in) :: time_steps(:)
double precision             :: y(nstate+naux,size(time_steps))

double precision :: y0(nstate), y1(nstate)

integer          :: n, nstep
double precision :: dt, t0, t1, tt  ! time step info

double precision :: X_0       ! initial biomass
double precision :: O2_0      ! initial O2 concentration
double precision :: DIC_0     ! initial O2 concentration

double precision :: A(nstate,nstate)

double precision :: O2,DIC ,H, pH, CO3, HCO3, CO2,& ! PR=rate of photosynthesis or respiration
                    PM,P,R,PQ,KM1,KM2,fccm, f_photo, X, mu, kLa(2), pcl_0(3), N_0, P_0, TA_0, &
                    PR,CR,LR

double precision ::   nPn, nCb, nLb,  mass_0, pP0, pC0, pL0

integer :: i,nc, m,msub

nstep = size(time_steps)

y = 0.d0



! unpack the parameters
! initial parameters
!X_0     = theta( 1)
O2_0    = theta( 1)   ! oxygen
DIC_0   = theta( 2)
TA_0    = theta( 3)
N_0     = theta( 4)
P_0     = theta( 5)

mass_0  = theta( 6)*1d3 ! mg/L to ug/L
pP0     = theta( 7)
pC0     = theta( 8 )
pL0     = 1. - pP0 - pC0


nPn = mass_0*pP0/dot_product(prot_stoich,atomic_mass)
nCb = mass_0*pC0/dot_product(carb_stoich,atomic_mass)
nLb = mass_0*pC0/dot_product(lipd_stoich,atomic_mass)


pcl_0   = (/ nPn, nCb, nLb /)



muc     = theta(nmu0 :nmu1 )
Pc      = theta(nP0  :nP1  )*(24.)   ! rate of production in uM/day from uM/hr
Rc      = theta(nR0  :nR1  )*(24.)    
kla1c   = log(2.)/theta(nkl10:nkl11)*60*24.    ! convert to hour1 from halflife in min
kla2c   = log(2.)/theta(nkl20:nkl21)*60*24.
km1c    = theta(nkm10:nkm11)
km2c    = theta(nkm20:nkm21)
fccmc   = theta(nccm0:nccm1)
tac     = theta(nta0:nta1)
tauPc   = log(2.)/theta(ntauP0:ntauP1)*60.*24.
tauRc   = log(2.)/theta(ntauR0:ntauR1)*60.*24.
flerc   = theta(nfler0:nfler1)
PQdc    = theta(nPQd0:nPQd1)
PQnc    = theta(nPQn0:nPQn1)
NCdc    = theta(nNCd0:nNCd1)

PRdc    = theta(nPRd0:nPRd1)*24.
PRnc    = theta(nPRn0:nPRn1)*24.
CRdc    = theta(nCRd0:nCRd1)*24.
CRnc    = theta(nCRn0:nCRn1)*24.
LRdc    = theta(nLRd0:nLRd1)*24.
LRnc    = theta(nLRn0:nLRn1)*24.

! load the intial conditaions
y0 = 0.d0
y0(    1) = O2_0
y0(    2) = DIC_0
y0(    3) = TA
y0(    4) = N_0
y0(    5) = P_0
y0(    6) = 1.0
y0(    7) = 1.0
y0(8: 10) = pcl_0
y0(   11) = 1.d0


call CC_solve_DIC_Talk( DIC_0*1d-6, TA*1.d-6, TK, S, CO2, H,  pH=pH, HCO3=HCO3,CO3=CO3, &
                const=10 )
CO2 =   CO2*1d6
HCO3 = HCO3*1d6
CO3 =   CO3*1d6
y =0.d0

t0 = time_steps(1)
call calc_XPR(t0,y0, PR,CR,LR, pH,CO2,HCO3,CO3 )
y(        :nstate   ,1) = y0
y(nstate+1:nstate+ 7,1) =(/ pH,CO2,HCO3,CO3, PR,CR,LR /)
y(nstate+8          ,1) = dot_product( atomic_mass, prot_stoich )*y0( 8)
y(nstate+9          ,1) = dot_product( atomic_mass, carb_stoich )*y0( 9)
y(nstate+10         ,1) = dot_product( atomic_mass, lipd_stoich )*y0(10)

do n = 2, nstep
   t0 = time_steps(n-1)+eps
   t1 = time_steps(n  )-eps
   dt = t1-t0


   !kLa=0
   !kLa(1) = lint_1D( t0, kla1t,kla1c ,2 )
   !if( nkl21-nkl20>-1) kLa(2) = lint_1D( t0, kla2t,kla2c,1 )   
   !write(27,*) t0,kLa

   !write(24,*) t0,ext_gas(t0)
   !msub=7
   !do m = 1,msub-1
   !   tt = ((msub-m)*t0+m*t1)/dble(msub)
   !   write(24,*) tt,ext_gas(tt)
   !enddo
   !write(24,*) t1,ext_gas(t1)

   call RK4 ( y0, t0, dt, dy_PR,  y1 )
   !y1 =  magnus4_nl_real( y0, t0, dt, dA_PR )
 

   call calc_XPR(t0,y0, PR,CR,LR, pH,CO2,HCO3,CO3)
   y(         :nstate ,n) = y1
   y(nstate+1:nstate+7,n) =(/ pH,CO2,HCO3,CO3, PR,CR,LR /)
   y(nstate+8         ,n) = dot_product( atomic_mass, prot_stoich )*y1( 8)
   y(nstate+9         ,n) = dot_product( atomic_mass, carb_stoich )*y1( 9)
   y(nstate+10        ,n) = dot_product( atomic_mass, lipd_stoich )*y1(10)

   y0 = y1     ! Don't delete!
enddo


return
end function sim_PR_de

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
subroutine calc_XPR( t,y,  PR, CR, LR, pH,CO2, HCO3, CO3 )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double precision,intent(in) :: t, y(:) 
double precision,intent(out)::  PR,CR,LR , pH,CO2, HCO3, CO3

double precision :: t_dawn, t_dusk, X, O2,DIC,TA, I ,S,TK, KM,fler,t_resp,t_photo,&
                    fphoto,  fccm , H, KM1,KM2, P,R,PM, PQ, PCL(3)


t_dawn = 24*(t-floor(t))-8.
t_dusk = 24*(t-floor(t))-20.

if( t_dusk <  0. ) t_dusk=12.
if( t_dawn <  0. ) t_dawn=12.
if( t_dusk >12.0 ) t_dusk=12.
if( t_dawn >12.0 ) t_dawn=12.


O2  = y(1)
DIC = y(2)
TA  = y(3)
PCL = y(8:10)

I  = ext_light(t)              ! irradiance at top of column (uEin)
TK = ext_temp(t)+273.15d0
S  = ext_sal( t )

K1f = ext_K1f( t )
K2f = ext_K2f( t )

call CC_solve_DIC_Talk( DIC*1d-6, TA*1.d-6, TK, S, CO2, H,  pH=pH,HCO3=HCO3,CO3=CO3, &
               K1_f=K1f, K2_f=K2f, const=10 )
CO2  =  CO2*1d6
HCO3 = HCO3*1d6
CO3  =  CO3*1d6

fler   = 0.d0 ; if(nfler1-nfler0 > -1) fler   = spline_hc( t,flert, flerc )
t_resp = 1d6  ; if(ntauR1-ntauR0 > -1) t_resp = spline_hc( t,tauRt, tauRc )
t_photo= 1d6  ; if(ntauP1-ntauP0 > -1) t_photo= spline_hc( t,tauPt, tauPc )

!if( nmu1-nmu0 >-1 )then
!   if( mu_t0 <=t .and. t<=mu_t0+mu_dt ) mu = spline_hc( t, mut, muc )
!endif

if( t_dawn < 10.0 )then
   R = -spline_hc( t,Rt,Rc)!*(1.0+fler*(1.-exp(-t_dawn*t_resp)))
else
   R = -spline_hc( t,Rt,Rc)!*(1.0+fler*exp(-t_dusk*t_resp))
endif

P  = 0.d0
PQ = 1.d0

if ( I>0.1) then
   KM1=0.
   KM2=0.
   PM   = spline_hc( t, Pt,Pc)
   fccm = 0.
   if(nccm1-nccm0 >-1) fccm = spline_hc( t, fccmt,fccmc )
   if(nkm11-nkm10 >-1) KM1  = spline_hc( t, kM1t ,kM1c  )
   if(nkm21-nkm20 >-1) KM2  = spline_hc( t, kM2t ,kM2c )
   fphoto= (1.-exp(-t_dawn*t_photo))
   CR=PM*( (1-fccm)*HCO3/(HCO3+KM1) +  fccm*CO2/(CO2+KM2)  )*fphoto - R  ! carb synthesis rate

   PR=0.d0 ; if(nPRd1-nPRd0 >-1) PR  = spline_hc( t, PRdt ,PRdc  )
   LR=0.d0 ; if(nLRd1-nLRd0 >-1) LR  = spline_hc( t, LRdt ,LRdc  )
else
   CR=R
   PR=0.d0 ; if(nPRn1-nPRn0 >-1) PR  = spline_hc( t, PRnt ,PRnc  )
   LR=0.d0 ; if(nLRn1-nLRn0 >-1) LR  = spline_hc( t, LRnt ,LRnc )
endif

!print *,PRdt
!print *,PRdc
!stop 4
return
end subroutine calc_XPR

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function dy_PR( t, y ) result( dy )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double precision,intent(in) :: t,  y(:)
double precision            ::    dy(size(y))

double precision :: DIC,TA, X,O2, H, pH, HCO3, CO2, CO3,& ! PR=rate of photosynthesis or respiration
                     I,spg(2), O2_H(2), CO2_H(2), kLa(2),  PM,P,R,PQ, dTA, NC, PR,CR,LR, Pn,Cb,Ld, &
                      Ne,Pe, Ni,Pi

double precision :: t_dawn, t_dusk , mu=0.d0

t_dawn = t-floor((t)/24.)*24-8.
t_dusk = t-floor((t)/24.)*24-20.

O2 = y( 1)
DIC= y( 2)
TA = y( 3)

Ne = y( 4)
Pe = y( 5)

Ni = y( 6)
Pi = y( 7)

Pn = y( 8)
Cb = y( 9)
Ld = y(10)

I       = ext_light(t)              ! irradiance at top of column (uEin)
spg     = ext_gas  (t)              ! sparging rate


! calculate the Henry's law level of dissolved gases for the input streams
S  = ext_sal( t )
TK = ext_temp(t)+273.15d0
O2_H  = xO2  * K0_O2 ( TK, S ) * rho_sw( TK, S ) * 1000.
CO2_H = xCO2 * K0_CO2( TK, S ) * rho_sw( TK, S ) * 1000.

kLa = 0.d0
kLa(1) = spline_hc( t, kla1t,kla1c )
kLa(1) = lint_1D( t, kla1t,kla1c ,3 )
if( nkl21-nkl20>-1) kLa(2) = lint_1D( t, kla2t,kla2c,3 )

! calculate the growth, P and R rates
call calc_XPR( t, y,    PR,CR,LR,  pH, CO2,HCO3,CO3 )

dTA = 0.
if( nta1-nta0>-1) dTA = lint_1D( t, tat, tac,2 )
NC  = 0.
if( nNCd1-nNCd0>-1) NC = spline_hc( t, NCdt, NCdc )
if( I < 0.1 ) NC= 0

dy =  0.
dy(1) = Pn*(react_mat( 1,1)*PR + react_mat( 1,2)*CR + react_mat( 1,3)*LR)            +  dot_product( kLa       * spg, ( O2_H- O2) ) 
dy(2) = Pn*(react_mat( 2,1)*PR + react_mat( 2,2)*CR + react_mat( 2,3)*LR) +0.0*dTA   +  dot_product( kLa*0.893 * spg, (CO2_H-CO2) ) 
dy(3) = Pn* react_mat( 3,4)*PR                                            +    dTA

dy(4) = Pn*(react_mat( 6,1)*PR + react_mat( 6,2)*CR + react_mat( 6,3)*LR) 
dy(5) = Pn*(react_mat( 7,1)*PR + react_mat( 7,2)*CR + react_mat( 7,3)*LR)

dy(6) = 0.
dy(7) = 0.

dy( 8) = Pn*(react_mat( 8,1)*PR + react_mat( 8,2)*CR + react_mat( 8,3)*LR) 
dy( 9) = Pn*(react_mat( 9,1)*PR + react_mat( 9,2)*CR + react_mat( 9,3)*LR) 
dy(10) = Pn*(react_mat(10,1)*PR + react_mat(10,2)*CR + react_mat(10,3)*LR)
return
end function dy_PR

!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!function dA_PR( t, y ) result( dA )
!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!double precision,intent(in) :: t,  y(:)
!double precision            ::    dA(size(y),size(y))
!
!
!double precision :: DIC,TA, X,O2, PR ,H, pH, HCO3, CO2, CO3,& ! PR=rate of photosynthesis or respiration
!                     I,spg(2), O2_H(2), CO2_H(2), kLa(2),  PM,P,R,PQ, dTA, NC
!
!double precision :: t_dawn, t_dusk , mu=0.d0
!
!integer :: ny
!
!ny = size(y)
!
!
!O2 = y(1)
!DIC= y(2)
!TA = y(3)
!
!I       = ext_light(t)              ! irradiance at top of column (uEin)
!spg     = ext_gas  (t)              ! sparging rate
!
!
!! calculate the Henry's law level of dissolved gases for the input streams
!S  = ext_sal( t )
!TK = ext_temp(t)+273.15d0
!O2_H  = xO2  * K0_O2 ( TK, S ) * rho_sw( TK, S ) * 1000.
!CO2_H = xCO2 * K0_CO2( TK, S ) * rho_sw( TK, S ) * 1000.
!
!kla=0.d0
!!kLa(1) = spline_hc( t, kla1t,kla1c )
!kLa(1) = lint_1D( t, kla1t,kla1c ,2 )
!if( nkl21-nkl20>-1) kLa(2) = lint_1D( t, kla2t,kla2c,1 )
!
!
!! calculate the growth, P and R rates
!call calc_XPR( t, y,  PR,CR,LR,    pH, CO2,HCO3,CO3)
!
!dTA = 0.
!if( nta1-nta0>-1) dTA = spline_hc( t, tat, tac )
!NC  = 0.
!if( nNCd1-nNCd0>-1) NC = spline_hc( t, NCdt, NCdc )
!
!
!dA = 0
!dA(1,1 ) =           - dot_product( kLa       , spg              )
!dA(1,2 ) =  (P+R)/DIC
!dA(1,ny) =           + dot_product( kLa       * spg, ( O2_H    ) )
!
!
!dA(2,2 ) = -(P+R)/DIC - dot_product( kLa*0.893 , spg)*CO2/DIC
!dA(2,ny) =            + dot_product( kLa*0.893 * spg, (CO2_H    ) )
!
!if(I>0.1) dA(3,ny) =  NC*(P+R)
!
!
!return
!end function dA_PR

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
subroutine initialise()
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


call set_size( mut,   0 ) ; call set_size( muc,   0 )
call set_size( Pt,    0 ) ; call set_size( Pc,    0 )
call set_size( Rt,    0 ) ; call set_size( Rc,    0 )
call set_size( kla1t, 0 ) ; call set_size( kla1c, 0 )
call set_size( kla2t, 0 ) ; call set_size( kla2c, 0 )
call set_size( fccmt, 0 ) ; call set_size( fccmc, 0 )
call set_size( km1t,  0 ) ; call set_size( km1c,  0 )
call set_size( km2t,  0 ) ; call set_size( km2c,  0 )
call set_size( tat,   0 ) ; call set_size( tac,   0 )
call set_size( tauPt, 0 ) ; call set_size( tauPc, 0 )
call set_size( tauRt, 0 ) ; call set_size( tauRc, 0 )
call set_size( flert, 0 ) ; call set_size( flerc, 0 )
call set_size( PQdt , 0 ) ; call set_size( PQdc , 0 )
call set_size( PQnt , 0 ) ; call set_size( PQnc , 0 )
call set_size( NCdt , 0 ) ; call set_size( NCdc , 0 )

call set_size( PRdt , 0 ) ; call set_size( PRdc , 0 )
call set_size( PRnt , 0 ) ; call set_size( PRnc , 0 )
call set_size( CRdt , 0 ) ; call set_size( CRdc , 0 )
call set_size( CRnt , 0 ) ; call set_size( CRnc , 0 )
call set_size( LRdt , 0 ) ; call set_size( LRdc , 0 )
call set_size( LRnt , 0 ) ; call set_size( LRnc , 0 )

return
end subroutine initialise

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
subroutine set_spline_var( varname, Cp )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! sets the spline control points and saves the variables for the spline setup
character*(*)   , intent(in)   :: varname
double precision,intent(in) :: Cp(:)      ! control points of spline

integer :: nc, na


select case( trim(varname))
case("mu")   ; mut    = Cp
case("P")    ; Pt     = Cp 
case("R")    ; Rt     = Cp   
case("kla1") ; kla1t  = Cp    
case("kla2") ; kla2t  = Cp   
case("fccm" ); fccmt  = Cp   
case("km1" ) ; km1t   = Cp   
case("km2" ) ; km2t   = Cp   
case("dta" ) ; tat    = Cp
case("tauR") ; tauRt  = Cp
case("tauP") ; tauPt  = Cp   
case("fler") ; flert  = Cp   
case("PQd")  ; PQdt   = Cp   
case("PQn")  ; PQnt   = Cp
case("NCd")  ; NCdt   = Cp   

case("PRd")  ; PRdt   = Cp
!case("PRn")  ; PRnt   = Cp
!case("CRd")  ; CRdt   = Cp
!case("CRn")  ; CRnt   = Cp
!case("LRd")  ; LRdt   = Cp
!case("LRn")  ; LRnt   = Cp
end select
             ! no internal res inits 
nc = nstate -1  - 2 

print *,trim(varname)
call calc_offset(nc,mut  ,  nmu0,  nmu1) 
print *,"mu",nc,nmu0,nmu1
call calc_offset(nc,Pt   ,   nP0,   nP1) 
print *,"P",nc,nP0,nP1
call calc_offset(nc,Rt   ,   nR0,   nR1) 
call calc_offset(nc,kla1t, nkl10, nkl11)
call calc_offset(nc,kla2t, nkl20, nkl21)
call calc_offset(nc,fccmt, nccm0, nccm1)
call calc_offset(nc, km1t, nkm10, nkm11)
call calc_offset(nc, km2t, nkm20, nkm21)

call calc_offset(nc,  tat,  nta0,  nta1)

call calc_offset(nc,tauPt,ntauP0,ntauP1)
call calc_offset(nc,tauRt,ntauR0,ntauR1)
call calc_offset(nc,flert,nfler0,nfler1)

call calc_offset(nc, PQdt, nPQd0, nPQd1)
call calc_offset(nc, PQnt, nPQn0, nPQn1)
call calc_offset(nc, NCdt, nNCd0, nNCd1)

call calc_offset(nc, PRdt, nPRd0, nPRd1)
call calc_offset(nc, LRdt, nLRd0, nLRd1)

call calc_offset(nc, PRnt, nPRn0, nPRn1)
call calc_offset(nc, LRnt, nLRn0, nLRn1)



!na=size(mut  );if(na.eq.0) nmu0 = nc ; nmu1  = nc+na ; nc = nmu1
!na=size(Pt   );nP0   = nc+na ; nP1   = nc+na ; nc = nP1
!na=size(Rt   );nR0   = nc+na ; nR1   = nc+na ; nc = nR1
!na=size(kla1t);nkl10 = nc+na ; nkl11 = nc+na ; nc = nkl11
!na=size(kla2t);nkl20 = nc+na ; nkl21 = nc+na ; nc = nkl21
!na=size(fccmt);nccm0 = nc+na ; nccm1 = nc+na ; nc = nccm1
!na=size(km1t );nkm10 = nc+na ; nkm11 = nc+na ; nc = nkm11
!nkm20 = nc+1 ; nkm21 = nc+size(km2t)  ; nc = nkm21
!nta0  = nc+1 ; nta1  = nc+size(tat )  ; nc = nta1
!ntauP0= nc+1 ; ntauP1= nc+size(tauPt) ; nc = ntauP1
!ntauR0= nc+1 ; ntauR1= nc+size(tauRt) ; nc = ntauR1
!nfler0= nc+1 ; nfler1= nc+size(flert) ; nc = nfler1
!nPQd0 = nc+1 ; nPQd1 = nc+size(PQdt ) ; nc = nPQd1 
!nPQn0 = nc+1 ; nPQn1 = nc+size(PQnt ) ; nc = nPQn1 
!nNCd0 = nc+1 ; nNCd1 = nc+size(NCdt ) ; nc = nNCd1
!
!nPRd0 = nc+1 ; nPRd1 = nc+size(PRdt ) ; nc = nPRd1
!nPRn0 = nc+1 ; nPRn1 = nc+size(PRnt ) ; nc = nPRn1
!nCRd0 = nc+1 ; nCRd1 = nc+size(CRdt ) ; nc = nCRd1
!nCRn0 = nc+1 ; nCRn1 = nc+size(CRnt ) ; nc = nCRn1
!nLRd0 = nc+1 ; nLRd1 = nc+size(LRdt ) ; nc = nLRd1
!nLRn0 = nc+1 ; nLRn1 = nc+size(LRnt ) ; nc = nLRn1


return
end subroutine set_spline_var

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
subroutine calc_offset( noff,vec, n0,n1)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
integer, intent(inout)      :: noff
double precision,intent(in) :: vec(:)
integer, intent(out)        :: n0,n1

integer  :: nv

nv=size(vec)
if( nv == 0 )then
   n0 = noff+0 
   n1 = noff-1
else
   n0 = noff + 1
   n1 = noff + nv
   noff = noff+nv
endif
return
end subroutine calc_offset
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
subroutine set_size( arr, n )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double precision, allocatable,intent(inout) :: arr(:)
integer,intent(in) :: n
if( allocated(arr) )then
   if( size(arr) .ne. n )then
      deallocate(arr)
      allocate(arr(n))
   endif
else
   allocate(arr(n))
endif
end subroutine set_size

end module model_mod

