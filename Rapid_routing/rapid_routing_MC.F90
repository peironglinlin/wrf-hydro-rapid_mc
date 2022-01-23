!----------------------------------------------------------------
!------Peirong Lin, April 2017 for Muskingum-Cunge RAPID--------
!----------------------------------------------------------------
module rapid_routing_MC

IMPLICIT NONE

CONTAINS


!--- update MC routing params every time step----
subroutine routing_params_MC(comid,bw,ss,bs,ll,nchan,hinit,tp,acat,&
                llcat,ncat,slp,tlag,&
                ZV_Qext,ZV_QoutinitR, &
                ZV_C1,ZV_C2,ZV_C3,MC_C4,ZM_A,ZV_k,ZV_x,&
                ZV_QoutR,ZV_QoutbarR)
use rapid_var, only : JS_riv_tot,IS_riv_tot, &
                ZV_QoutprevR,IV_riv_loc1,IV_riv_index,&
                JS_R,IS_R,ZS_dtR,ZS_TauR,ierr,ZV_one,ZS_one,&
                ZM_Net,ZV_b,JS_riv_bas,IM_index_up,IS_riv_bas,&
                BS_opt_influence,IV_nbup,&
                ZV_babsmax,ZV_bhat,ZV_QoutRabsmin,ZV_QoutRabsmax,&
                ZV_b,ZV_babsmax,ZV_bhat,ksp,rank,&
                ZV_Qlat_r,ZV_Qstore,ZV_Qstoreprev,&
                IS_opt_ov,beta, & !vars for overland flow
                speed1,speed2lb,speed2b,speed3b,speed3db !for speed up certain rivers


implicit none

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscksp.h"
#include "finclude/petscpc.h"
#include "finclude/petscviewer.h"
#include "finclude/petsclog.h"

!---variables in/out--
integer, dimension(IS_riv_tot),intent(in) :: comid,tp
real, dimension(IS_riv_tot),intent(in) :: bw,ss,bs,ll,nchan,hinit,acat !channel
real, dimension(IS_riv_tot),intent(in) :: llcat,ncat,tlag        !overland
real, dimension(IS_riv_tot),intent(inout) :: slp !inout for modifying zero slp
Vec,intent(in) :: ZV_Qext,ZV_QoutinitR
Vec,intent(inout) :: ZV_C1,ZV_C2,ZV_C3,MC_C4,ZM_A,ZV_k,ZV_x
!---variables for calc. in this subroutine------------
Vec,intent(out) :: ZV_QoutR,ZV_QoutbarR
PetscInt :: IS_localsize,JS_localsize
PetscScalar, pointer :: ZV_QoutR_p(:),ZV_QoutprevR_p(:),ZV_QoutinitR_p(:),     &
                        ZV_QoutbarR_p(:),ZV_Qext_p(:),ZV_C1_p(:),ZV_C2_p(:),   &
                        ZV_C3_p(:),ZV_b_p(:),MC_C4_p(:),                       &
                        ZV_babsmax_p(:),ZV_QoutRabsmin_p(:),ZV_QoutRabsmax_p(:)
!lpr for overland
PetscScalar, pointer :: ZV_Qlat_r_p(:),ZV_Qstore_p(:),ZV_Qstoreprev_p(:) 
real :: shapefn,Qj,error        !intermediate variables
integer :: maxiter              !intermediate variables
real, dimension(IS_riv_tot) :: z,h   !intermediate: can output later
real, dimension(IS_riv_tot) :: Tw,Ck   !intermediate: can output later
real, dimension(IS_riv_tot) :: tch,tov  !channel and overland flow travel time
PetscScalar,dimension(:), allocatable :: Km,x,C1,C2,C3,C4
real :: D, D1
real :: area_t,wp_t,hr_t  !hr: hydraulic radius
real, dimension(IS_riv_tot) :: ZV_qp
PetscScalar, pointer :: ZV_qp_p(:)
real :: qnow !temporary variable
real :: tmpperc !temporary variable
!real :: beta !calibration factor for tlag
integer :: idebug,inegative  !integer for river reach ID for debugging
integer :: tmpid
real :: tmpfrac !tmp variable for fraction of lateral flow
real :: m !C=m*V， Q=a*A^m


idebug=1013 !change this value for debugging


allocate(Km(IS_riv_tot))
allocate(x(IS_riv_tot))
allocate(C1(IS_riv_tot))
allocate(C2(IS_riv_tot))
allocate(C3(IS_riv_tot))
allocate(C4(IS_riv_tot))


print *,'******** Updating MC routing_params ********'
call VecGetLocalSize(ZV_Qext,IS_localsize,ierr)
call VecSet(ZV_QoutbarR,0*ZS_one,ierr)     !Qoutbar=0   ;segment fault when use 0*ZS_one
call VecCopy(ZV_QoutinitR,ZV_QoutprevR,ierr)               !QoutprevR=QoutinitR
call VecGetArrayF90(ZV_Qext,ZV_Qext_p,ierr)
call VecGetArrayF90(ZV_QoutprevR,ZV_QoutprevR_p,ierr)
!---for overland
call VecGetArrayF90(ZV_Qlat_r,ZV_Qlat_r_p,ierr)
call VecGetArrayF90(ZV_Qstore,ZV_Qstore_p,ierr)
call VecGetArrayF90(ZV_Qstoreprev,ZV_Qstoreprev_p,ierr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do JS_localsize = 1,IS_localsize

!!!----calculate Tov-----------------------
if(slp(JS_localsize).le.0. .or. slp(JS_localsize).ge.1.) then
        slp(JS_localsize)=0.01
end if
if(isnan(ZV_Qext_p(JS_localsize)).or.ZV_Qext_p(JS_localsize).le.0.) then
tov(JS_localsize)=0.
else
tov(JS_localsize)=(llcat(JS_localsize)*ncat(JS_localsize))**0.6&
        /(((slp(JS_localsize))**0.3)*((ZV_Qext_p(JS_localsize))**0.4))
end if
!-----------overland flow considerations--------------------------
!calculate Qlat,r after consider overland time lag for runoff
if(IS_opt_ov.eq.0) then
ZV_Qlat_r_p(JS_localsize)=ZV_Qext_p(JS_localsize)  !m3/s
tmpperc=1.
ZV_Qstore_p(JS_localsize)=0.
ZV_Qstoreprev_p(JS_localsize)=0.
else !If user choose to turn on OVERLAND time delay
!ZV_Qlat_r_p(JS_localsize)=(ZV_Qext_p(JS_localsize)*ZS_TauR+ZV_Qstoreprev_p(JS_localsize))*&
!        (1.-exp(-1.*beta/(tov(JS_localsize)/3600.)/ncat(JS_localsize)))/ZS_TauR
!consider time lag for water stored in the overland plane
ZV_Qlat_r_p(JS_localsize)=(ZV_Qext_p(JS_localsize)*ZS_TauR+ZV_Qstoreprev_p(JS_localsize))*&
        (1.-exp(-1.*(tlag(JS_localsize)*beta/tov(JS_localsize))))/ZS_TauR
tmpperc=ZV_Qlat_r_p(JS_localsize)/ZV_Qext_p(JS_localsize)
ZV_Qstore_p(JS_localsize)=ZV_Qext_p(JS_localsize)*ZS_TauR+ZV_Qstoreprev_p(JS_localsize)-&
       ZV_Qlat_r_p(JS_localsize)*ZS_TauR  
ZV_Qstoreprev_p(JS_localsize)=ZV_Qstore_p(JS_localsize) !+ZV_Qstoreprev_p(JS_localsize)!update Qstoreprev
end if
!-------------------------------------------------------------------



!!!ZV_qp(JS_localsize) = ZV_QoutprevR_p(JS_localsize)+ZV_Qlat_r_p(JS_localsize)!ZV_Qext_p(JS_localsize)
call VecGetArrayF90(ZV_b,ZV_b_p,ierr)
!please modify back as this might be the correct one to use
!ZV_qp(JS_localsize) = ZV_b_p(JS_localsize)+ZV_Qlat_r_p(JS_localsize)  !NWM1
!ZV_qp(JS_localsize) = ZV_QoutprevR_p(JS_localsize)+ZV_Qext_p(JS_localsize)!NWM
ZV_qp(JS_localsize) = ZV_QoutprevR_p(JS_localsize)+ZV_Qlat_r_p(JS_localsize)


!---add a few lines to control the calculation of m--
!       -- threshold: 0.1 as arbitary--
if(ZV_Qlat_r_p(JS_localsize)/ZV_QoutprevR_p(JS_localsize) .le. 0.1) then
m=1.0   !flow travel using its velocity not flood wave celerity
else
m=5./3.
end if
print *,'       *********** m = ',m,' ****************'


!**Note: all c1,c2,c3,c4 calculations done only when ZV_qp>0, otherwise, they equal 0*****
!if(ZV_qp(JS_localsize).gt.0.) then
!!-----------------
z(JS_localsize) = 1./ss(JS_localsize)
h(JS_localsize) = sqrt(ZV_qp(JS_localsize))*0.1 !initial guess of h
!find a depth (h) which makes the Qi=quc (error<0.01)
!LPR note: Newton-Raphson iteration method
error = 1.0  
maxiter = 0
do while (error.gt.0.0001.and.maxiter<100.)
        maxiter=maxiter+1
        !---trapezoidal channel shape function
        shapefn=SHAPE1(bw(JS_localsize),z(JS_localsize),h(JS_localsize))
        Qj=FLOW(nchan(JS_localsize),bs(JS_localsize),&
                bw(JS_localsize),h(JS_localsize),z(JS_localsize))
        !h(JS_localsize)=h(JS_localsize)-(1-ZV_qp(JS_localsize)/Qj)/shapefn
        !error=abs((Qj-ZV_qp(JS_localsize))/ZV_qp(JS_localsize)) 
        h(JS_localsize)=h(JS_localsize)-(1-ZV_QoutprevR_p(JS_localsize)/Qj)/shapefn
        error=abs((Qj-ZV_QoutprevR_p(JS_localsize))/ZV_QoutprevR_p(JS_localsize))
        !Now use ZV_qp to find h; but this seems not reasonable -- not used 
end do
if(JS_localsize.eq.idebug) then
print *,'------------- DEBUG INFORMATION FOR OVERLAND FLOW-----------------'
print *,'       ****** Processing COMID = ',comid(idebug),' ***************' 
print *,'       beta (adjustment factor)=',beta
print *,'       n_catchment=',ncat(JS_localsize)
print *,'       Tov=',tov(JS_localsize),' ;tlag=',tlag(JS_localsize),'-----'
print *,'       ZV_qp = ',ZV_qp(idebug)
print *,'       ZV_b_p = ',ZV_b_p(idebug)
print *,'       ZV_QoutprevR = ',ZV_QoutprevR_p(idebug)
print *,'       ZV_Qext = ',ZV_Qext_p(idebug)
print *,'	ZV_Qstoreprev=',ZV_Qstoreprev_p(JS_localsize)
print *,'       ZV_Qlat_r = ',ZV_Qlat_r_p(idebug)
print *,'   --- Percentage Qlat reaching channels ',tmpperc
print *,'   --- MC params: bw = ',bw(idebug),' m; ss=',ss(idebug)
print *,'   --- LPR CHECK: h (after iteration) =',h(idebug), 'm '
print *,'-------------------------------------------------------------------'
end if
call VecRestoreArrayF90(ZV_b,ZV_b_p,ierr)

Tw(JS_localsize)=bw(JS_localsize)+2*z(JS_localsize)*h(JS_localsize)
area_t=(bw(JS_localsize)+z(JS_localsize)*h(JS_localsize))*h(JS_localsize)
wp_t=bw(JS_localsize)+2*h(JS_localsize)*sqrt(1+(z(JS_localsize)**2)) !different from WRF-Hydro (h to z)
hr_t=area_t/wp_t   !hr contains the information from Q
!wave celerity on Manning's equation (beta=5/3),following WRF-Hydro, pg 287 Chow, Mdt, Mays
!Ck(JS_localsize)=sqrt(bs(JS_localsize))/nchan(JS_localsize)*(5./3.)*(hr_t**0.667)
!!!!wave celerity based on McDonnell and Beven (2014, WRR, Table 1(1))

!!!!************LPR SPEEDING UP CERTAIN RIVER REACHES *******************
tmpid=comid(JS_localsize)
!-------------------------------------------------------------------------
!******** Speed1: 38 reaches for Order 1 Streams (*********************
if(tmpid.eq.1628079 .or. tmpid.eq.1628081 .or. tmpid.eq.1628083 &
.or. tmpid.eq.1628085 .or. tmpid.eq.1628091 .or. tmpid.eq.1628093 .or.tmpid.eq.1628095 &
.or. tmpid.eq.1628101 .or. tmpid.eq.1628103 .or. tmpid.eq.1628105 .or. tmpid.eq.1628113 &
.or. tmpid.eq.1628115 .or. tmpid.eq.1628143 .or. tmpid.eq.1628151 .or. tmpid.eq.1628155 &
.or. tmpid.eq.1628159 .or. tmpid.eq.1628163 .or. tmpid.eq.1628165 .or. tmpid.eq.1628171 &
.or. tmpid.eq.1628173 .or. tmpid.eq.1628179 .or. tmpid.eq.1628181 .or. tmpid.eq.1628183 &
.or. tmpid.eq.1628185 .or. tmpid.eq.1628187 .or. tmpid.eq.1628193 .or. tmpid.eq.1628201 &
.or. tmpid.eq.1628203 .or. tmpid.eq.1628205 .or. tmpid.eq.1628213 .or. tmpid.eq.1628225 &
.or. tmpid.eq.1628231 .or. tmpid.eq.1628235 .or. tmpid.eq.1628257 .or. tmpid.eq.1628519 &
.or. tmpid.eq.1628523 .or. tmpid.eq.1630201 .or. tmpid.eq.1630233) then
Ck(JS_localsize)=speed1*sqrt(bs(JS_localsize))/nchan(JS_localsize)*&
         m*hr_t**0.667
!-------------------------------------------------------------------------
!******** Speed2: 11 reaches for Order 2 Little Blanco *************************
else if(tmpid.eq.1628229 .or. tmpid.eq.1628233 .or. tmpid.eq.1628237  &
.or. tmpid.eq.1628255 .or. tmpid.eq.1628259 .or. tmpid.eq.1628261 .or.tmpid.eq.1630203 &
.or. tmpid.eq.1630205 .or. tmpid.eq.1630213 .or. tmpid.eq.1633013 .or. tmpid.eq.1633027) then
Ck(JS_localsize)=speed2lb*sqrt(bs(JS_localsize))/nchan(JS_localsize)*&
         m*hr_t**0.667
!-------------------------------------------------------------------------
!******** Speed3: 14 reaches for Order 2 Blanco ********************
else if(tmpid.eq.1628097 .or. tmpid.eq.1628099 .or. tmpid.eq.1628111 .or. tmpid.eq.1628117 &
.or. tmpid.eq.1628121 .or. tmpid.eq.1628123 .or. tmpid.eq.1628131 .or. tmpid.eq.1628135 &
.or. tmpid.eq.1628133 .or. tmpid.eq.1628137 .or. tmpid.eq.1628139 .or. tmpid.eq.1628141 &
.or. tmpid.eq.1628509 .or. tmpid.eq.1628529  !lpr more needed
.or. tmpid.eq.1628169 .or. tmpid.eq.1330582 .or. tmpid.eq.1628125) then
Ck(JS_localsize)=speed2b*sqrt(bs(JS_localsize))/nchan(JS_localsize)*&
         m*hr_t**0.667
!-------------------------------------------------------------------------
!******** Speed4: 14 reaches for Order 3 Blanco ***********************
else if(tmpid.eq.1628127 .or. tmpid.eq.1628145 .or. tmpid.eq.1628149 .or. tmpid.eq.1628157 &
.or. tmpid.eq.1628161 .or. tmpid.eq.1628167 .or.tmpid.eq.1628175 &
.or. tmpid.eq.1628177 .or. tmpid.eq.1628199 .or. tmpid.eq.1628207 &
.or. tmpid.eq.1628515 .or. tmpid.eq.1628527 .or. tmpid.eq.1628631 .or. tmpid.eq.1628533) then
Ck(JS_localsize)=speed3b*sqrt(bs(JS_localsize))/nchan(JS_localsize)*&
         m*hr_t**0.667
!-------------------------------------------------------------------------
!******** Speed4: 11 reaches for Order 3 Downstream Blanco ****************
else if(tmpid.eq.1628227 .or. tmpid.eq.1628245 .or. tmpid.eq.1628253 &
.or. tmpid.eq.1630217 .or. tmpid.eq.1630223 .or. tmpid.eq.1630235 .or. tmpid.eq.1630239 &
.or. tmpid.eq.1630245 .or. tmpid.eq.1630249 .or. tmpid.eq.1630251 .or. tmpid.eq.1630253) then 
Ck(JS_localsize)=speed3db*sqrt(bs(JS_localsize))/nchan(JS_localsize)*&
         m*hr_t**0.667
!***********************************************************************
!Ck(JS_localsize)=speeding*sqrt(bs(JS_localsize))/nchan(JS_localsize)*&
!         (5./3.)*hr_t**0.667
!***********************************************************************
else
Ck(JS_localsize)=sqrt(bs(JS_localsize))/nchan(JS_localsize)*m*(hr_t**0.667)
end if


D1= ZV_qp(JS_localsize)/(2*Tw(JS_localsize)*bs(JS_localsize)) !pg295
x(JS_localsize)=0.5-(D1/(Ck(JS_localsize)*ll(JS_localsize)))  !pg295
!x(JS_localsize)=0.01 !!!LPRLPRLPR!!!!! REMEMBER TO DELETE THIS LATER!!!!


!lpr: deal with erraneous places-----
if(x(JS_localsize).le.0..or.isnan(x(JS_localsize))) then
    x(JS_localsize)=0.3 !if <0, set to 0.3
end if
if(Ck(JS_localsize).eq.0. .or. isnan(Ck(JS_localsize))) then
    Ck(JS_localsize) = 0.35*5./3.
end if


tch(JS_localsize)=ll(JS_localsize)/Ck(JS_localsize)
Km(JS_localsize)=tch(JS_localsize) !+tov(JS_localsize) !ll(JS_localsize)/Ck(JS_localsize)
D = Km(JS_localsize)*(1-x(JS_localsize))+ZS_dtR/2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---lpr add some controls to deal with Connectors and short reaches
!if(tp(JS_localsize).eq.0 .or. ll(JS_localsize).le.50.) then
if(ll(JS_localsize).le.50.) then  !for short reaches: Qin=Qout
C1(JS_localsize)=0.
C2(JS_localsize)=1.     !Qin=Qout
C3(JS_localsize)=0.
else
C1(JS_localsize)=(-Km(JS_localsize)*x(JS_localsize)+ZS_dtR/2)/D
C2(JS_localsize)=(ZS_dtR/2+Km(JS_localsize)*x(JS_localsize))/D
C3(JS_localsize)=(Km(JS_localsize)*(1-x(JS_localsize))-ZS_dtR/2)/D
end if
C4(JS_localsize)=ZV_Qlat_r_p(JS_localsize)*ZS_dtR/D

!else !*** Note: if ZV_qp<=0, Qout=0
!C1(JS_localsize)=0.
!C2(JS_localsize)=0.
!C3(JS_localsize)=0.
!C4(JS_localsize)=0.
!end if

end do


!---print out information after C1, C2, C3, C4 are calculated for each channel---
print *,'   --- LPR CHECK: Ck = ',Ck(idebug),' m/s'
print *,'       --- LPR CHECK: Tw = ',Tw(idebug),' m'
print *,'   --- LPR CHECK: Tch = ',tch(idebug),'s; Tov = ',tov(idebug),'s'
print *,'   --- LPR CHECK: Km = ',Km(idebug),'; *** X = ',x(idebug)
print *,'   --- C1 = ',C1(idebug),'; C2=',C2(idebug),' ; C3=',C3(idebug),' ; C4=',C4(idebug)
print *,'   --- average celerity = ', real(sum(Ck)/size(Ck)), ' -----'




if (rank==0) then
call VecSetValues(ZV_C1,IS_riv_bas,IV_riv_loc1,&
               C1(IV_riv_index),INSERT_VALUES,ierr)
call VecAssemblyBegin(ZV_C1,ierr)
call VecAssemblyEnd(ZV_C1,ierr)
end if 

if (rank==0) then
call VecSetValues(ZV_C2,IS_riv_bas,IV_riv_loc1,&
               C2(IV_riv_index),INSERT_VALUES,ierr)
call VecAssemblyBegin(ZV_C2,ierr)
call VecAssemblyEnd(ZV_C2,ierr)
end if 

if (rank==0) then
call VecSetValues(ZV_C3,IS_riv_bas,IV_riv_loc1,&
               C3(IV_riv_index),INSERT_VALUES,ierr)
call VecAssemblyBegin(ZV_C3,ierr)
call VecAssemblyEnd(ZV_C3,ierr)
end if 

if (rank==0) then
call VecSetValues(MC_C4,IS_riv_bas,IV_riv_loc1,&
               C4(IV_riv_index),INSERT_VALUES,ierr)
call VecAssemblyBegin(MC_C4,ierr)
call VecAssemblyEnd(MC_C4,ierr)
end if 

if (rank==0) then
call VecSetValues(ZV_k,IS_riv_bas,IV_riv_loc1,&
               Km(IV_riv_index),INSERT_VALUES,ierr)
call VecAssemblyBegin(ZV_k,ierr)
call VecAssemblyEnd(ZV_k,ierr)
end if 

if (rank==0) then
call VecSetValues(ZV_x,IS_riv_bas,IV_riv_loc1,&
               x(IV_riv_index),INSERT_VALUES,ierr)
call VecAssemblyBegin(ZV_x,ierr)
call VecAssemblyEnd(ZV_x,ierr)
end if 

!---after calculate C1, C2, C3, update A linear system matrix---
call MatCopy(ZM_Net,ZM_A,DIFFERENT_NONZERO_PATTERN,ierr)   !A=Net
call MatDiagonalScale(ZM_A,ZV_C1,ZV_one,ierr)              !A=diag(C1)*A
call MatScale(ZM_A,-ZS_one,ierr)                           !A=-A
call MatShift(ZM_A,ZS_one,ierr)                  
print *,'****** Muskingum-Cunge parameter UPDATED! ********'



call VecGetArrayF90(ZV_C1,ZV_C1_p,ierr)
call VecGetArrayF90(ZV_C2,ZV_C2_p,ierr)
call VecGetArrayF90(ZV_C3,ZV_C3_p,ierr)
call VecGetArrayF90(MC_C4,MC_C4_p,ierr)



print *,'************* ROUTING STARTED!!! ******************'
!time loop
do JS_R=1,IS_R
call VecAXPY(ZV_QoutbarR,ZS_one/IS_R,ZV_QoutprevR,ierr)    !ZV_QoutbarR = ZV_QoutbarR+ZV_QoutprevR/dt
call VecGetArrayF90(ZV_QoutbarR,ZV_QoutbarR_p,ierr)
call MatMult(ZM_Net,ZV_QoutprevR,ZV_b,ierr)                !b2=Net*Qoutprev
call VecGetArrayF90(ZV_b,ZV_b_p,ierr)  !lpr note: ZV_b = quc, the inflow coming from upstream
call VecGetArrayF90(ZV_QoutprevR,ZV_QoutprevR_p,ierr)


do JS_localsize=1,IS_localsize
     ZV_b_p(JS_localsize)=ZV_b_p(JS_localsize)*ZV_C2_p(JS_localsize)           &
                         +(ZV_C1_p(JS_localsize)+ZV_C2_p(JS_localsize))        &
                         !*ZV_Qext_p(JS_localsize)                              &
                         *ZV_Qlat_r_p(JS_localsize)  &!if not turn on overland, ZV_Qlat_r=ZV_Qext_p
                         +ZV_C3_p(JS_localsize)*ZV_QoutprevR_p(JS_localsize)                                            
end do
print *,'   --- ZV_b_p(JS_localsize) = ',ZV_b_p(idebug)


call VecGetArrayF90(ZV_QoutR,ZV_QoutR_p,ierr)
inegative=0
do JS_riv_bas=1,IS_riv_bas
     ZV_QoutR_p(JS_riv_bas)=ZV_b_p(JS_riv_bas)                                 &
                            +sum(ZV_C1_p(JS_riv_bas)                           &
                                  *ZV_QoutR_p(IM_index_up(JS_riv_bas,1:        &
                                   IV_nbup(IV_riv_index(JS_riv_bas))))) 
     if (ZV_QoutR_p(JS_riv_bas).lt.0) then
     print *,'      !!! ************** Q<0 calculated after routing => set to zero *******!!!'
     !print *,'     --- COMID = ',comid(JS_localsize)
     !print *,'     --- ZV_QoutR_p(JS_riv_bas)=',ZV_QoutR_p(JS_riv_bas)
     !print *,'     --- C1 = ',ZV_C1_p(idebug),'; C2=',ZV_C2_p(idebug),&
     !   ' ;C3=',ZV_C3_p(idebug),';C4=',ZV_C1_p(JS_localsize)+ZV_C2_p(JS_localsize)
     !print *,'     --- Qup,p=',sum(ZV_QoutR_p(IM_index_up(JS_riv_bas,1:        &
     !                              IV_nbup(IV_riv_index(JS_riv_bas)))))
     !print *,'     --- Qdn,p=',ZV_QoutprevR_p(JS_localsize)
     inegative=inegative+1
     ZV_QoutR_p(JS_riv_bas)=0.
     end if
end do
!print *,'ZV_QoutR (Qout) = ',ZV_QoutR_p(1026)
!Taking into account the knowledge of how many upstream locations exist.
!Similar to exact preallocation of network matrix
print *,'**** AFTER ROUTING *******************************'
print *,'   !!!****** inegative = ',inegative,' *******!!!'
print *,'   --- ZV_QoutR_p(Qout) = ',ZV_QoutR_p(idebug)
call VecRestoreArrayF90(ZV_QoutR,ZV_QoutR_p,ierr)
call VecRestoreArrayF90(ZV_QoutprevR,ZV_QoutprevR_p,ierr)
call VecRestoreArrayF90(ZV_b,ZV_b_p,ierr)
!end if


!-------------------------------------------------------------------------------
!Calculation of babsmax, QoutRabsmin and QoutRabsmax
!-------------------------------------------------------------------------------
if (BS_opt_influence) then

call VecGetArrayF90(ZV_b,ZV_b_p,ierr)
call VecGetArrayF90(ZV_babsmax,ZV_babsmax_p,ierr)
do JS_localsize=1,IS_localsize
     if (ZV_babsmax_p(JS_localsize)<=abs(ZV_b_p(JS_localsize))) then
         ZV_babsmax_p(JS_localsize) =abs(ZV_b_p(JS_localsize))
     end if
end do
call VecRestoreArrayF90(ZV_b,ZV_b_p,ierr)
call VecRestoreArrayF90(ZV_babsmax,ZV_babsmax_p,ierr)

call VecGetArrayF90(ZV_QoutR,ZV_QoutR_p,ierr)
call VecGetArrayF90(ZV_QoutRabsmin,ZV_QoutRabsmin_p,ierr)
call VecGetArrayF90(ZV_QoutRabsmax,ZV_QoutRabsmax_p,ierr)
do JS_localsize=1,IS_localsize
     if (ZV_QoutRabsmin_p(JS_localsize)>=abs(ZV_QoutR_p(JS_localsize))) then
         ZV_QoutRabsmin_p(JS_localsize) =abs(ZV_QoutR_p(JS_localsize))
     end if
     if (ZV_QoutRabsmax_p(JS_localsize)<=abs(ZV_QoutR_p(JS_localsize))) then
         ZV_QoutRabsmax_p(JS_localsize) =abs(ZV_QoutR_p(JS_localsize))
     end if
end do
call VecRestoreArrayF90(ZV_QoutR,ZV_QoutR_p,ierr)
call VecRestoreArrayF90(ZV_QoutRabsmin,ZV_QoutRabsmin_p,ierr)
call VecRestoreArrayF90(ZV_QoutRabsmax,ZV_QoutRabsmax_p,ierr)

end if
!-------------------------------------------------------------------------------
!Reset previous
!-------------------------------------------------------------------------------
call VecCopy(ZV_QoutR,ZV_QoutprevR,ierr)              !Qoutprev=Qout
!-------------------------------------------------------------------------------
!End temporal loop
!-------------------------------------------------------------------------------
end do

call VecRestoreArrayF90(ZV_C1,ZV_C1_p,ierr)
call VecRestoreArrayF90(ZV_C2,ZV_C2_p,ierr)
call VecRestoreArrayF90(ZV_C3,ZV_C3_p,ierr)
call VecRestoreArrayF90(ZV_Qext,ZV_Qext_p,ierr)

end subroutine routing_params_MC






REAL FUNCTION SHAPE1(bw,z,h)
real :: bw,z,h
real :: sh1,sh2,sh3
        !---trapezoidal channel shape
        sh1=(bw+2*z*h)*(5*bw+6*h*sqrt(1+z**2))
        sh2=4*z*h**2*sqrt(1+z**2)
        sh3=(3*h*(bw+z*h))*(bw+2*h*sqrt(1+z**2))
        if(sh3.eq.0) then
                SHAPE1=0.
        else
                SHAPE1=(sh1+sh2)/sh3
        end if
END FUNCTION SHAPE1


!http://ocw.usu.edu/Biological_and_Irrigation_Engineering/
!Irrigation___Conveyance_Control_Systems/6300__L16_ChannelCrossSections.pdf
REAL FUNCTION FLOW(n,So,bw,h,z)
real :: n,So,bw,z,h
real :: WP,AREA
        WP=bw+2*h*sqrt(1+z**2) !modified from bw+2*h*sqrt(1+h**2); tried
        AREA=(bw+z*h)*h
        if(WP.lt.0) then
                print *,"   !!!****** lpr ERROR IN FLOW: Muskingum-Cunge *******!!!"
                print *,"   !!!****** lpr ERROR: wetted paramter is zero *******!!!"
                print *,'   !!!****** MODEL STOP AT RAPID MUSKINGUM-CUNGE ROUTING ***'
                FLOW=0.00001
                stop
        else
                FLOW=(1/n)*sqrt(So)*(AREA**(5./3.)/(WP**(2./3.)))
        end if

END FUNCTION FLOW


























end module rapid_routing_MC 



