!----------------------------------------------------------------
!------Peirong Lin, April 2017 for Muskingum-Cunge RAPID--------
!----------------------------------------------------------------
module rapid_routing_MC_HAND

IMPLICIT NONE

CONTAINS


!--- update MC routing params every time step----
subroutine routing_params_MC_hand(comid,bw,ss,bs,ll,nchan,hinit,tp,acat,&
                llcat,ncat,slp,tlag,&
                vh,qh,twh,&
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
                speed1 !for speed up certain rivers


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
real, dimension(IS_riv_tot,83),intent(in) :: vh,qh,twh
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
real :: qmax,qmin,vmax,vmin,twmax,twmin,vnow,twnow
integer :: tmp

idebug=1013 !my lucky river !change this value for debugging


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
qnow=ZV_qp(JS_localsize) 
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
!print *,'   --- MC params: bw = ',bw(idebug),' m; ss=',ss(idebug)
!print *,'   --- LPR CHECK: h (after iteration) =',h(idebug), 'm '
print *,'-------------------------------------------------------------------'
end if
call VecRestoreArrayF90(ZV_b,ZV_b_p,ierr)



!**************SEARCH FOR HAND HR **************************************
tmp=1 !---searching Geometry and Velocity matching Qnow---
if(qnow.le.qh(JS_localsize,1)) then
    if(JS_localsize.eq.idebug) then
    print *,'       !!! *** predicted Q <=0 **** STOP!!!'
    stop
    end if
    qmax=qh(JS_localsize,2)
    qmin=qh(JS_localsize,1)    
    vmax=vh(JS_localsize,2)
    vmin=vh(JS_localsize,1)
    twmax=twh(JS_localsize,2)
    twmin=twh(JS_localsize,1)
    !! lpr: might have some problem for the smallest number, check later
    qnow=(qmax-qmin)*0.1 !qmin !(qmax-qmin)*0.1       !if Q=0, assume 10% of lowest value
    vnow=(vmax-vmin)*0.1 !vmin !(vmax-vmin)*0.1       !if Q=0, assume 10% of lowest value
    twnow=(twmax-twmin)*0.1 !twmin !(twmax-twmin)*0.1
else if(qnow.ge.qh(JS_localsize,83)) then
    if(JS_localsize.eq.idebug) then
    print *,'       !!! *** predicted Q > max(handQ) ****!!!'
    end if
    qnow=qh(JS_localsize,83)
    vnow=vh(JS_localsize,83)
    twnow=twh(JS_localsize,83)
else  !Q>0.and.Q<Qmax
!---find V and TW based on HAND V-Q and TW-Q relationship---
    if(JS_localsize.eq.idebug) then
    print *,'   --- NOW SEARCHING V and TW based on HAND -------'
    end if
do tmp = 2,83
if(qnow.le.qh(JS_localsize,tmp)) then
    qmax=qh(JS_localsize,tmp)
    qmin=qh(JS_localsize,tmp-1)
    vmax=vh(JS_localsize,tmp)
    vmin=vh(JS_localsize,tmp-1)
    twmax=twh(JS_localsize,tmp)
    twmin=twh(JS_localsize,tmp-1)
    vnow=(vmax-vmin)*((qnow-qmin)/(qmax-qmin))+vmin !linear interpolation
    twnow=(twmax-twmin)*((qnow-qmin)/(qmax-qmin))+twmin !linear interpolation
    EXIT
else
    CONTINUE
end if
end do
end if

!use NWM channel roughness (divide by 0.05 because Xing used 0.05 everywhere);
!beta is 5/3 for simplicity
!***********************************************************************
Ck(JS_localsize)=vnow*0.05/nchan(JS_localsize)*(5./3.) 
!***********************************************************************
D1= ZV_qp(JS_localsize)/(2*twnow*bs(JS_localsize)) !pg295
x(JS_localsize)=0.5-(D1/(Ck(JS_localsize)*ll(JS_localsize)))  !pg295
!lpr: deal with erraneous places-----
if(x(JS_localsize).le.0..or.isnan(x(JS_localsize))) then
    x(JS_localsize)=0.3 !if <0, set to 0.3
end if
if(Ck(JS_localsize).eq.0. .or. isnan(Ck(JS_localsize))) then
    Ck(JS_localsize) = 0.35*5./3.
end if
if(JS_localsize.eq.idebug) then
print *,'   --- the velocity ranks ',tmp,' in the HAND table ---'
print *,'   --- Tw(now)=',twnow
print *,'   --- CELERITY Ck = ',Ck(JS_localsize)
end if


tch(JS_localsize)=ll(JS_localsize)/Ck(JS_localsize)
Km(JS_localsize)=tch(JS_localsize) !+tov(JS_localsize) !ll(JS_localsize)/Ck(JS_localsize)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
D = Km(JS_localsize)*(1-x(JS_localsize))+ZS_dtR/2
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


end do


!---print out information after C1, C2, C3, C4 are calculated for each channel---
print *,'   --- LPR CHECK: Ck = ',Ck(idebug),' m/s'
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

end subroutine routing_params_MC_hand


end module rapid_routing_MC_HAND



