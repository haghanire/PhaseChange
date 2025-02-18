!**********************************************************************
MODULE SHARE_TEM
    USE SHARE, ONLY:NX
    REAL(8)::T1(-1:Nx+2),T_half(-1:Nx+2),T_halfx(1:Nx-2)
    REAL(8)::rx,ry,dx,dy,dt
    REAL(8)::Axx, Bxx,Cxx,RHSx(1:Nx-2)
    
    END MODULE SHARE_TEM
!!!*******************************************

SUBROUTINE TEMPERATURE_FD2    !!!The sharp interface macroscopic energy equation-Euler method
USE SHARE
USE SHARE_TEM
IMPLICIT NONE

INTEGER::INB(Nx)    ! 0: Not INB    1: INB_Phase 1      2: INB_Phase 2       3: INB_Phase 3      4: PI
REAL(8),PARAMETER::C_I=0.5d0
REAL(8),PARAMETER::T_I=tr
REAL(8)::Theta_alphaV,T0V,T1V,T2V,qV
REAL(8)::Theta_alphaL,T0L,T1L,T2L,qL
REAL(8)::c3(-1:nx+2),c1(-1:nx+2),tem(-1:nx+2),t_sat,t_gas


!print*,'Tsat',T_SAT
c1=c
c3=1.0-c1
tem=t
t_sat=tR
t_gas=tL

    INB=0
    DO X=2,Nx-1
        IF (C3(X)==C_I) THEN
            INB(X)=4
        ELSE IF ( ((C3(X)-C_I)*(C3(X+1)-C_I))<0.0d0 ) THEN
           IF (C1(X)>0.5d0) INB(X)=1
           IF (C3(X)>0.5d0) INB(X)=3
        ELSE IF ( ((C3(X)-C_I)*(C3(X-1)-C_I))<0.0d0 ) THEN
           IF (C1(X)>0.5d0) INB(X)=1
           IF (C3(X)>0.5d0) INB(X)=3
        END IF
    END DO


    do i=2,nx
        if  (c1(i)>c_cutoff) then
            tem(i)=t_sat
        end if
    end do
    
T1=Tem
T1(-1:ic)=T_GAS
T1(Nx)=T_SAT

DO X=1,Nx
    IF (INB(X)==4) T1(X)=T_I
  !  print*,'INB4'
END DO

DO X=1,Nx-1    !!! Nx-2 : Number of equation in x-direction
    
    rx=Cs2*(tau_tv-0.5)

    Axx=rx
    Bxx=1.0d0-2.0d0*rx
    Cxx=rx
END DO
  
!!! x-direction
Do X=2,Nx-1
     RHSx(X-1)= Axx*T1(X-1)+Bxx*T1(X)+Cxx*T1(X+1)-U(X)*0.5d0*(T1(X+1)-T1(X-1))
     
     IF (INB(X)==3 ) THEN      !!! Vapor Phase
       Theta_alphaV=(C3(X)-C_I)/(C3(X)-C3(X+1))
       T1V=T1(X)
       T2V=T1(X-1)
       T0V=(2.0d0*T_I+(Theta_alphaV-1.0d0)*T2V)/(1.0d0+Theta_alphaV)
       RHSx(X-1)= Axx*T1(X-1)+Bxx*T1(X)+Cxx*T0V  -U(X)*0.5d0*(T0V-T1(X-1))
     END IF
     
     IF (INB(X)==1 ) THEN     !!! Liquid Phase
        Theta_alphaV=(C3(X-1)-C_I)/(C3(X-1)-C3(X))
        Theta_alphaL=1.0d0-Theta_alphaV
        T1L=T1(X)
        T2L=T1(X+1)
        T0L=(2.0d0*T_I+(Theta_alphaL-1.0d0)*T2L)/(1.0d0+Theta_alphaL)
        RHSx(X-1)= Axx*T0L+Bxx*T1(X)+Cxx*T1(X+1)  -U(X)*0.5d0*(T1(X+1)-T0L)
      END IF  
END DO
    

    
    Do X=2,Nx-1
        Tem(X)=RHSx(X-1)
    END DO    

Tem(-1:ic)=T_GAS
Tem(Nx)=T_SAT

DO X=1,Nx
    IF (INB(X)==4) Tem(X)=T_I
END DO

    do i=2,nx
        if  (c1(i)>c_cutoff) then
            tem(i)=t_sat
        end if
    end do
    
    
    t=tem
    END 

 