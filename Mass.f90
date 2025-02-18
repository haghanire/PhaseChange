
    SUBROUTINE MASS_DOT2   !!!sharp-interface treatment of interfacial mass flux
    USE SHARE
    IMPLICIT NONE
    INTEGER::INB(Nx)    ! 0: Not INB    1: INB_Phase 1      2: INB_Phase 2       3: INB_Phase 3
    REAL(8),PARAMETER::C_I=0.5d0
    REAL(8),PARAMETER::T_I=TR
    REAL(8)::X_I,Y_I
    
    REAL(8)::Theta_alphaV,T0V,T1V,T2V,qV
    REAL(8)::Theta_alphaL,T0L,T1L,T2L,qL
    REAL(8)::edotnablaT_V  !     temperature gradient at PI on the vapor side
    REAL(8)::edotnablaT_L  !     temperature gradient at PI on the liquid side
    REAL(8)::mdotprime     !     interfacial mass flux
    REAL(8)::r_lwf         !     nondimensional argument in linear weighting function
    REAL(8)::lwf           !     linear weighting function
    REAL(8)::c3(-1:nx+2),c1(-1:nx+2),tem(-1:nx+2),t_sat,t_gas,alpha_B,rho3,h_fg,rho1
    
c1=c
c3=1.0-c1
tem=t
t_sat=tR
t_gas=tL
alpha_B=Cs2*(tau_tv-0.5)
rho1=RhoL
rho3=rhoG
h_fg=hfg
    C3=1.0d0-C1
    INB=0
    DO X=2,Nx-1
        IF (C3(X)==0.5d0) THEN
            INB(X)=4
        ELSE IF ( ((C3(X)-C_I)*(C3(X+1)-C_I))<0.0d0 ) THEN
           IF (C1(X)>0.5d0) INB(X)=1
           IF (C3(X)>0.5d0) INB(X)=3
        ELSE IF ( ((C3(X)-C_I)*(C3(X-1)-C_I))<0.0d0 ) THEN
           IF (C1(X)>0.5d0) INB(X)=1
           IF (C3(X)>0.5d0) INB(X)=3
        END IF
    END DO 
    
    Mass=0.0d0
    
    DO X=2,NX-2
            IF (INB(X)==3 ) THEN
  
            Theta_alphaV=(C3(X)-C_I)/(C3(X)-C3(X+1))
            T1V=Tem(X)
            T2V=Tem(X-1)
            T0V=(2.0d0*T_I+(Theta_alphaV-1.0d0)*T2V)/(1.0d0+Theta_alphaV)
            edotnablaT_V=( (1.0d0+2.0d0*Theta_alphaV)*T0V-4.0d0*Theta_alphaV*T1V-(1.0d0-2.0d0*Theta_alphaV)*T2V )*0.5d0
            
            Theta_alphaL=1.0d0-Theta_alphaV
            T1L=Tem(X+1)
            T2L=Tem(X+2)
            T0L=(2.0d0*T_I+(Theta_alphaL-1.0d0)*T2L)/(1.0d0+Theta_alphaL)
            
            edotnablaT_L=( (1.0d0+2.0d0*Theta_alphaL)*T0L-4.0d0*Theta_alphaL*T1L-(1.0d0-2.0d0*Theta_alphaL)*T2L )*0.5d0
            qV=-alpha_B*Rho3*edotnablaT_V
            qL=-alpha_B*Rho1*edotnablaT_L
            mdotprime=(qV-qL)/H_fg
            
            r_lwf=Theta_alphaV   
            lwf=1.0d0-r_lwf    !!! linear weighting function
            IF (r_lwf>=1.0d0) lwf=0.0d0
            MASS(X)=mdotprime*lwf
            !!!
            r_lwf=1.0d0-Theta_alphaV
            lwf=1.0d0-r_lwf
            IF (r_lwf>=1.0d0) lwf=0.0d0
            MASS(X+1)=mdotprime*lwf 
            END IF
            
            IF (INB(X)==4) THEN !!!  P_I    Backward second order
            T0V=T_I      
            T1V=Tem(X-1)
            T2V=Tem(X-2)
            edotnablaT_V=( 3.0d0*T0V-4.0d0*T1V+T2V )*0.5d0
            
            T0L=T_I      
            T1L=Tem(X+1)
            T2L=Tem(X+2)
            edotnablaT_L=( 3.0d0*T0L-4.0d0*T1L+T2L )*0.5d0
            qV=-alpha_B*Rho3*edotnablaT_V
            qL=-alpha_B*Rho3*edotnablaT_L
            mdotprime=(qV-qL)/H_fg
            
            MASS(X)=mdotprime          
        END IF
       END DO
    
    
    END SUBROUTINE MASS_DOT2