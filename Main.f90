!-------------------------------------- 
!           1D Stefan Problem (D1Q3)
!           Interface Capturing : Conservative Allen-Cahn (Fakhari et al 2015)
!           Interface Capturing : Cahn-Hilliard (Safari et al 2013)    
!           Hydrodynamics       :   
!--------------------------------------     
    
    MODULE SHARE

    PARAMETER ( nx = 300 )

    REAL(8):: h(0:2,0:nx+1),heq(0:2,0:nx+1),g(0:2,0:nx+1),geq(0:2,0:nx+1)
    REAL(8):: gama(0:2,0:nx+1)
    REAL(8):: htemp(0:2,0:nx+1),gtemp(0:2,0:nx+1)
      
    REAL(8):: vis(1:nx),tau(1:nx)
    REAL(8):: u(1:nx),T(-1:nx+2)

    REAL(8):: c(-1:nx+2),D2c(1:nx),D1Cc(1:nx),D1Mc(1:nx)
    REAL(8):: muphi(-1:nx+2),D2muphi(1:nx),D1Cmuphi(1:nx), D1Mmuphi(0:nx)
    REAL(8):: p(-1:nx+2), D1CP(1:nx),D1MP(1:nx)
    REAL(8):: Rho(-1:nx+2), D1CRho(1:nx), D1MRho(1:nx), D2Rho(1:nx)
      
    REAL(8):: SourceG(0:nx+1),SourceH(0:nx+1)
    REAL(8):: C_S,muG,muL,tc,R,Z,xcen,ycen,err
    REAL(8), PARAMETER :: tauL=0.5d0,tauG=0.5d0
    REAL(8), PARAMETER :: D = 6.d0
    REAL(8), PARAMETER :: sigma = 0.0001d0
    REAL(8), PARAMETER :: beta = 12.d0*sigma / D
    REAL(8), PARAMETER :: kapa = 1.5*D*sigma
      
    INTEGER, PARAMETER ::e(0:2)   =[0, 1, -1]
    REAL(8), PARAMETER:: Wa(0:2)=[2.d0/3.d0, 1.d0/6.d0, 1.d0/6.d0]
      
    REAL(8), PARAMETER :: RhoL = 1.0
    REAL(8), PARAMETER :: RhoG =1.0d0
    REAL(8), PARAMETER :: visL = 0.01d0
    REAL(8), PARAMETER :: visG = 0.01d0
	REAL(8), PARAMETER :: Cs2 = 1.d0/3.d0
	INTEGER, PARAMETER :: ic =15
	INTEGER, PARAMETER :: tfinal = 3000000

	REAL(8), PARAMETER :: m =0.01   !!! 0.0002d0/beta
    REAL(8), PARAMETER::MBETA=0.01d0
    REAL(8), PARAMETER::Mobility=MBETA/Beta
    
	REAL(8), PARAMETER :: Ja = 0.02d0
	REAL(8), PARAMETER :: Pe = 100.
	
    
    REAL(8), PARAMETER :: tau_tl = 1.0
    REAL(8), PARAMETER :: tau_tv = 1.0d0 ! 3.0d0*0.001d0+0.5d0
    REAL(8), PARAMETER :: alfa = 2.5d-5 ! (tau_tv-0.5)/3.0
    
	REAL(8), PARAMETER :: cx = 0.0001! 5d-5 !5d-5
	REAL(8), PARAMETER :: ct =  (tau_tv-0.5)*cx**2/(3*alfa)   !1d-4/15.
	REAL(8), PARAMETER :: hfg = 10.d3*(ct/cx)**2 !1.0/ Ja !10.d3*(ct/cx)**2
	REAL(8), PARAMETER :: tr = 0.0d0 !hfg*Ja !hfg*Ja
    REAL(8), PARAMETER :: tl = (tr + hfg*Ja) !1.0d0
    REAL(8), PARAMETER :: c_cutoff = 0.5d0
	REAL(8), PARAMETER :: dt=1.d0

    INTEGER:: delx,dely,i,j,k,cons,REAL_intloc,x
    REAL(8)::intloc
    REAL(8)::tau_t,u_interface,mass(-1:nx+1)

    REAL(8)::C1_P(0:2,1:Nx)
    REAL(8)::TForceX(1:Nx)
    ENDmodule SHARE

!********************************************************************************************

PROGRAM main

    USE SHARE
	IMPLICIT NONE

	INTEGER :: tt,tt1,ii,jj,loc,loc_old,swap,count
	REAL :: rn,time
	REAL(8) ::  a1,sum,sum_num,umax,sumRho
	LOGICAL :: written

    PRINT*, 'Hfg=',hfg
    PRINT*, 'cx=',cx
    PRINT*, 'ct=',ct
    
    OPEN (10,file='result.plt')
    OPEN (11,file='x interface.plt')
    OPEN (12,file='velocity interface.plt')
    
	PRINT*, sqrt(2*kapa*beta)/6. , beta,ic
	PRINT*, 'cx = ', cx
	PRINT*, 'ct = ', ct
	PRINT*, 'cu = ', cx/ct
	PRINT*, 'Total length = ', nx*cx
	PRINT*, 'h_fg = ', hfg
	PRINT*, 'tl = ', Tl
    PRINT*, 'alpha*Rho*Delta T/Hfg=', ((tau_tv-0.5)/3.0d0)*Rhog*(tl-tr)/hfg
	CALL initialize

	err = 0.1
	tt = 0

  !  a1 =0.6201  !Ja=1.0
   a1 =0.0997 !Ja=0.02
    
 !   a1 = 0.706*sqrt(Ja)*(1-0.21*((0.506*Ja)**(0.93-1.5*Ja)))
    
    DO tt = 1, 10000   

        IF (ISNAN(sumRho)) PAUSE
    CALL CAC (0)
  !  CALL CH (0)
  !  PRINT*, maxval(u), sum(rho(0:nx)),tt
 	ENDDO
	 

   tt=0
    DO while  (loc_old<250+ic)  ! tt=0,1500000
     !   print*,'ttt=',tt
        IF (ISNAN(sumRho)) PAUSE
        
     CALL CAC (1)
  ! CALL CH (1)
    
        umax=0.0
        sumRho=0
		
		DO i=0,nx
        IF (abs(u(i))>umax) then 
        umax=abs(u(i))
        ENDIF 
        sumRho=sumRho+Rho(i)
        ENDDO
        
	    CALL interfaceloc

        IF (MOD(tt,100)==0) then
				 
        PRINT*, umax,"  ",sumRho,"  ",tt,"   ",intloc

    !    IF ((real_intloc.ne.loc_old)) then
            WRITE(11,*) tt*ct,(intloc-ic)*cx,2*a1*sqrt(alfa*ct*tt) 
            WRITE(12,*) tt*ct,u(250)*cx/ct,(1.0-(RhoG/RhoL))*a1*sqrt(alfa/(ct*tt+1e-32)) 
     !   ENDIF 
                
        loc_old = real_intloc         
        ENDIF 
        
        IF (MOD(tt,100000)==0) CALL results

    tt=tt+1
	ENDDO 

    END PROGRAM main
    
    !*********************************************************************************************
    SUBROUTINE CAC (inval)
    USE SHARE
    IMPLICIT NONE
    INTEGER,INTENT(IN)::inval
    
    CALL COLLISION_H_CAC !S
    CALL COLLISION_G_CAC           
    CALL STREAMING
    CALL boundary
    CALL MACROSCOPIC_H
    CALL MACROSCOPIC_G_CAC (inval)
    
    END SUBROUTINE
    !*********************************************************************************************
    SUBROUTINE CH (inval)
    USE SHARE
    IMPLICIT NONE
    INTEGER,INTENT(IN)::inval
    
    CALL COLLISION_H_CH
    CALL COLLISION_G_CH           
    CALL STREAMING
    CALL boundary
    CALL MACROSCOPIC_H
    CALL MACROSCOPIC_G_CH (inval)
    
    END SUBROUTINE
    !*********************************************************************************************

 SUBROUTINE initialize
    USE SHARE
    IMPLICIT NONE

     DO i=-1,nx+2
        c(i) = 0.5d0+0.5d0*dtanh(2.0d0*REAL(i-ic)/D)
        Rho(i)=RhoG+c(i)*(RhoL-RhoG)
        p(i)=0.d0
        u(i)=0.d0
     ENDDO
     
     DO i=0,nx+1
          h(0:2,i)=Wa(0:2)*c(i)
          g(0:2,i)=0.d0
          END DO

   DO i=-1,ic-1
		T(i) = Tl  !((Tr-Tl)/ic)*i+Tl
	ENDDO

    DO i=ic,nx+1
		T(i) = Tr
    ENDDO

    CALL laplace(Rho,D2Rho)
    CALL gradient_C(Rho,D1CRho)
    CALL gradient_M(Rho,D1MRho)
    
    CALL laplace(c,D2c)
    CALL gradient_C(c,D1Cc)
    CALL gradient_M(c,D1Mc)
    
    CALL POTENTIAL
    
    DO i=1,Nx
	vis(i)=visG+c(i)*(visL-visG)
	tau(i)=0.5d0 !!!3.0d0*vis(:)/Rho(:)
	tau_t=1.0d0
    END DO
    C1_P=0.0d0
    CALL EQULIBRIUM 
    ENDSUBROUTINE initialize
!*********************************************************************************************
 SUBROUTINE COLLISION_H_CAC
    USE SHARE
    IMPLICIT NONE
    REAL(8)::tmp1,Temp1H(0:2), Temp2H(0:2)

    CALL EQULIBRIUM_H_CAC
    
    DO i=1,nx
        SourceH(i) = -mass(i) / (RhoL)
    ENDDO

    DO i=1,nx
    DO k = 0,2
        
        tmp1=3.0d0*M/(3.0d0*M+0.5d0)
        Temp1H(k)=-wa(k)*(1.0d0+3.0d0*tmp1*( dble(e(k))*U(i) ))*(Mass(i)/RhoL)
        Temp2H(k)=Temp1H(k)-C1_P(k,i)
        C1_P(k,i)=Temp1H(k)
      
        
        h(k,i) = h(k,i) -(h(k,i)-heq(k,i))/(3.0d0*M+0.5d0) +temp1H(k)+0.5d0*temp2H(k)
    ENDDO
    ENDDO

    END SUBROUTINE
!*********************************************************************************************
 SUBROUTINE COLLISION_H_CH
    USE SHARE
    IMPLICIT NONE
    REAL(8)::tmp1,tmp2,temp1H(0:2), temp2H(0:2), temp3H(0:2), FsX, tempfor(0:2)
    REAL(8)::EdotNablarho(0:2),EdotNablaC(0:2),EdotNablaP(0:2)
    
    CALL EQULIBRIUM_H_CH           
    
    DO i=1,nx
    DO k = 0,2

    SourceH(i) = -mass(i) / (RhoL)    
    EdotNablaRho(k)=0.25d0*(-Rho(i+2*e(k)) +5.0d0*Rho(i+e(k))-3.0d0*Rho(i)-Rho(i-e(k)) )
    EdotNablaC  (k)=0.25d0*(-  c(i+2*e(k)) +5.0d0*  c(i+e(k))-3.0d0*  c(i)-  c(i-e(k)) )
    EdotNablap  (k)=0.25d0*(-  p(i+2*e(k)) +5.0d0*  p(i+e(k))-3.0d0*  p(i)-  p(i-e(k)) )
    
    FsX=muphi(i)*D1Mc(i)
    tempfor(k)=muphi(i)*EdotNablaC(k)
    
    temp1H(k)=( EdotNablaC(k) -3.0d0*(C(i)/Rho(i))*(EdotNablap(k)-tempfor(k)) )*gama(k,i)
    temp2H(k)=( D1Mc(i) -3.0d0*(C(i)/Rho(i))*(D1Mp(i)-FsX) )*gama(k,i)*u(i)
    
    
    tmp1=0.0d0
    tmp2=0.0d0
    
    IF ( i+e(k)==0 .or. i+e(k)==nx+1 ) then
        tmp1=D2muphi(i)*gama(k,i)*Mobility
        tmp2=SourceH(i)*gama(k,i)
    ELSE
        tmp1=D2muphi(i+e(k))*gama(k,i+e(k))*Mobility
        tmp2=SourceH(i+e(k))*gama(k,i+e(k))
    ENDIF
        
        temp3H(k)=(D2muphi(i)*gama(k,i)*Mobility+SourceH(i)*gama(k,i))*0.5d0+(tmp1+tmp2)*0.5d0
        h(k,i) = heq(k,i)+temp1H(k)-temp2H(k)+temp3H(k)

    ENDDO
    ENDDO

    END SUBROUTINE
    
 !*********************************************************************************************
 SUBROUTINE COLLISION_G_CH 
    USE SHARE
    IMPLICIT NONE
    REAL(8)::tmp1,tmp2, temp1G(0:2), temp2G(0:2), temp3G(0:2), tempfor(0:2),FsX
    REAL(8)::EdotNablarho(0:2),EdotNablaC(0:2),EdotNablaP(0:2)


    CALL EQULIBRIUM_G
    DO i=1,nx
        SourceG(i) = mass(i)*(1.d0/RhoG -1.d0/RhoL)   ! Sourceterm for momentum
    ENDDO   

    
    DO i=1,nx
    DO k = 0,2
        EdotNablaRho(k)=0.25d0*(-Rho(i+2*e(k)) +5.0d0*Rho(i+e(k))-3.0d0*Rho(i)-Rho(i-e(k)) )
        EdotNablaC  (k)=0.25d0*(-  C(i+2*e(k)) +5.0d0*  C(i+e(k))-3.0d0*  C(i)- C (i-e(k)) )
        FsX=muphi(i)*D1Mc(i)
        tempfor(k)=muphi(i)*EdotNablaC(k)
        
        temp1G(k)=((gama(k,i)-Wa(k))*EdotNablaRho(k)/3.0d0+tempfor(k)*gama(k,i))
        temp2G(k)=((gama(k,i)-Wa(k))*D1Mrho(i)/3.0d0+FsX*gama(k,i))*u(i)
        temp3G(k)=Cs2*Rho(i)*SourceG(i)*Wa(k)
        g(k,i) = g(k,i)-(g(k,i)-geq(k,i))/(tau(i)+0.5D0)+temp1G(k)-temp2G(k)+temp3G(k)
    ENDDO
    ENDDO
    
    END
 !*********************************************************************************************
 SUBROUTINE COLLISION_G_CAC
    USE SHARE
    IMPLICIT NONE
    REAL(8)::tmp1,tmp2, temp1G(0:2), temp2G(0:2), temp3G(0:2), tempfor(0:2),FsX
    REAL(8)::EdotNablarho(0:2),EdotNablaC(0:2),EdotNablaP(0:2)


    !CALL EQULIBRIUM_G
    DO i=1,nx
        SourceG(i) = mass(i)*(1.d0/RhoG -1.d0/RhoL)   ! Sourceterm for momentum
    ENDDO   

    
    DO i=1,nx
    DO k = 0,2
        FsX=muphi(i)*D1cc(i)
        temp1G(k)=(gama(k,i)-Wa(k))+p(i)*Wa(k)
        Temp2G(k)=Wa(k)*3.0d0*( E(k)*TForcex(i) )/Rho(i)  
        Geq(k,i)=Temp1G(k)-0.5d0* Temp2G(k) 
        
        g(k,i) = g(k,i)-(g(k,i)-geq(k,i))/(tau(i)+0.5D0)+Temp2G(k)+SourceG(i)*Wa(k)
    ENDDO
    ENDDO
    
    END
    
!**********************************************************************    
SUBROUTINE EQULIBRIUM_H_CAC
    USE SHARE
    IMPLICIT NONE

    DO k = 0,2
    DO i=1,nx
    heq(k,i) = c(i)*gama(k,i)+Wa(k)*12.d0*M*c(i)*(1.d0-c(i))*e(k)/D              
    ENDDO
    ENDDO

    END SUBROUTINE EQULIBRIUM_H_CAC
    

!**********************************************************************    
SUBROUTINE EQULIBRIUM_H_CH
    USE SHARE
    IMPLICIT NONE
    REAL(8)::temp1H(0:2), temp2H(0:2), temp3H(0:2)
    REAL(8)::EdotNablaRho(0:2), EdotNablac(0:2), EdotNablap(0:2), tempfor(0:2), FsX
               
    
    DO i=1,nx
    DO k = 0,2
    EdotNablaRho(k)=0.5d0*(Rho(i+e(k))-Rho(i-e(k)) )
    EdotNablaC  (k)=0.5d0*(C  (i+e(k))- C (i-e(k)) )
    EdotNablap  (k)=0.5d0*(p  (i+e(k))- p (i-e(k)) )
    
    FsX=muphi(i)*D1Cc(i)
    tempfor(k)=muphi(i)*EdotNablaC(k)
    
    temp1H(k)=( EdotNablaC(k) -3.0d0*(C(i)/Rho(i))*(EdotNablap(k)-tempfor(k)) )*gama(k,i)
    temp2H(k)=( D1Cc(i) -3.0d0*(C(i)/Rho(i))*(D1Cp(i)-FsX) )*gama(k,i)*u(i)
    
    heq(k,i)=c(i)*gama(k,i)-0.5d0*(temp1H(k)-temp2H(k))             
    ENDDO
    ENDDO

    END SUBROUTINE EQULIBRIUM_H_CH 
    
!**********************************************************************    
 
SUBROUTINE EQULIBRIUM_G
    USE SHARE
    IMPLICIT NONE
    REAL(8)::tmp1,tmp2,temp1G(0:2),temp2G(0:2),temp3G(0:2),FsX, tempfor(0:2)
    REAL(8)::EdotNablarho(0:2),EdotNablaC(0:2),EdotNablaP(0:2)
    
    DO i=1,nx
    SourceG(i) = mass(i)*(1.d0/RhoG -1.d0/RhoL)   ! Sourceterm for momentum
    ENDDO
    
    DO i=1,nx
    DO k = 0,2

    EdotNablaRho(k)=0.5d0*(Rho(i+e(k))-Rho(i-e(k)) )
    EdotNablaC  (k)=0.5d0*(C  (i+e(k))- C (i-e(k)) )
    FsX=muphi(i)*D1Cc(i)
    tempfor(k)=muphi(i)*EdotNablaC(k)
  
    temp1G(k)=( (gama(k,i)-Wa(k))*EdotNablaRho(k)/3.0d0  + tempfor(k)*gama(k,i) )
    temp2G(k)=( (gama(k,i)-Wa(k))*D1crho(i)/3.0d0  + FsX*gama(k,i) )*u(i)
    temp3G(k)=Cs2*Rho(i)*SourceG(i)*Wa(k)
    geq(k,i) = Rho(i)*(1.d0/3.d0)*gama(k,i)+(P(i)-Rho(i)/3.d0)*Wa(k)-0.5d0*(temp1G(k)-temp2G(k))-0.5d0*temp3G(k)
    ENDDO
    ENDDO

END SUBROUTINE EQULIBRIUM_G   
    
 !**********************************************************************  
       
 SUBROUTINE STREAMING
    USE SHARE
    IMPLICIT NONE
    INTEGER::xe,xw
    
    DO i=1,NX
    Xe=i+1
    Xw=i-1
    IF (i==NX) Xe=0
    IF (i==0) Xw=NX

    htemp(0 ,i) = h(0,i)
    htemp(1,xe) = h(1,i)
    htemp(2,xw) = h(2,i)
    
    gtemp(0 ,i) = g(0,i)
    gtemp(1,xe) = g(1,i)
    gtemp(2,xw) = g(2,i)

    ENDDO
    h=htemp
    g=gtemp

    ENDSUBROUTINE STREAMING
!**********************************************************************

    SUBROUTINE MACROSCOPIC_H
    USE SHARE
    IMPLICIT NONE


    DO i=1,nx
    c(i)=h(0,i)+h(1,i)+h(2,i)
    Rho(i)=RhoG+c(i)*(RhoL-RhoG)
    ENDDO
    c(0)=c(2) 
    c(-1)=c(3)              
    c(nx+1)=2*c(nx)-c(nx-1) 
    c(nx+2)=2*c(nx+1)-c(nx)
    Rho( 0)=RhoG+c( 0)*(RhoL-RhoG)
    Rho(-1)=RhoG+c(-1)*(RhoL-RhoG)
    Rho(nx+1)=RhoG+c(nx+1)*(RhoL-RhoG)
    Rho(nx+2)=RhoG+c(nx+2)*(RhoL-RhoG)
    
    
    CALL laplace(c,D2c)
    CALL gradient_C(c,D1Cc)
    CALL gradient_M(c,D1Mc)
    
    CALL laplace(Rho,D2Rho)
    CALL gradient_C(Rho,D1CRho)
    CALL gradient_M(Rho,D1MRho)
    
    CALL POTENTIAL
    
    DO i=1,Nx
	vis(i)=visG+c(i)*(visL-visG)
	tau(i)=0.5d0 !!!3.0d0*vis(:)/Rho(:)
	tau_t=1.0d0
    END DO

    ENDSUBROUTINE MACROSCOPIC_H
    
    !**********************************************************************************************
    SUBROUTINE POTENTIAL 
    USE SHARE
    IMPLICIT NONE
    
    DO i=1,nx
    muphi(i)=4.d0*beta*c(i)*(c(i)-1.0d0)*(c(i)-0.5d0)-kapa*D2c(i)
    ENDDO
    
    muphi(0)=muphi(2)
    muphi(-1)=muphi(3)
    muphi(nx+1)=2*muphi(nx)-muphi(nx-1) 
    muphi(nx+2)=2*muphi(nx+1)-muphi(nx) 
    
    CALL LAPLACE(muphi,D2muphi)    
    
    END 

     !**********************************************************************************************
    SUBROUTINE EQULIBRIUM 
    USE SHARE
    IMPLICIT NONE   
    
    DO i=1,nx
    DO k = 0,2
    gama(k,i) = Wa(k)*(1.0D0+3.0D0*(e(k)*u(i))+4.5D0*(e(k)*u(i))**2-1.5D0*(u(i)**2))
    ENDDO
    ENDDO
    
    END 
    !**********************************************************************************************
    SUBROUTINE MACROSCOPIC_G_CH (inval)
	USE SHARE
    IMPLICIT NONE
    INTEGER,intent(in)::inval
    
    DO i=1,nx
    SourceG(i) = mass(i)*inval*(1.d0/RhoG -1.d0/RhoL)   ! Sourceterm for momentum
    SourceH(i) = -mass(i)*inval / (RhoL)
    ENDDO
    
    DO i=1,nx
    u(i)=(1.d0/Rho(i))*((g(1,i)-g(2,i))/Cs2+muphi(i)*D1Cc(i)/2.d0)
    p(i)=g(0,i)+g(1,i)+g(2,i)+(u(i)*D1CC(i)*(RhoL-RhoG)*Cs2)/2.d0 + Cs2*Rho(i)*SourceG(i)*0.5d0
    ENDDO
    
    p (0)=p(2)
    p(-1)=p(3)
    p(nx+1)=2*p(nx)-p(nx-1) 
    p(nx+2)=2*p(nx+1)-p(nx) 
    
    CALL GRADIENT_C(p,D1Cp)  
    CALL GRADIENT_M(p,D1Mp) 
               
    CALL EQULIBRIUM
    CALL TEMPERATURE_FD2
    Mass=0.0d0
    IF (inval==1) CALL  MASS_DOT2 ! mass is + for evaporation  
    
    ENDSUBROUTINE  

    !**********************************************************************************************
    SUBROUTINE MACROSCOPIC_G_CAC (inval)
	USE SHARE
    IMPLICIT NONE
    INTEGER,intent(in)::inval
    REAL(8)::Fp, Fs, Fmu,tmp1, Temp1G(0:2), Temp2G(0:2)
    
    DO i=1,nx
    SourceG(i) = mass(i)*inval*(1.d0/RhoG -1.d0/RhoL)   ! Sourceterm for momentum
    SourceH(i) = -mass(i)*inval / (RhoL)
    ENDDO
    
    DO i=1,nx
        Fs=muphi(i)*D1Cc(i)
        p(i)=g(0,i)+g(1,i)+g(2,i)   !+ mass(i)*inval*(1.d0/RhoG -1.d0/RhoL)   !p*
        Fp=-p(i)*D1CRho(i)/3.0d0    !p*
        Temp1G(0:2)=P(i)*wa(0:2)+(gama(0:2,i)-Wa(0:2))  !geq
        tmp1=e(1)*e(1)*(g(1,i)-Temp1G(1))+e(2)*e(2)*(g(2,i)-Temp1G(2))
        Fmu=-(tau(i)/(tau(i)+0.5))*tmp1*D1CRho(i)
        TForceX(i)=Fs+Fp+Fmu
        u(i)=g(1,i)-g(2,i)+0.5d0*TForceX(i)/Rho(i)
    ENDDO
    
    p (0)=p(2)
    p(-1)=p(3)
    p(nx+1)=2*p(nx)-p(nx-1) 
    p(nx+2)=2*p(nx+1)-p(nx) 
    
    CALL GRADIENT_C(p,D1Cp)  
    CALL GRADIENT_M(p,D1Mp) 
               
    CALL EQULIBRIUM
    CALL TEMPERATURE_FD2
    Mass=0.0d0
    IF (inval==1) CALL  MASS_DOT2 ! mass is + for evaporation  
    
    ENDSUBROUTINE 
    
!*********************************************************************************************

SUBROUTINE boundary
    USE SHARE
    IMPLICIT NONE

    !Left Boundary
    h(1,1)=h(2,1)
    g(1,1)=g(2,1)
    !Right Boundary
	h(2,nx)=h(2,nx-1)
    g(2,nx)=g(2,nx-1)   

    ENDSUBROUTINE boundary

!****************************************************************************

SUBROUTINE GRADIENT_C(B,D1C)
	USE SHARE
    IMPLICIT NONE
	REAL(8),INTENT(IN) ::B (-1:nx+2)
	REAL(8),INTENT(OUT)::D1C(1:nx)

    DO i=1,nx

	D1C(i) = ( B(i+1)- B(i-1) ) / 2.d0

	ENDDO

    END
!**********************************************************************

SUBROUTINE GRADIENT_M(B,D1C)
	USE SHARE 
	IMPLICIT NONE
	REAL(8),INTENT(IN) ::B (-1:Nx+2)
	REAL(8),INTENT(OUT)::D1C(1:nx)
	REAL(8)::tempCD,tempBD
    INTEGER::Xe, Xee, Xw, Xww
    
    DO x=1,nx	
    Xe=X+1
    Xee=X+2
    Xw=X-1
    Xww=X-2

		  
    ! X
	tempCD= (B(Xe)-B(Xw))/2.d0 
	tempBD= (-B(Xee)+B(Xww))/4.0d0+B(Xe)-B(Xw)

	D1C(x)=(TempCD+TempBD)/2.0d0 
		   
    END DO

END    
!*********************************************************************************************

SUBROUTINE LAPLACE (B,D2)
	USE SHARE

	REAL(8),INTENT(IN) ::B (-1:nx+2)
	REAL(8),INTENT(OUT)::D2(1:nx)

	DO i=1,nx
	    D2(i)=B(i+1)+B(i-1)-2*B(i)
	ENDDO
END

!*********************************************************************************************

SUBROUTINE interfaceloc
    USE SHARE
    IMPLICIT NONE

		DO i=1,nx-1
			IF ((c(i+1)>0.5).and.(c(i)<0.5)) then
				intloc = i+(0.5-c(i))/(c(i+1)-c(i))
				REAL_intloc = i + 1 
			ENDIF
		ENDDO

ENDSUBROUTINE interfaceloc

!**********************************************************************

SUBROUTINE RESULTS
    USE SHARE
    IMPLICIT NONE
    INTEGER::Ny=50
    REAL(8)::tempvel=0.0d0

    WRITE(10,*)'Variables = X, Y, Ux, Uy, C1, Ph, Tem, Mass'
	WRITE(10,*) 'Zone I=',Nx,',J=',Ny,',F=Point'
     
    DO j=1,Ny
    DO i=1,Nx
    WRITE(10,*) i, j,  u(i), tempvel, c(i), P(i), T(i), Mass(i)
    ENDDO
    ENDDO
    
END  