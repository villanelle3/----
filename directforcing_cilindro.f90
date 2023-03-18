module directforcing
use constantes
use geometria

real(pr),parameter,dimension(6)::alfa=(/ 0.0_pr, -0.691750960670_pr, -1.727127405211_pr, -0.694890150986_pr,  &
&  -1.039942756197_pr, -1.531977447611_pr /), beta=(/ 0.122_pr, 0.477263056358_pr, 0.381941220320_pr, &
&   0.447757195744_pr, 0.498614246822_pr, 0.186648570846_pr /), c=(/ 0.0_pr, 0.122_pr, 0.269115878630_pr, & 
&   0.447717183551_pr, 0.749979795490_pr, 0.898555413085_pr /) !Allampalli 6 estagios

integer(pint)::ii

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! DIRECTFORCING - (Wang et al.,2008)	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine forcaDF
use constantes
use geometria
implicit none

complex(pr),dimension(nx,ny)::fxd,fyd,auxfxd,auxfyd,F0x,F0y,ZBx,ZBy
real(pr)::epslon,sfx,sfy,Cdc,Clc,ruido(ny),theta,gama,u_rand(nx,ny),v_rand(nx,ny),maxfx0,maxfy0,maxfx,maxfy,Cdc2,Clc2
integer(pint)::it,q
integer::nseedx,nseedy
real::rand2


call ZFFT2D(fx,NX,NY,0)
call ZFFT2D(fy,NX,NY,0)

fx=0.0_pr
fy=0.0_pr
front(:)%SFxl=0.0_pr
front(:)%SFyl=0.0_pr
front(:)%STT=0.0_pr

sfx=0.0_pr
sfy=0.0_pr

F0x=0.0_pr
F0y=0.0_pr
maxfx=1.0_pr
maxfx0=1.0_pr
maxfy=1.0_pr
maxfy0=1.0_pr
  
epslon=1.0d-4
it=0

ZBx=phi*(ufis-Qtx) 		!Zona de amortecimentao (Buffer Zone)
ZBy=phi*(vfis-Qty)		!Zona de amortecimentao (Buffer Zone)
call ZFFT2D(ZBx,NX,NY,-1)
call ZFFT2D(ZBy,NX,NY,-1)
call projecao(ZBx,ZBy)
u=u-ZBx
v=v-ZBy
	
do while ((((maxfx/maxfx0)>=epslon).or.((maxfy/maxfy0)>=epslon)).and.(it<nit))
	fxd=0.0_pr
	fyd=0.0_pr
	DO ifi=1,nfi
		call interpolacao
		!calculo da forca lagrangiana 
		front(ifi)%Fxl=(front(ifi)%ulf-front(ifi)%ul)
		front(ifi)%Fyl=(front(ifi)%vlf-front(ifi)%vl)
		
		!calculo do ultimo ponto da geometria
		front(ifi)%Fxl(front(ifi)%nl)=front(ifi)%Fxl(1)
		front(ifi)%Fyl(front(ifi)%nl)=front(ifi)%Fyl(1)

		!calculo da forcas geradas pelo fluido - arrasto, sustentacao, fluido-estrutura
		front(ifi)%SFxl=front(ifi)%SFxl+sum(front(ifi)%Fxl(1:front(ifi)%nl-1))
		front(ifi)%SFyl=front(ifi)%SFyl+sum(front(ifi)%Fyl(1:front(ifi)%nl-1))
		
		call distribuicao(fxd,fyd)
	ENDDO 

!******************************************************************************************!
!*************** 	Condicoes de contorno da regiao de forcagem	*******************!
!******************************************************************************************!

 	do i=iimin-zp,iimin 
 		do j=jjmin,jjmax	
 			fxd(i,j)=real(umax-ufis(i,j))
			fyd(i,j)=real(0.0_pr-vfis(i,j))
 		enddo
 	enddo
! 
    do i=64,nx 
    	fxd(i,64)=real(0.0_pr-ufis(i,64))     
        fyd(i,64)=real(0.0_pr-vfis(i,64))
     enddo 
  
! open(13,file='FORCA.dat',status='unknown')
! ! 1	format(3(1x,E14.7))
! 	write(13,*)'variables= "x","y","fx","fy"'
!	write(13,*)'zone i=',nx,' j=',ny
! 	do j=1,ny
! 		do i=1,nx
 !			write(13,*) x(i),y(j),real(fxd(i,j)),real(fyd(i,j))
 !		enddo
 !	enddo
 !close(13)

	maxfx=maxval(real(fxd)-real(F0x))
	maxfy=maxval(real(fyd)-real(F0y))
	if ((it==1).and.(nstep>=1)) then
		maxfx0=maxfx
		maxfy0=maxfy
	endif
	F0x=real(fxd)
	F0y=real(fyd)
	it=it+1
	
!******************************************************************************************!
!***************    imposicao da forca euleriana sobre o campo estimado *******************!
!******************************************************************************************!

	call ZFFT2D(fxd,NX,NY,-1)
	call ZFFT2D(fyd,NX,NY,-1)
	fx=fx+fxd
	fy=fy+fyd
	call projecao(fxd,fyd)

	u=u+fxd
	v=v+fyd
	
	ufis=u
	vfis=v
	call ZFFT2D(ufis,NX,NY,1)
	call ZFFT2D(vfis,NX,NY,1)
	ufis=real(ufis)
	vfis=real(vfis)
end do 
fx=fx+ZBx
fy=fy+ZBy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! CALCULO DO COEFICIENTE DE ARRASTO E SUSTENTAÇAO	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if ((t>=tmin-0.05d0*dt).and.(t<=tmin+0.05d0*dt)) then
	!if(mod(nstep,ts)==0) then
	do ifi=1,nfi
		front(ifi)%Cd= -2.0_pr*front(ifi)%SFxl*rho/dt*front(ifi)%ds*dx/(rho*umax**2.0_pr*front(ifi)%raio)
		front(ifi)%Cl= -2.0_pr*front(ifi)%SFyl*rho/dt*front(ifi)%ds*dy/(rho*umax**2.0_pr*front(ifi)%raio)
	enddo
	
	Cdc=front(1)%Cd
	Clc=front(1)%Cl
	
	!Cdc2=front(2)%Cd
	!Clc2=front(2)%Cl

	write(22,*) t, Cdc
	write(23,*) t, Clc
	!write(24,*) t, Cdc, Clc
	!write(25,*) t, Cdc2
	!write(26,*) t, Clc2
	
	
endif

end subroutine forcaDF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! INTERPOLACAO DIRECTFORCING	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine interpolacao
use constantes
implicit none
integer(pint):: npoints,npoints2,q,ii,jj,ir,jr
real(pr):: acumx,acumy

front(ifi)%ul=0.0_pr
front(ifi)%vl=0.0_pr

npoints=4
npoints2=npoints/2

do q=1,front(ifi)%nl-1
	ir=int(front(ifi)%xl(q)/dx)
	jr=int(front(ifi)%yl(q)/dy)
	acumx=0.0_pr
	acumy=0.0_pr
	do jj=1,npoints
		do ii=1,npoints
			acumx=Dij(front(ifi)%xl(q),front(ifi)%yl(q),ii-npoints2+(ir),jj-npoints2+(jr))*ufis(ii-npoints2+(ir),jj-npoints2+(jr))*dx*dy
			acumy=Dij(front(ifi)%xl(q),front(ifi)%yl(q),ii-npoints2+(ir),jj-npoints2+(jr))*vfis(ii-npoints2+(ir),jj-npoints2+(jr))*dx*dy
			front(ifi)%ul(q)=front(ifi)%ul(q)+acumx
			front(ifi)%vl(q)=front(ifi)%vl(q)+acumy
		enddo
	enddo
end do

front(ifi)%ul(front(ifi)%nl)=front(ifi)%ul(1)
front(ifi)%vl(front(ifi)%nl)=front(ifi)%vl(1)
end subroutine interpolacao

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!	  DISTRIBUICAO  DIRECTFORCING 	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine distribuicao(fxd,fyd)
use constantes
implicit none
integer(pint):: npoints,npoints2,q,ii,jj,ir,jr
real(pr):: acumx,acumy
complex(pr),dimension(nx,ny)::fxd,fyd

npoints=4
npoints2=npoints/2
do q=1,front(ifi)%nl-1
	ir=int(front(ifi)%xl(q)/dx)
 	jr=int(front(ifi)%yl(q)/dy)
 	acumx=0.0_pr
 	acumy=0.0_pr
	do jj=1,npoints
		do ii=1,npoints
			acumx=Dij(front(ifi)%xl(q),front(ifi)%yl(q),ii-npoints2+(ir),jj-npoints2+(jr))*front(ifi)%Fxl(q)*dx*front(ifi)%ds
			acumy=Dij(front(ifi)%xl(q),front(ifi)%yl(q),ii-npoints2+(ir),jj-npoints2+(jr))*front(ifi)%Fyl(q)*dx*front(ifi)%ds

			fxd(ii-npoints2+int(ir),jj-npoints2+int(jr))=fxd(ii-npoints2+(ir),jj-npoints2+(jr))+acumx
			fyd(ii-npoints2+int(ir),jj-npoints2+int(jr))=fyd(ii-npoints2+(ir),jj-npoints2+(jr))+acumy
		enddo
	enddo
enddo
! write(*,*)'Fxl',sum(Fxl(1:nl-1))
! write(*,*)'fx',sum(real(fx))
! write(*,*)'Fyl',sum(Fyl(1:nl-1))
! write(*,*)'fy',sum(real(fy))

! call ZFFT2D(fx,NX,NY,-1)
! call ZFFT2D(fy,NX,NY,-1)
! !fx=sigma*fx  !Filtragem
! !fy=sigma*fy
end subroutine distribuicao

function Dij(xk,yk,xii,yjj)
use constantes
implicit none
real(pr)::xk,yk,rx,ry,xx,yy,Dij
integer(pint)::xii,yjj
	rx=abs((xk-x(xii))/dx)
	ry=abs((yk-y(yjj))/dy)
	Dij = (f(rx)*f(ry))/(dx*dy)
return
end function Dij

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! 	FUNCAO PESO  - INTERPOLACAO E DISTRIBUICAO	 !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function f(r)
use constantes
real(pr)::r,f
!!Wang et al. e Lima e Silva
! 	if (r<1.0d0)	then
! 		f=(3.0d0-2.0d0*r+sqrt(1.0d0+4.0d0*r-4.0d0*r*r))*0.125d0
! 	elseif ((r>=1.0d0).and.(r<2.0d0)) then
! 		f=(5.0d0-2.0d0*r-sqrt(-7.0d0+12.0d0*r-4.0d0*r*r))*0.125d0
! 	else
! 		f=0.0d0
! 	end if
! 
! !!Su et al. funcao chapeu
! 	if (r<=1.0_pr)	then
! 		f=(1.0_pr-r)
! 	else
! 		f=0.0_pr
! 	end if

!Peskin et al. (Griffith e Peskin, 2005)
	if (r<1.0_pr)	then
		f=1.0_pr-0.5_pr*r-(r*r)+0.5_pr*r*r*r
	elseif ((r>=1.0_pr).and.(r<2.0_pr)) then
		f=1.0_pr-11.0_pr/6.0_pr*r+(r*r)-1.0_pr/6.0_pr*r*r*r
	else
		f=0.0_pr
	end if
end function f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! 		FUNCAO PARA GERAR RUIDO			 !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function rand2 (idum)
!***********************************************************************************************!
!"long" period ( > 2.e+18 ) random number generator of l'ecuyer with bays-durham shuffle and add
! safeguards. returns a uniform random deviate between 0. and 1 (exclusive of the endpoint values).
! call with idum, a make negative integer to initialize  thereafter, do not alter idum between succesive
! deviates  in a sequence.  rnmax should approximate the largest
! floating value that  is less than 1.  ref: numerical recepies -  pp 271
!***********************************************************************************************!
	integer:: idum,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
	real:: rand2,am,eps,rnmx
	parameter (im1=2147483563 , im2=2147483399 , am=1./im1)
	parameter (imm1=im1-1 , ia1=40014 , ia2=40692)
	parameter (iq1=53668 , iq2=52774 , ir1=12211 , ir2=3791)
	parameter (ntab=64,ndiv=1+imm1/ntab,eps=1.2e-7,rnmx=1.-eps)
	integer:: idum2,j,k,iv(ntab),iy
	save iv,iy,idum2
	data idum2 /123456789/, iv /ntab*0/ , iy /0/
!***********************************************************************************************!
	if (idum.le.0) then    !initialize.
		idum  = max(-idum,1)
		idum2 = idum
		do j=ntab+8,1,-1
			k    = idum/iq1
			idum = ia1*(idum - k*iq1)  - k*ir1
			if (idum.lt.0)   idum = idum + im1
			if (j.le.ntab)   iv(j) = idum
		enddo
	endif
	k    = idum/iq1
	idum = ia1*(idum - k*iq1) - k*ir1
	if (idum.lt.0)  idum = idum + im1
	k     = idum2/iq2
	idum2 = ia2*(idum2 - k*iq2) - k*ir2
	if (idum2.lt.0) idum2 = idum2 + im2
	j     = 1 + iy/ndiv
	iy    = iv(j) - idum2
	iv(j) = idum
	if (iy.lt.1) iy = iy + imm1
!===========================
	rand2 = min(am*iy , rnmx)
!===========================
return
end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!	 		ROTACAO DO CORPO RIGIDO	 	!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine torque
implicit none
integer(pint)::q
real(pr)::F_aux(3),iner(3,3)
      do q=1,front(ifi)%nl-1
		F_aux(1)=front(ifi)%Fxl(q)
		F_aux(2)=front(ifi)%Fyl(q)
		F_aux(3)=0.0_pr
		iner(1,:)=(/     0.0_pr,           0.0_pr,         front(ifi)%rr(2,q)/)
		iner(2,:)=(/     0.0_pr,           0.0_pr,        -front(ifi)%rr(1,q)/)
		iner(3,:)=(/-front(ifi)%rr(2,q),   front(ifi)%rr(1,q),        0.0_pr/)

		front(ifi)%TT=matmul(iner,F_aux)
		front(ifi)%STT=front(ifi)%STT+front(ifi)%TT(3)
	enddo

end subroutine torque

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! 	   	CALCULO DA REGIAO DE FORCAGEM	     	 !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cc_periodico_forca
implicit none
complex(pr),dimension(nx,ny)::fxd,fyd,auxfxd,auxfyd
        auxfxd(1,:)=fxd(1,:)
	auxfyd(1,:)=fyd(1,:)
	fxd(1,:)=fxd(1,:)+fxd(nx,:)
	fxd(nx,:)=fxd(nx,:)+auxfxd(1,:)
	fyd(1,:)=fyd(1,:)+fyd(nx,:)
	fyd(nx,:)=fyd(nx,:)+auxfyd(1,:)
end subroutine cc_periodico_forca

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!	    	GRAU DE LIBERDADE	    	 	 !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine gld_1
implicit none
integer(pint)::q
character(9)::n
real*8::theta,gama
do ifi=1,nfi    
     theta=(2.0_pr*pi)/(front(ifi)%nl-1)
     gama=theta
	call RK46fsi_y
	    do q=1,front(ifi)%nl-1 							
		    front(ifi)%yl(q)=(front(ifi)%ycent+front(ifi)%raio*sin(gama))
		    gama=gama+theta		
	    enddo	
	    front(ifi)%yl(front(ifi)%nl)=front(ifi)%yl(1)
	    front(ifi)%vlf=front(ifi)%vll
enddo	
if ((t>=(tmin-tsave)-0.05_pr*dt).and.(t<=(tmin-tsave)+0.05_pr*dt)) then
    write(485,*) t,(front(1)%ycent-(Ly*0.5d0))/(2.0d0*front(1)%raio)	!posicao do centro eixo y
    write(486,*) t,front(1)%vll	      		    			   ! velocidade centro para o eixo y
    write(487,*) front(1)%vll,front(1)%SFyl*front(1)%ds*dy*rho/dt	 ! forca de amortecimento
endif
pause
end subroutine gld_1


subroutine gld_2
implicit none
integer(pint)::q
character(9)::n
real*8::theta,gama
do ifi=1,nfi    
     theta=(2.0_pr*pi)/(front(ifi)%nl-1)
     gama=theta
	call RK46fsi_x
	 do q=1,front(ifi)%nl-1 							
	      front(ifi)%xl(q)=front(ifi)%xcent+front(ifi)%raio*cos(gama)
	      gama=gama+theta		
	  enddo	
	    front(ifi)%xl(front(ifi)%nl)=front(ifi)%xl(1)
	    front(ifi)%ulf=front(ifi)%ull
enddo	
if ((t>=(tmin-tsave)-0.05_pr*dt).and.(t<=(tmin-tsave)+0.05_pr*dt)) then
  write(490,*) t,(front(1)%xcent-front(1)%xcento)/(2.0d0*front(1)%raio)	!posicao do centro eixo x
  write(491,*) t,front(1)%ull	      		        ! velocidade centro para o eixo x
  write(492,*) (front(1)%xcent-front(1)%xcento)/(2.0d0*front(1)%raio),(front(1)%ycent-(Ly*0.5d0))/(2.0d0*front(1)%raio)	 ! x versus y
  write(493,*) t,sqrt(((front(1)%ycent-(Ly*0.5d0))/(2.0d0*front(1)%raio))**2 +((front(1)%xcent-(front(1)%xcent-(Lx*0.5d0))/(2.0d0*front(1)%raio))/(2.0d0*front(1)%raio))**2)	! posicao final
endif
end subroutine gld_2

subroutine gld_3
implicit none
integer(pint)::q
character(9)::n
real*8::theta,gama
do ifi=1,nfi    
     theta=(2.0_pr*pi)/(front(ifi)%nl-1)
     gama=theta
	call RK46fsi_xy
	    do q=1,front(ifi)%nl-1 							
		    front(ifi)%yl(q)=(front(ifi)%ycent+front(ifi)%raio*sin(gama))
		    front(ifi)%xl(q)=front(ifi)%xcent+front(ifi)%raio*cos(gama)
		    gama=gama+theta		
	    enddo	
	    front(ifi)%yl(front(ifi)%nl)=front(ifi)%yl(1)
	    front(ifi)%vlf=front(ifi)%vll
	    front(ifi)%xl(front(ifi)%nl)=front(ifi)%xl(1)
	    front(ifi)%ulf=front(ifi)%ull
enddo	
if ((t>=(tmin-tsave)-0.05_pr*dt).and.(t<=(tmin-tsave)+0.05_pr*dt)) then
  write(485,*) t,(front(1)%ycent-(Ly*0.5d0))/(2.0d0*front(1)%raio)	!posicao do centro eixo y
  write(486,*) t,front(1)%vll	      		        ! velocidade centro para o eixo y
  write(490,*) t,(front(1)%xcent-front(1)%xcento)/(2.0d0*front(1)%raio)	!posicao do centro eixo x
  write(491,*) t,front(1)%ull	      		        ! velocidade centro para o eixo x
  write(492,*) (front(1)%xcent-front(1)%xcento)/(2.0d0*front(1)%raio),(front(1)%ycent-(Ly*0.5d0))/(2.0d0*front(1)%raio)	 ! x versus y
  write(493,*) t,sqrt(((front(1)%ycent-(Ly*0.5d0))/(2.0d0*front(1)%raio))**2 +((front(1)%xcent-(front(1)%xcent-(Lx*0.5d0))/(2.0d0*front(1)%raio))/(2.0d0*front(1)%raio))**2)	! posicao final
endif
end subroutine gld_3


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!	EULER PARA SOLUCAO DA ESTRUTURA		 !!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine eulerfsi_y
      real*8::delta_y,Clest

      !********************************************************************!
      !*****************  ADMENSIONALIZACAO DAS VARIAVEIS *****************!
      !********************************************************************!
      dt=dt*(umax/(front(1)%raio*2.0d0))
      Clest=-(2.0_pr*front(1)%SFyl*rho/dt*front(1)%ds*dy)/(rho*umax**2*2.0_pr*front(1)%raio)
      front(1)%vll=front(1)%vll/umax
      delta_y=(front(1)%ycent-front(1)%ycento)/(front(1)%raio*2.0d0)
      front(1)%vll= ((2.0_pr*Clest)/(pi*massa)-(4.0_pr*pi*cy)/Ur*front(1)%vll-(2.0_pr*pi/Ur)**2*delta_y)*dt+front(1)%vll
      
      !********************************************************************!
      !*****************  DMENSIONALIZACAO DAS VARIAVEIS *****************!
      !********************************************************************!
      front(1)%vll=front(1)%vll*umax
	dt=dt*((front(1)%raio*2.0d0)/umax) 
      front(1)%ycent=front(1)%ycent+front(1)%vll*dt
end subroutine eulerfsi_y

subroutine eulerfsi_x
      real*8::Cdest,delta_x
      !********************************************************************!
      !*****************  ADMENSIONALIZACAO DAS VARIAVEIS *****************!
      !********************************************************************!
      dt=dt*(umax/(front(1)%raio*2.0d0))
      Cdest=-(2.0_pr*front(1)%SFxl*rho/dt*front(1)%ds*dy)/(rho*umax**2*2.0_pr*front(1)%raio)
      front(1)%ull=front(1)%ull/umax
      delta_x=(front(1)%xcent-front(1)%xcento)/(front(1)%raio*2.0d0)
      front(1)%ull= ((2.0_pr*Cdest)/(pi*massa)-(4.0_pr*pi*cx)/Ur*front(1)%ull-(2.0_pr*pi/Ur)**2*delta_x)*dt+front(1)%ull

      !********************************************************************!
      !*****************  DMENSIONALIZACAO DAS VARIAVEIS *****************!
      !********************************************************************!
      front(1)%ull=front(1)%ull*umax
      dt=dt*((front(1)%raio*2.0d0)/umax) 
      front(1)%xcent=front(1)%xcent+front(1)%ull*dt
end subroutine eulerfsi_x

subroutine eulerfsi_xy
	real*8::delta_y,Clest,Cdest,delta_x
	!********************************************************************!
	!*****************  ADMENSIONALIZACAO DAS VARIAVEIS *****************!
	!********************************************************************!
	dt=dt*(umax/(front(1)%raio*2.0d0))
	
	Cdest=-(2.0_pr*front(1)%SFxl*rho/dt*front(1)%ds*dy)/(rho*umax**2*2.0_pr*front(1)%raio)
	front(1)%ull=front(1)%ull/umax
	delta_x=(front(1)%xcent-front(1)%xcento)/(front(1)%raio*2.0d0)
	front(1)%ull= ((2.0_pr*Cdest)/(pi*massa)-(4.0_pr*pi*cx)/Ur*front(1)%ull-(2.0_pr*pi/Ur)**2*delta_x)*dt+front(1)%ull

	Clest=-(2.0_pr*front(1)%SFyl*rho/dt*front(1)%ds*dy)/(rho*umax**2*2.0_pr*front(1)%raio)
	front(1)%vll=front(1)%vll/umax
	delta_y=(front(1)%ycent-front(1)%ycento)/(front(1)%raio*2.0d0)
	front(1)%vll= ((2.0_pr*Clest)/(pi*massa)-(4.0_pr*pi*cy)/Ur*front(1)%vll-(2.0_pr*pi/Ur)**2*delta_y)*dt+front(1)%vll
	
	!********************************************************************!
	!*****************  DMENSIONALIZACAO DAS VARIAVEIS *****************!
	!********************************************************************!
	front(1)%vll=front(1)%vll*umax
	front(1)%ull=front(1)%ull*umax
	dt=dt*((front(1)%raio*2.0d0)/umax) 
	front(1)%ycent=front(1)%ycent+front(1)%vll*dt
	front(1)%xcent=front(1)%xcent+front(1)%ull*dt
end subroutine eulerfsi_xy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!	RUNGE KUTTA 6 PASSOS OTMIZADO 	Berland et al., 2007		 !!!!!!!
!!!!!!!!!!!!!!!!!!		INTERACAO FLUIDO-ESTSRUTURA			!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine RK46fsi_x
implicit none
	real(pr)::k1x,k1xp,Cdest
	real(pr)::RHSxfip,RHSxfiv
	integer(pint)::q
	character(9)::n
	real*8::theta,gama,delta_x
	k1x=0.0_pr
	k1xp=0.0_pr
	RHSxfip=0.0d0
	RHSxfiv=0.0d0
	!********************************************************************!
	!*****************  ADMENSIONALIZACAO DAS VARIAVEIS *****************!
	!********************************************************************!
	Cdest=-(2.0_pr*front(1)%SFxl*rho/dt*front(1)%ds*dx)/(rho*umax**2*2.0_pr*front(1)%raio)
	dt=dt*(umax/(front(1)%raio*2.0d0))
	front(1)%ull=front(1)%ull/umax
	front(1)%xcent=front(1)%xcent/(front(1)%raio*2.0d0)
	front(1)%xcento=front(1)%xcento/(front(1)%raio*2.0d0)
	    do ii=1,6 
		  delta_x=(front(1)%xcent-front(1)%xcento)
		  RHSxfiv=(2.0_pr*Cdest)/(pi*massa)-(4.0_pr*pi*cx)/Ur*front(1)%ull-(2.0_pr*pi/Ur)**2*delta_x
		  k1x=alfa(ii)*k1x+dt*RHSxfiv
		  front(1)%ull=front(1)%ull+k1x*beta(ii)
		  RHSxfip=front(1)%ull
		  k1xp=alfa(ii)*k1xp+dt*RHSxfip
		  front(1)%xcent=front(1)%xcent+k1xp*beta(ii)
	    enddo
	!********************************************************************!
	!*****************  DMENSIONALIZACAO DAS VARIAVEIS *****************!
	!********************************************************************!
	front(1)%ull=front(1)%ull*umax
	dt=dt*((front(1)%raio*2.0d0)/umax) 
	front(1)%xcent=front(1)%xcent*(front(1)%raio*2.0d0)
	front(1)%xcento=front(1)%xcento*(front(1)%raio*2.0d0)

	if ((t>=(tmin-tsave)-0.05_pr*dt).and.(t<=(tmin-tsave)+0.05_pr*dt)) then
			write(22,*) t,Cdest
	endif
end subroutine RK46fsi_x

subroutine RK46fsi_y
implicit none
      real(pr)::k1y,k1yp,Clest
      real(pr)::RHSyfip,RHSyfiv
      integer(pint)::q
      character(9)::n
      real*8::theta,gama,delta_y
      k1y=0.0_pr
      k1yp=0.0_pr
      RHSyfip=0.0d0
      RHSyfiv=0.0d0
      !********************************************************************!
      !*****************  ADMENSIONALIZACAO DAS VARIAVEIS *****************!
      !********************************************************************!
      Clest=-(2.0_pr*front(1)%SFyl*rho/dt*front(1)%ds*dy)/(rho*umax**2*2.0_pr*front(1)%raio)
      dt=dt*(umax/(front(1)%raio*2.0d0))
      front(1)%vll=front(1)%vll/umax
      front(1)%ycent=front(1)%ycent/(front(1)%raio*2.0d0)
      front(1)%ycento=front(1)%ycento/(front(1)%raio*2.0d0)
	
	  do ii=1,6 
		delta_y=(front(1)%ycent-front(1)%ycento)
		RHSyfiv=(2.0_pr*Clest)/(pi*massa)-(4.0_pr*pi*cy)/Ur*front(1)%vll-(2.0_pr*pi/Ur)**2*delta_y
		k1y=alfa(ii)*k1y+dt*RHSyfiv
		front(1)%vll=front(1)%vll+k1y*beta(ii)
		RHSyfip=front(1)%vll
		k1yp=alfa(ii)*k1yp+dt*RHSyfip
		front(1)%ycent=front(1)%ycent+k1yp*beta(ii)
	  enddo
	  
      !********************************************************************!
      !*****************  DMENSIONALIZACAO DAS VARIAVEIS *****************!
      !********************************************************************!
      front(1)%vll=front(1)%vll*umax
      dt=dt*((front(1)%raio*2.0d0)/umax) 
      front(1)%ycent=front(1)%ycent*(front(1)%raio*2.0d0)
      front(1)%ycento=front(1)%ycento*(front(1)%raio*2.0d0)
	  if ((t>=(tmin-tsave)-0.05_pr*dt).and.(t<=(tmin-tsave)+0.05_pr*dt)) then
		  write(23,*) t,Clest
	  endif
end subroutine RK46fsi_y

subroutine RK46fsi_xy
implicit none
      real(pr)::k1y,k1yp,Clest,k1x,k1xp,Cdest
      real(pr)::RHSyfip,RHSyfiv,RHSxfip,RHSxfiv
      integer(pint)::q
      character(9)::n
      real*8::theta,gama,delta_y,delta_x
      k1y=0.0_pr
      k1yp=0.0_pr
      RHSyfip=0.0d0
      RHSyfiv=0.0d0
      k1x=0.0_pr
      k1xp=0.0_pr
      RHSxfip=0.0d0
      RHSxfiv=0.0d0
      !********************************************************************!
      !*****************  ADMENSIONALIZACAO DAS VARIAVEIS *****************!
      !********************************************************************!

      Clest=-(2.0_pr*front(1)%SFyl*rho/dt*front(1)%ds*dy)/(rho*umax**2*2.0_pr*front(1)%raio)
      Cdest=-(2.0_pr*front(1)%SFxl*rho/dt*front(1)%ds*dx)/(rho*umax**2*2.0_pr*front(1)%raio)
      dt=dt*(umax/(front(1)%raio*2.0d0))
      
      front(1)%vll=front(1)%vll/umax
      front(1)%ull=front(1)%ull/umax
      
      front(1)%ycent=front(1)%ycent/(front(1)%raio*2.0d0)
      front(1)%ycento=front(1)%ycento/(front(1)%raio*2.0d0)
      
      front(1)%xcent=front(1)%xcent/(front(1)%raio*2.0d0)
      front(1)%xcento=front(1)%xcento/(front(1)%raio*2.0d0)
	
	  do ii=1,6 
		delta_y=(front(1)%ycent-front(1)%ycento)
		RHSyfiv=(2.0_pr*Clest)/(pi*massa)-(4.0_pr*pi*cy)/Ur*front(1)%vll-(2.0_pr*pi/Ur)**2*delta_y
		k1y=alfa(ii)*k1y+dt*RHSyfiv
		front(1)%vll=front(1)%vll+k1y*beta(ii)
		RHSyfip=front(1)%vll
		k1yp=alfa(ii)*k1yp+dt*RHSyfip
		front(1)%ycent=front(1)%ycent+k1yp*beta(ii)
		
		delta_x=(front(1)%xcent-front(1)%xcento)
		RHSxfiv=(2.0_pr*Cdest)/(pi*massa)-(4.0_pr*pi*cx)/Ur*front(1)%ull-(2.0_pr*pi/Ur)**2*delta_x
		k1x=alfa(ii)*k1x+dt*RHSxfiv
		front(1)%ull=front(1)%ull+k1x*beta(ii)
		RHSxfip=front(1)%ull
		k1xp=alfa(ii)*k1xp+dt*RHSxfip
		front(1)%xcent=front(1)%xcent+k1xp*beta(ii)
	  enddo 
      !********************************************************************!
      !*****************  DMENSIONALIZACAO DAS VARIAVEIS *****************!
      !********************************************************************!
      front(1)%vll=front(1)%vll*umax
      front(1)%ull=front(1)%ull*umax
      
      dt=dt*((front(1)%raio*2.0d0)/umax) 
      
      front(1)%xcent=front(1)%xcent*(front(1)%raio*2.0d0)
      front(1)%xcento=front(1)%xcento*(front(1)%raio*2.0d0)
      
      front(1)%ycent=front(1)%ycent*(front(1)%raio*2.0d0)
      front(1)%ycento=front(1)%ycento*(front(1)%raio*2.0d0)
      if ((t>=(tmin-tsave)-0.05_pr*dt).and.(t<=(tmin-tsave)+0.05_pr*dt)) then
		      write(22,*) t,Cdest
		      write(23,*) t,Clest
      endif
end subroutine RK46fsi_xy

end module directforcing