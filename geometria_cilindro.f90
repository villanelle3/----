module geometria
use constantes
real(pr)::epsp,epspl,epsw,epswl,rhoff,Fxcp,Fycp,xli,yli,dil
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 	GERA GEOMETRIA 		 !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine gera_geometria
use constantes
implicit none

real*8::theta,gama,pos
real(pr),dimension(nx,ny)::ul,vl,wwp
integer::q
!******************************************************************************************!
!****************** Posicao das condicoes de contorno INDICES DOS DOMINIO FISICO  *********!
!******************************************************************************************!
iimin=int(1.8d0*nx/Lx)
iimax=nx+1
jjmin=1
jjmax=ny
!******************************************************************************************!
!****************** Posicao da fronteira imersa COINCIDENTES  *****************************!
!******************************************************************************************!
imin=iimin+64
imax=imin+8
jmin=ny/2
jmax=int(ny/2)

!******************************************************************************************!
!****************** Posicao das sondas CASO ESCOAMENTO SOBRE CILINDRO  ********************!
!******************************************************************************************!

px1=iimin-zp
px2=iimin
px3=int(front(1)%xcent/dx)+int(0.5d0/dx)
px4=int(front(1)%xcent/dx)+int(1.0d0/dx)
px5=int(front(1)%xcent/dx)+int(1.5d0/dx)
px6=int(front(1)%xcent/dx)-int(0.5d0/dx)
px7=int(front(1)%xcent/dx)-int(1.0d0/dx)
px8=int(front(1)%xcent/dx)-int(1.5d0/dx)

py1=Ny*0.5d0
py2=Ny*0.5d0-int(1.5d0/dy)
py3=Ny*0.5d0+int(1.5d0/dy)

!******************************************************************************************!
!******************	Gera RETANGULOS  EX: Couette circular	*******************!
!******************************************************************************************!

front(1)%altura=1.0d0

if (start==0) then
    front(1)%xcent=64*dx+2.0d0
    front(2)%xcent=64*dx
    front(1)%ycent=2.0d0 !+front(1)%altura*0.5d0	!ycent=posicao da reta na coordenda y
    front(2)%ycent=2.0d0	!ycent=posicao da reta na coordenda y
      else
	front(1)%xcent=x(64)+1.0d-2		!posicao x do centro da esfera
	open(31,file='restart_2.dat',status='old')
		      read(31,*)t
		      read(31,*)nstep
		      read(31,*)nn
		      read(31,*)tmin
		      read(31,*)front(1)%yin
		      read(31,*)front(1)%vll
		    do j=1,ny
			    do i=1,nx
				    read(31,*) x(i),y(j),up(i,j),vp(i,j),pp(i,j),wwp(i,j)
			    enddo
		    enddo
	 close(31)
endif

!******************************************************************************************!
!***************	QUEDA DE PARTICULA Parametros das forcas de repulsao	***********!
!******************************************************************************************!

 
front(1)%raio=1.0d0		
front(2)%raio=6.0d0

do ifi=1,nfi
	front(ifi)%nl=int(front(ifi)%raio/(min(dx,dy)))

	allocate(front(ifi)%xl(front(ifi)%nl),front(ifi)%yl(front(ifi)%nl),front(ifi)%ulf(front(ifi)%nl),front(ifi)%vlf(front(ifi)%nl),front(ifi)%ul(front(ifi)%nl),front(ifi)%vl(front(ifi)%nl),front(ifi)%Fxl(front(ifi)%nl),front(ifi)%Fyl(front(ifi)%nl),front(ifi)%ulp(front(ifi)%nl),front(ifi)%vlp(front(ifi)%nl),front(ifi)%rro(3,front(ifi)%nl),front(ifi)%omega_aux(3,front(ifi)%nl),front(ifi)%iner(3,front(ifi)%nl),front(ifi)%rr(3,front(ifi)%nl), stat=status)
	if (status/=0) then
		write(*,*)'A alocacao das variaveis lagrangianas falhou!'
		stop
	end if
	!inicializacao das variaveis
	front(ifi)%ulf=0.0_pr
	front(ifi)%vlf=0.0_pr
	front(ifi)%ul=0.0_pr
	front(ifi)%vl=0.0_pr
	front(ifi)%xl=0.0_pr
	front(ifi)%yl=0.0_pr
enddo
do ifi=1,nfi
	front(ifi)%xl=0.0_pr
	front(ifi)%yl=0.0_pr

	do q=1,front(ifi)%nl-1
		if (ifi==1) then
			front(ifi)%yl(q)=front(ifi)%ycent+(q-1)*real(front(ifi)%raio/front(ifi)%nl)
			front(ifi)%xl(q)=front(ifi)%xcent
		else
			front(ifi)%xl(q)=front(ifi)%xcent+q*real(front(ifi)%raio/front(ifi)%nl)
			front(ifi)%yl(q)=front(ifi)%ycent
		endif
	enddo

	front(ifi)%xl(front(ifi)%nl)=front(ifi)%xl(1)
	front(ifi)%yl(front(ifi)%nl)=front(ifi)%yl(1)
	front(ifi)%ds=front(ifi)%raio/front(ifi)%nl
enddo
front(1)%omega=0.0_pr
front(2)%omega=0.0_pr


raio1=front(1)%raio
raio2=front(2)%raio

omega1=front(1)%omega
omega2=front(2)%omega


!******************************************************************************************!
!******************** 			GERACAO DA TORRE SOBRE TERRENO	*******************!
!******************************************************************************************!	

open(47,file='GEOMETRIA.dat',status='unknown')!,form="unformatted")
	do ifi=1,nfi
		do q=1,front(ifi)%nl
			write(47,*) front(ifi)%xl(q),front(ifi)%yl(q)
		enddo
	enddo
close(47)
end subroutine gera_geometria


!******************************************************************************************!
!******************** 			GERACAO DA VIGA		*******************!
!******************************************************************************************!	

! if (start==0) then
! 	front(1)%ycent=Ly*0.5d0			!posicao y do centro da esfera
! 	front(1)%xcent=x(iimin)+15.0d0*(2.0d0*front(1)%raio)		!posicao x do centro da esfera
!       else
! 	front(1)%xcent=x(iimin)+15.0d0*(2.0d0*front(1)%raio)		!posicao x do centro da esfera
! 	open(31,file='restart_2.dat',status='replace')
! 		  read(31,*)front(1)%ycent
! 		  read(31,*)front(1)%vll
! 	      close(31)
! endif

! do ifi=1,nfi	
! 	do q=1,front(ifi)%nl-1 
! 		if (q*(min(dx,dy)) <=front(ifi)%altura) then
! 			front(ifi)%yl(q)=front(ifi)%ycent+q*(min(dx,dy))
! 			front(ifi)%xl(q)=front(ifi)%xcent
! 		else if (q*(min(dx,dy)) >front(ifi)%altura .and. q*(min(dx,dy)) <=front(ifi)%altura+front(ifi)%largura) then
! 			front(ifi)%xl(q)=front(ifi)%xcent-front(ifi)%yl(q)+q*(min(dx,dy))
! 			front(ifi)%yl(q)=front(ifi)%ycent+front(ifi)%altura
! 		else if (q*(min(dx,dy)) >front(ifi)%altura+front(ifi)%largura .and. q*(min(dx,dy)) <=2.0d0*front(ifi)%altura+front(ifi)%largura) then
! 			front(ifi)%xl(q)=front(ifi)%xcent+front(ifi)%largura
! 			front(ifi)%yl(q)=front(ifi)%ycent+q*(min(dx,dy))-(front(ifi)%altura+front(ifi)%largura)
! 		else
! 			front(ifi)%xl(q)=front(ifi)%xcent+q*(min(dx,dy))-(2.0d0*front(ifi)%altura+front(ifi)%largura)
! 			front(ifi)%yl(q)=front(ifi)%ycent
! 		endif
! 	enddo
! 	front(ifi)%xl(front(ifi)%nl)=front(ifi)%xl(1)
! 	front(ifi)%yl(front(ifi)%nl)=front(ifi)%yl(1)
! 	! front(ifi)%ds=front(ifi)%raio*theta

! 	front(ifi)%xcento=front(ifi)%xcent				! posicao inicial do centro da esfera X
! 	front(ifi)%ycento=front(ifi)%ycent 					! posicao inicial do centro da esfera Y
	
! 	front(ifi)%rro(1,:)=front(ifi)%xl-front(ifi)%xcent
! 	front(ifi)%rro(2,:)=front(ifi)%yl-front(ifi)%ycent
! 	front(ifi)%rro(3,:)=0.0_pr
! enddo

! front(1)%omega=0.0_pr
! raio1=front(1)%raio
! omega1=front(1)%omega


end module geometria