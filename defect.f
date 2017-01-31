!If you decide to use this code, please cite "Fusion Science and Technology, 66, 235 (2014)" Thanks!
program main  

implicit none

integer, parameter :: num=100   !the timestep of radiativeatom
integer, parameter :: numb=100 !the timestep of referenceatom
integer, parameter :: M=12583   !the number of atoms in one cycle
!integer, parameter :: n=300000
doubleprecision,parameter :: first=1.89
integer :: i,j,k,s,Na,en,nk,nr,stat,r_id,p_id,dn,nvt,n1,n2,n3,nat,n4,n5,n6,nit,n7,n8,n9,fn
integer :: cvacancysum=0 
integer :: sivacancysum=0
integer :: cantisitesum=0, siantisitesum=0, Cinterstitialsum=0,siinterstitialsum=0
type :: reference	
 	integer,allocatable :: reference_Element(:)
	doubleprecision,allocatable :: posx(:)
	doubleprecision,allocatable :: posy(:)
	doubleprecision,allocatable :: posz(:)
end type reference

type :: radiative
	integer,allocatable :: radiative_element(:)
	doubleprecision,allocatable :: Ptempx(:)
	doubleprecision,allocatable :: Ptempy(:)
	doubleprecision,allocatable :: Ptempz(:)
end type radiative

doubleprecision,allocatable :: dist(:)  
       
doubleprecision :: min=0.0

integer :: error1=0,error2=0  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
en=20      !the number of cycle in displacement
dn=10
fn=30

allocate(dist(m))
allocate(radiative_Element(m),stat=error1)

allocate(reference_Element(m),stat=error2)

allocate(Ptempx(m))
allocate(Ptempy(m))
allocate(Ptempz(m))
allocate(Posx(m))
allocate(Posy(m))
allocate(Posz(m))
!allocate(TIMESTEP(en))

nk=en+1
nr=dn+2
!!!!!!!!!!!!!!!!read data from the dump.irradiated
open(unit=en+1,file="dump.final")
open(unit=en+5,file="radiativenatoms")
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do k=1,num            
     read(nk,*,iostat=stat) 
	 if(stat/=0) then
	   write(*,*) k
	   stop
	 end if  
     read(nk,*,iostat=stat)  !TIMESTEP(k)
!     write(*,*) Timestep(k)
	 if(stat/=0) then
	   write(*,*) k
	   stop
	 end if 
     do i=1,7      
      read(nk,*)
     end do  	
    do j=1,m
     read(nk,*,iostat=stat) r_id,radivative%radiative_element(j),radivative%Ptempx(j),radivative%ptempy(j),radivative%ptempz(j)
	if(error1/=0) then
		write(*,*) "wrong!"
		stop
	end if
!	139 format(1x,i8,i2,3f16.6) 
	if(stat/=0) then
		write(*,141) k,j,r_id,radivative%radiative_element(j),radivative%Ptempx(j),radivative%ptempy(j),radivative%ptempz(j)
		stop
	end if    
       if(k==num) then
       		write(en+5,139) r_id,radivative%radiative_element(j),radivative%Ptempx(j),radivative%ptempy(j),radivative%ptempz(j)	
       end if

!136 format(1x,A,3f16.9)
139 format(1x,i8,i2,3f16.6)
141 format(1x,i2,i8,i10,i2,3f16.6)
    end do
 end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!read data from the dump.start
open(unit=dn+6,file="referenceatoms")    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     do j=1,m
       read(dn+6,*)   p_id,reference%reference_element(j),reference%posx(j),reference%posy(j),reference%posz(j)
	if(error2/=0) then
		write(*,*) "wrong!"
		stop
	end if
    end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
nvt=fn+1
n1=fn+2
n2=fn+3
n3=fn+4
nat=fn+5
n4=fn+6
n5=fn+7
n6=fn+8
nit=fn+9
n7=fn+10
n8=fn+11
n9=fn+12
open(unit=fn+1,file="data_v")
open(unit=fn+2,file="data_vx")
open(unit=fn+3,file="data_vy")
open(unit=fn+4,file="data_vz")
open(unit=fn+5,file="data_a")
open(unit=fn+6,file="data_ax")
open(unit=fn+7,file="data_ay")
open(unit=fn+8,file="data_az")
open(unit=fn+9,file="data_i")
open(unit=fn+10,file="data_ix")
open(unit=fn+11,file="data_iy")
open(unit=fn+12,file="data_iz")
!!!!!!choose a reference atom!!!!!!!!!!
do i=1,m
	do j=1,m
		dist(j)=sqrt((posx(i)-ptempx(j))**2.0d0+(posy(i)-ptempy(j))**2.0d0+(posz(i)-ptempz(j))**2.0d0)
	end do
	min=dist(1)
	do j=2,m
		if(min>dist(j)) then
			min=dist(j)
		end if
	end do
	if(min>0.5*first) then
		if(reference_element(i)==1) then
			Cvacancysum=cvacancysum+1
			write(*,"('cvacancy is ',f8.4, f8.4, f8.4)") posx(i),posy(i),posz(i)
			write(nvt,"(1x,i2)") reference_element(i)
			write(n1,"(1x,f8.4)") posx(i)
			write(n2,"(1x,f8.4)") posy(i)
			write(n3,"(1x,f8.4)") posz(i)
		else
			sivacancysum=sivacancysum+1
			write(*,"('sivacancy is ',f8.4, f8.4, f8.4)") posx(i),posy(i),posz(i)
			write(nvt,"(1x,i2)") reference_element(i)
			write(n1,"(1x,f8.4)")posx(i)
			write(n2,"(1x,f8.4)")posy(i)
			write(n3,"(1x,f8.4)")posz(i)
		end if
	end if
end do
min=0.0
write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!"
do i=1,m
	do j=1,m
		dist(j)=sqrt((posx(i)-ptempx(j))**2.0d0+(posy(i)-ptempy(j))**2.0d0+(posz(i)-ptempz(j))**2.0d0)
	end do
	min=dist(1)
	s=reference_element(1)
	do j=2,m
		if(min>dist(j)) then
			min=dist(j)
			s=radiative_element(j)
		end if
	end do
	if(min<=0.5*first) then
		if(reference_element(i)/=s) then
			if(reference_element(i)==1) then 
				cantisitesum=cantisitesum+1
				write(*,"('cantisite is ',f8.4, f8.4, f8.4)") posx(i),posy(i),posz(i)
				write(nat,"(1x,i2)") reference_element(i)
				write(n4,"(1x,f8.4)")posx(i)
				write(n5,"(1x,f8.4)")posy(i)
				write(n6,"(1x,f8.4)")posz(i)
			else
				siantisitesum=siantisitesum+1
				write(*,"('siantisite is ',f8.4, f8.4, f8.4)") posx(i),posy(i),posz(i)
				write(nat,"(1x,i2)") reference_element(i)
				write(n4,"(1x,f8.4)")posx(i)
				write(n5,"(1x,f8.4)")posy(i)
				write(n6,"(1x,f8.4)")posz(i)
			end if
		end if
	end if
end do
!!!!!!!choose an irradiated atom!!!!!!!!
write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!"
do i=1,m
	do j=1,m
		dist(j)=sqrt((posx(j)-ptempx(i))**2.0d0+(posy(j)-ptempy(i))**2.0d0+(posz(j)-ptempz(i))**2.0d0)
	end do
	min=dist(1)
	s=reference_element(1)
	do j=2,m
		if(min>dist(j)) then
			min=dist(j)
		end if
	end do
	if(min>0.5*first .and. min<first) then
		if(radiative_element(i)==1) then
			cinterstitialsum=cinterstitialsum+1
			write(*,"('cinterstitial is ',f8.4, f8.4, f8.4)") ptempx(i),ptempy(i),ptempz(i)
			write(nit,"(1x,i2)") radiative_element(i)
			write(n7,"(1x,f8.4)")ptempx(i)
			write(n8,"(1x,f8.4)")ptempy(i)
			write(n9,"(1x,f8.4)")ptempz(i)
		else
			siinterstitialsum=siinterstitialsum+1
			write(*,"('siinterstitial is ',f8.4, f8.4, f8.4)") ptempx(i),ptempy(i),ptempz(i)
			write(nit,"(1x,i2)") radiative_element(i)
			write(n7,"(1x,f8.4)")ptempx(i)
			write(n8,"(1x,f8.4)")ptempy(i)
			write(n9,"(1x,f8.4)")ptempz(i)
		end if
	end if
end do
write(*,*) "the number of C vacancy is "
write(*,"(I3)") cvacancysum
write(*,*) "the number of Si vacancy is "
write(*,"(I3)") sivacancysum
write(*,*) "the number of C interstitial is "	
write(*,"(I3)") cinterstitialsum
write(*,*) "the number of Si interstitial is "	
write(*,"(I3)") siinterstitialsum
write(*,*) "the number of C antisite is "
write(*,"(I3)") cantisitesum
write(*,*) "the number of Si antisite is "
write(*,"(I3)") siantisitesum	
end 
