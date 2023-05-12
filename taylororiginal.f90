program taylor0
use tipoak
use funtzioak
use mcf_interpoli
real(kind=dp)::t,ta,tb,h,xres,errorea,b,c
real(kind=dp),dimension(2)::posi,velo,for,xni,yni
real(kind=dp),dimension(100)::xn,yn
integer::i,j,k,p
integer::n
real(kind=dp),parameter::pi=acos(-1.0_dp),mu=81.45_dp/82.45_dp,mup=1.0_dp/82.45_dp !C-ren kalkulurako
n=1000

!proiektila
b=0.0_dp
ta=0.0_dp
tb=10.0_dp
t=ta
h=0.1
open(unit=17,file="proiekt0.dat",status="replace",action="write")

do p=1,5
k=0
j=0
posi=(/0.0_dp,0.0_dp/)
velo=(/20.0_dp,20.0_dp/)

do i=1,100
j=j+1
xn(i)=posi(1)
yn(i)=posi(2)

write(unit=17,fmt="(3f20.13)")t,posi
write(unit=17,fmt=*)
write(unit=17,fmt=*)
call fg(velo(1),velo(2),for(1),for(2),b)
posi=posi + h*velo + for*h**2/2
velo=velo + for*h

t=t+h
if (k==1) then
exit
end if
if (posi(2)<0) then
k=1
end if

end do
xni=(/xn(j-1),xn(j)/)
yni=(/yn(j-1),yn(j)/)
call POLINT(yni,xni,2,0.0_dp,xres,errorea)
print*,xres
b=b+0.5_dp
end do
!magnetikoa
n=1000
t=0.0_dp

posi=(/1.0_dp,0.0_dp/)
velo=(/0.0_dp,2*acos(-1.0_dp)/)
h=0.01_dp
open(unit=18,file="mag0.dat",status="replace",action="write")
do i=1,601
write(unit=18,fmt="(3f20.13)")t,posi

call fmag(velo(1),velo(2),for(1),for(2))
posi=posi + h*velo + for*h**2/2
call fmag(velo(1),velo(2),for(1),for(2))
velo=velo + for*h
t=t+h

end do

!harmonikoa
t=0.0_dp
tb=10.0_dp
n=5000
h=(tb-t)/5000.0_dp
posi=(/1.0_dp,0.0_dp/)
velo=(/200.0_dp*pi,0.0_dp/)

open(unit=19,file="harm0.dat",status="replace",action="write")
do i=1,n+1
write(unit=19,fmt="(3f20.13)")t,posi
call fharm(t,posi(1),velo(1),for(1))
posi=posi + h*velo + for*h**2/2
call fharm(t,posi(1),velo(1),for(1))
velo=velo + for*h
t=t+h
end do

!grabitatorioa
t=0.0_dp

tb=6.192169331396_dp
n=1000000
h=(tb-t)/1000000.0_dp

posi=(/1.2_dp,0.0_dp/)
velo=(/0.0_dp,-1.049357509830_dp/)

open(unit=20,file="orb0.dat",status="replace",action="write")

do i=1,1000000
c=(velo(1)**2.0_dp+velo(2)**2.0_dp-posi(1)**2.0_dp-posi(2)**2.0_dp)/2&
-mu/((posi(1)+mup)**2.0_dp+posi(2)**2.0_dp)**(1.0_dp/2.0_dp)&
-mup/((posi(1)-mu)**2.0_dp+posi(2)**2.0_dp)**(1.0_dp/2.0_dp)
write(unit=20,fmt="(4f20.13)")t,posi,c

call forbi(posi(1),posi(2),velo(1),velo(2),for(1),for(2))
posi=posi + h*velo + for*h**2/2
call forbi(posi(1),posi(2),velo(1),velo(2),for(1),for(2))
velo=velo + for*h
t=t+h
end do

end program taylor0
