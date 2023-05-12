program taylor
use tipoak
use funtzioak
use mcf_interpoli
real(kind=dp)::xres,err,t,ta,tb,vint,xint,yint,fint,h,c,b,theta
real(kind=dp),dimension(2)::posi,velo,velint,for,forint
integer::i,j,k,p,q
integer::n
real(kind=dp),dimension(100)::xn,yn
real(kind=dp),dimension(2)::xni,yni
real(kind=dp),parameter::pi=acos(-1.0_dp),mu=81.45_dp/82.45_dp,mup=1.0_dp/82.45_dp !C-ren kalkulurako

print*,sin(0.5)
n=1000
k=0
j=1
!proiektila
b=0.0_dp

open(unit=13,file="proiekt.dat",status="replace",action="write")

do p=1,5
t=0.0_dp
h=0.1
posi=(/0.0_dp,0.0_dp/)
velo=(/20.0_dp,20.0_dp/)

k=0
j=1
do i=1,100
k=k+1
xn(i)=posi(1)
yn(i)=posi(2)

write(unit=13,fmt="(3f20.13)")t,posi
write(unit=13,fmt=*)
write(unit=13,fmt=*)


call fg(velo(1),velo(2),for(1),for(2),b)
posi=posi + h*velo + for*h**2/2

velint=velo + for*h
call fg(velint(1),velint(2),forint(1),forint(2),b)

velo=velo + (for+forint)*h/2
t=t+h
if (j==0) then
exit
end if

if (posi(2)<0) then
j=0
end if

end do
!xni=(/xn(k-1),xn(k)/)
!yni=(/yn(k-1),yn(k)/)


!call polint(yni,xni,2,0.0_dp,xres,err)
!print*,xres
b=b+0.5_dp
end do



!proiektila angeluekin
b=0.1_dp

open(unit=30,file="proiektheta.dat",status="replace",action="write")

do q=1,5
theta=pi/12.0_dp
do p=1,7
t=0.0_dp
h=0.1_dp
posi=(/0.0_dp,0.0_dp/)
velo=(/30.0_dp*cos(theta),30.0_dp*sin(theta)/)

k=0
j=1
do i=1,100
k=k+1
xn(i)=posi(1)
yn(i)=posi(2)

write(unit=30,fmt="(3f20.13)")t,posi
write(unit=30,fmt=*)
write(unit=30,fmt=*)


call fg2(velo(1),velo(2),for(1),for(2),b)
posi=posi + h*velo + for*h**2/2

velint=velo + for*h
call fg2(velint(1),velint(2),forint(1),forint(2),b)

velo=velo + (for+forint)*h/2
t=t+h
if (j==0) then
exit
end if

if (posi(2)<0) then
j=0
end if

end do
xni=(/xn(k-1),xn(k)/)
yni=(/yn(k-1),yn(k)/)


call polint(yni,xni,2,0.0_dp,xres,err)
print*,xres,theta

theta=theta + pi/36.0_dp
end do
b=b+0.1_dp
j=1
k=0
end do

!magnetikoa
t=0.0_dp
j=1
tb=6.0_dp
posi=(/1.0_dp,0.0_dp/)
velo=(/0.0_dp,2*acos(-1.0_dp)/)
h=0.01_dp
open(unit=14,file="mag.dat",status="replace",action="write")
do i=1,6000

write(unit=14,fmt="(3f20.13)")t,posi


call fmag(velo(1),velo(2),for(1),for(2))
posi=posi + h*velo + for*h**2/2

velint=velo + for*h
call fmag(velint(1),velint(2),forint(1),forint(2))

velo=velo + (for+forint)*h/2
t=t+h
end do

!harmonikoa
n=5000_dp

t=0.0_dp
tb=10.0_dp
h=(tb-t)/5000.0_dp
posi=(/1.0_dp,0.0_dp/)
velo=(/200.0_dp*pi,0.0_dp/)

open(unit=15,file="harm.dat",status="replace",action="write")
do i=1,n+1
write(unit=15,fmt="(3f20.13)")t,posi



call fharm(t,posi(1),velo(1),for(1))
posi(1)=posi(1) + h*velo(1) + for(1)*h**2/2

velint(1)=velo(1) + for(1)*h
call fharm(t+h,posi(1),velint(1),forint(1))

velo(1)=velo(1) + (for(1)+forint(1))*h/2
t=t+h

if (t==tb) then

exit
end if

end do



!grabitatorioa
t=0.0_dp
tb=6.192169331396_dp
n=8000.0_dp

h=(tb-t)/n
posi=(/1.2_dp,0.0_dp/)
velo=(/0.0_dp,-1.049357509830_dp/)

open(unit=16,file="orb.dat",status="replace",action="write")
do i=1,8001
c=(velo(1)**2.0_dp+velo(2)**2.0_dp-posi(1)**2.0_dp-posi(2)**2.0_dp)/2&
-mu/((posi(1)+mup)**2.0_dp+posi(2)**2.0_dp)**(1.0_dp/2.0_dp)&
-mup/((posi(1)-mu)**2.0_dp+posi(2)**2.0_dp)**(1.0_dp/2.0_dp)


write(unit=16,fmt="(4f20.13)")t,posi,c
write(unit=16,fmt=*)
write(unit=16,fmt=*)

call forbi(posi(1),posi(2),velo(1),velo(2),for(1),for(2))
posi=posi + h*velo + for*h**2/2

velint=velo + for*h
call forbi(posi(1),posi(2),velint(1),velint(2),forint(1),forint(2))

velo=velo + (for+forint)*h/2
t=t+h


end do


end program taylor
