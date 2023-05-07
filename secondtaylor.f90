program taylor
use tipoak
use funtzioak
use mcf_interpoli
real(kind=dp)::xres,err,t,ta,tb,vint,xint,yint,fint,h,c
real(kind=dp),dimension(2)::posi,velo,velint,for,forint
integer::i,j,k
integer::n
real(kind=dp),dimension(:),allocatable::xn,yn,xni,yni
real(kind=dp),parameter::mu=81.45_dp/82.45_dp,mup=1.0_dp/82.45_dp
n=1000
allocate(xn(100))
allocate(yn(100))
k=0
j=1
!proiektila
h=0.1
posi=(/0.0_dp,0.0_dp/)
velo=(/20.0_dp,20.0_dp/)

open(unit=13,file="proiekt.dat",status="replace",action="write")

do i=1,100
k=k+1
xn(i)=posi(1)
yn(i)=posi(2)

write(unit=13,fmt="(3f20.13)")t,posi
write(unit=13,fmt=*)
write(unit=13,fmt=*)


call fg(velo(1),velo(2),for(1),for(2))
posi=posi + h*velo + for*h**2/2

velint=velo + for*h
call fg(velint(1),velint(2),forint(1),forint(2))

velo=velo + (for+forint)*h/2
t=t+h
if (j==0) then
exit
end if

if (posi(2)<0) then
j=0
end if

end do
allocate(xni(2))
allocate(yni(2))
xni=(/xn(k-1),xn(k)/)
yni=(/yn(k-1),yn(k)/)


call polint(yni,xni,2,0.0_dp,xres,err)
print*,xres

!magnetikoa
t=0.0
j=1

posi=(/1.0_dp,0.0_dp/)
velo=(/0.0_dp,2*acos(-1.0_dp)/)
h=0.01
open(unit=14,file="mag.dat",status="replace",action="write")
do i=1,1000
write(unit=14,fmt="(3f20.13)")t,posi


call fmag(velo(1),velo(2),for(1),for(2))
posi=posi + h*velo + for*h**2/2

velint=velo + for*h
call fmag(velint(1),velint(2),forint(1),forint(2))

velo=velo + (for+forint)*h/2
t=t+h
j=j+1
if (j>15000) then

exit
end if

end do

!harmonikoa
t=0.0
h=(tb-t)/n
posi=(/0.0_dp,0.0_dp/)
velo=(/0.0_dp,0.0_dp/)

open(unit=15,file="harm.dat",status="replace",action="write")
do i=1,n
write(unit=15,fmt="(3f20.13)")t,posi



call fharm(t,posi(1),posi(2),velo(1),velo(2),for(1),for(2))
posi=posi + h*velo + for*h**2/2

velint=velo + for*h
call fharm(t+h,posi(1),posi(2),velint(1),velint(2),forint(1),forint(2))

velo=velo + (for+forint)*h/2
t=t+h

if (t==tb) then

exit
end if

end do



!grabitatorioa
t=0.0
tb=100
n=10000
h=(tb-t)/n
posi=(/1.2_dp,0.0_dp/)
velo=(/0.0_dp,1.049357509830_dp/)

open(unit=16,file="orb.dat",status="replace",action="write")
do i=1,n
if (t>6.20) then
exit
end if
c=(velo(1)**2+velo(2)**2-posi(1)**2-posi(2)**2)/2-mu/sqrt((posi(1)+mup)**2+posi(2)**2)-mup/sqrt((posi(1)-mu)**2+posi(2)**2)
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
