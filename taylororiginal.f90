program taylor0
use tipoak
use funtzioak
use mcf_interpoli
real(kind=dp)::t,ta,tb,h,xres,errorea
real(kind=dp),dimension(2)::posi,velo,for,xni,yni
real(kind=dp),dimension(100)::xn,yn
integer::i,j,k
integer::n
n=1000
j=0
k=0
!proiektila
ta=0.0_dp
tb=10.0_dp
t=ta
h=0.1

posi=(/0.0_dp,0.0_dp/)
velo=(/20.0_dp,20.0_dp/)
open(unit=17,file="proiekt0.dat",status="replace",action="write")
do i=1,100
j=j+1
xn(i)=posi(1)
yn(i)=posi(2)

write(unit=17,fmt="(3f20.13)")t,posi
write(unit=17,fmt=*)
write(unit=17,fmt=*)
call fg(velo(1),velo(2),for(1),for(2))
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

!magnetikoa
n=1000
t=0.0

posi=(/1.0_dp,0.0_dp/)
velo=(/0.0_dp,2*acos(-1.0_dp)/)
h=0.01
open(unit=18,file="mag0.dat",status="replace",action="write")
do i=1,1000
write(unit=18,fmt="(3f20.13)")t,posi

call fmag(velo(1),velo(2),for(1),for(2))
posi=posi + h*velo + for*h**2/2
velo=velo + for*h
t=t+h
j=j+1
if (j>100000) then

exit
end if

end do

!harmonikoa
t=0.0
h=(tb-t)/n
posi=(/0.0_dp,0.0_dp/)
velo=(/0.0_dp,0.0_dp/)

open(unit=19,file="harm0.dat",status="replace",action="write")
do i=1,n
write(unit=19,fmt="(3f20.13)")t,posi
call fharm(t,posi(1),posi(2),velo(1),velo(2),for(1),for(2))
posi=posi + h*velo + for*h**2/2
velo=velo + for*h
t=t+h
end do

!grabitatorioa
t=0.0
tb=100
n=10000
h=(tb-t)/n
posi=(/1.2_dp,0.0_dp/)
velo=(/0.0_dp,1.049357509830_dp/)

open(unit=20,file="orb0.dat",status="replace",action="write")

do i=1,n
write(unit=20,fmt="(3f20.13)")t,posi

call forbi(posi(1),posi(2),velo(1),velo(2),for(1),for(2))
posi=posi + h*velo + for*h**2/2
velo=velo + for*h
t=t+h
end do

end program taylor0
