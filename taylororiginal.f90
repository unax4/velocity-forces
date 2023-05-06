program taylor0
use tipoak
use funtzioak
use m_interpoli_dp

real(kind=dp)::t,ta,tb,vint,xint,yint,fint,h,xres,errorea
real(kind=dp),dimension(2)::posi,velo,velint,for,forint
integer::i,j,k
integer::n
real(kind=dp),allocatable,dimension(:)::xn,yn
n=1000
j=1

!proiektila
allocate(yn(n/2))
allocate(xn(n/2))
ta=0.0_dp
tb=10.0_dp
t=ta
h=(tb-ta)/n

posi=(/0.0_dp,0.0_dp/)
velo=(/20.0_dp,20.0_dp/)

open(unit=17,file="proiekt0.dat",status="replace",action="write")

write(unit=17,fmt="(3f20.13)")t,posi
write(unit=17,fmt=*)
write(unit=17,fmt=*)
call fg(velo(1),velo(2),for(1),for(2))
posi=posi + h*velo + for*h**2/2
velo=velo + for*h

t=t+h
if (posi(2)<0) then

exit
end if
end do



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