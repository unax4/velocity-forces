program euler
use tipoak
real(kind=dp)::x,v,t,ta,tb,xmid,vmid,tmid,amid,h
integer::i,j,k,n
ta=0.0_dp
tb=3.0_dp
x=0.0_dp
v=10.0_dp
t=ta
n=100


h=(tb-ta)/n

open(unit=12,file="soluzioa.dat",status="replace",action="write")

do n=1,100
write(unit=12,fmt="(3f10.6)")t,x,v

tmid=t+h/2
vmid=v+f(t,x,v)*h/2
xmid=x+v*h/2
amid=f(tmid,xmid,vmid)


t=t+h
v=v+amid*h
x=x+vmid*h

end do

close(unit=12)
contains 
function  f(r,c,s)
real(kind=dp),intent(in)::s,r,c
real(kind=dp)::F
real(kind=dp)::b
b=1_dp
F=-b*s**2*sin(r)
end function F

end program euler
