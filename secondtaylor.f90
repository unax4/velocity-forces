program taylor
use tipoak
real(kind=dp)::x,y,v,t,f,ta,tb,vint,xint,yint,fint,h
integer::i,j,k
real,dimension(2)::q,s
integer,parameter::n=100


!bektore gabe
ta=0_dp
tb=10_dp
x=0
y=3_dp
v=15_dp
t=ta
h=(tb-ta)/n
f=force(t,x,v)

do i=1,n
write(unit=*,fmt="(3f12.4)")t,x,v

x=x+h*n*v+f*h**2/2
vint=v+f*h
fint=force(t+h,x,vint)

v=v+(f+fint)*h/2
t=t+h
f=force(t,x,v) !aire erresistentzia problementzat azken  jarri, f txikia bada indar kontserbatzaileekin konparatuz kendu eta f=fint)
end do



!forma bektoriala
q(1)=x
q(2)=v
s(1)=v
s(2)=force


contains
function force(temp,pos,vel)
real(kind=dp),intent(in)::temp,pos,vel
real(kind=dp)::force
real(kind=dp),parameter::c=-0.3_dp
force=c*vel
end function force

end program taylor
