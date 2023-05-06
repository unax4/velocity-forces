module funtzioak 
use tipoak
public::fg,fmag,fharm,forb

contains
subroutine fg(vx,vy,fx,fy)
real(kind=dp),intent(in)::vx,vy
real(kind=dp),intent(inout)::fx,fy
real(kind=dp),parameter::b=0.3_dp,m=2.0_dp,g=9.8_dp

fx=-b*vx
fy=-b*vy - g
end subroutine fg

subroutine fmag(vx,vy,fx,fy)
real(kind=dp),intent(in)::vx,vy
real(kind=dp),intent(out)::fx,fy
real(kind=dp),parameter::B=2*acos(-1.0_dp),q=1.0_dp,m=1.0_dp
fx=B*q*vy/m
fy=-B*q*vx/m
end subroutine fmag


subroutine fharm(t,x,y,vx,vy,fx,fy)
real(kind=dp),intent(in)::vx,vy
real(kind=dp),intent(out)::fx,fy
fx=-(4*pi**2+1)*x-2*vx+sqrt(9*pi**4+10*pi**2+1)*cos(pi*t+atan(2*pi/(3*pi**2+1)))
fy=-(4*pi**2+1)*y-2*vy+sqrt(9*pi**4+10*pi**2+1)*cos(pi*t+atan(2*pi/(3*pi**2+1)))
end subroutine fharm

subroutine forb(x,y,vx,vy,fx,fy)
real(kind=dp),intent(in)::vx,vy
real(kind=dp),intent(out)::fx,fy
real(kind=dp),parameter::pi=acos(-1.0_dp)
fx=0.0
fy=1.0_dp
end subroutine forb

end module funtzioak