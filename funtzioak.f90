module funtzioak 
use tipoak
public::fg,fmag,fharm,forb

contains
subroutine fg(vx,vy,fx,fy)
real(kind=dp),intent(in)::vx,vy
real(kind=dp),intent(out)::fx,fy
real(kind=dp),parameter::b=0.3_dp,m=2.0_dp,g=9.8_dp

fx=-b*vx
fy=-b*vy - g
end subroutine force

subroutine fmag(vx,vy,fx,fy)
real(kind=dp),intent(in)::vx,vy
real(kind=dp),intent(out)::fx,fy
real(kind=dp),parameter::B=acos(-1.0_dp)*2,q=-1.0_dp,m=2.0_dp
fx=B*q*vy/m
fy=-B*q*vx/m
end subroutine fmag


subroutine fharm(vx,vy,fx,fy)
real(kind=dp),intent(in)::vx,vy
real(kind=dp),intent(out)::fx,fy
fx=1.0_dp
fy=1.0_dp
end subroutine fharm

subroutine forb(vx,vy,fx,fy)
real(kind=dp),intent(in)::vx,vy
real(kind=dp),intent(out)::fx,fy
fx=1.0_dp
fy=1.0_dp
end subroutine forb

end module