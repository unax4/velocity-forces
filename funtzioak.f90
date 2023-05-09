module funtzioak 
use tipoak
public::fg,fmag,fharm,forbi

contains
subroutine fg(vx,vy,fx,fy,b)
real(kind=dp),intent(in)::vx,vy,b
real(kind=dp),intent(inout)::fx,fy
real(kind=dp),parameter::m=2.0_dp,g=9.8_dp

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


subroutine fharm(t,x,vx,fx)
real(kind=dp),intent(in)::vx,x,t
real(kind=dp),intent(out)::fx
real(kind=dp),parameter::pi=acos(-1.0_dp)
real(kind=dp), parameter:: c1=100.0_dp, c2=0.0_dp, w1=2*pi, gamba=1.0_dp, A=1.0_dp, w=pi, psi=0.0_dp
real (kind=dp):: w0, F, nu
w0= sqrt(w1**2+gamba**2)
F= A*sqrt((w0**2-w**2)**2+(2.0_dp*gamba*w)**2)
nu= psi - atan((2.0_dp*gamba*w)/(w0**2-w**2))

fx=-(w0**2)*x - 2.0_dp*gamba*vx + F*cos(w*t-nu)

end subroutine fharm

subroutine forbi(x,y,vx,vy,fx,fy)
real(kind=dp),intent(in)::x,y,vx,vy
real(kind=dp),intent(out)::fx,fy
real(kind=dp),parameter::mu=81.45_dp/82.45_dp,mup=1.0_dp/82.45_dp
fx=x+2.0_dp*vy-mu*(x+mup)/((x+mup)**2.0_dp+y**2.0_dp)**(3.0_dp/2.0_dp)&
        -mup*(x-mu)/((x-mu)**2.0_dp+y**2.0_dp)**(3.0_dp/2.0_dp)
     fy=y-2.0_dp*vx-mu*y/((x+mup)**2.0_dp+y**2.0_dp)**(3.0_dp/2.0_dp)-&
        mup*y/((x-mu)**2.0_dp+y**2.0_dp)**(3.0_dp/2.0_dp)
end subroutine forbi

end module funtzioak
