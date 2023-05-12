module RK_fun

public::fg, fmag, fharm, forbi,fg2

contains

function fg (t, q) result (s)

use tipoak

real (kind=dp), intent(in):: t
real (kind=dp), intent(in), dimension(:):: q
real (kind=dp), dimension(size(q)):: s
real (kind=dp), parameter:: b=2.0_dp, g=9.8_dp

!if (q(2)<=0.0_dp .and. q(4)<0.0_dp) then
        s=0.0_dp
!else

s(1)=q(3)
s(2)=q(4)
s(3)=-b*s(1)
s(4)=-b*s(2) - g

!end if
end function fg

function fg2 (t, q) result (s)

use tipoak

real (kind=dp), intent(in):: t
real (kind=dp), intent(in), dimension(:):: q
real (kind=dp), dimension(size(q)):: s
real (kind=dp):: v_mod
real (kind=dp), parameter:: beta=0.3_dp, g=9.8_dp

!if (q(2)<=0.0_dp .and. q(4)<0.0_dp) then
        s=0.0_dp
!else

s(1)=q(3)
s(2)=q(4)
v_mod=sqrt(s(1)**2+s(2)**2)
s(3)=-beta*v_mod*s(1)
s(4)=-beta*v_mod*s(2) - g

!end if
end function fg2


function fmag (t, q) result (s)

use tipoak

real (kind=dp), intent(in):: t
real (kind=dp), intent(in), dimension(:):: q
real (kind=dp), dimension(size(q)):: s
real (kind=dp), parameter:: pi= acos(-1.0_dp), eremua=pi*2.0_dp, karga=-1.0_dp, masa=1.0_dp

s(1)=q(3)
s(2)=q(4)
s(3)=(karga/masa)*s(2)*eremua
s(4)=-(karga/masa)*s(1)*eremua

end function fmag

function fharm (t, q) result (s)

use tipoak

real (kind=dp), intent(in):: t
real (kind=dp), intent(in), dimension(:):: q
real (kind=dp), dimension(size(q)):: s

real(kind=dp), parameter:: pi=acos(-1.0_dp), c1=100.0_dp, c2=0.0_dp, w1=2*pi, gamba=1.0_dp, A=1.0_dp, w=pi, psi=0.0_dp
real (kind=dp):: w0, F, nu
w0= sqrt(w1**2+gamba**2)
F= A*sqrt((w0**2-w**2)**2+(2.0_dp*gamba*w)**2)
nu= psi - atan((2.0_dp*gamba*w)/(w0**2-w**2))

s(1)=q(3)
s(2)=q(4)
s(3)=-(w0**2)*q(1) - 2.0_dp*gamba*s(1) + F*cos(w*t-nu)
!s(4)=-(w0**2)*q(2) - 2.0_dp*gamba*s(2) + F*cos(w*t-nu)
s(4)=0.0_dp

end function fharm

function forbi (t, q) result (s)

use tipoak

real (kind=dp), intent(in):: t
real (kind=dp), intent(in), dimension(:):: q
real (kind=dp), dimension(size(q)):: s

real(kind=dp), parameter:: pi=acos(-1.0_dp), mu=81.45_dp/82.45_dp, mup=1.0_dp/82.45_dp

s(1)=q(3)
s(2)=q(4)
!s(3)=q(1)+2.0_dp*s(2)-&
      !  mu*(q(1)+mup)/((q(1)+mup)**2.0_dp+q(2)**2.0_dp)**(3.0_dp/2.0_dp)-mup*(q(1)-mu)/((q(1)-mup)**2.0_dp+q(2)**2.0_dp)&
       ! **(3.0_dp/2.0_dp)
!s(4)= q(2) - 2.0_dp*s(1) - mu*q(2)/((q(1)+mup)**2.0_dp+q(2)**2.0_dp)**(3.0_dp/2.0_dp)&
        !- mup*q(2)/((q(1)-mup)**2.0_dp+q(2)**2.0_dp)**(3.0_dp/2.0_dp)

s(3)=q(1)+2.0_dp*s(2)-mu*(q(1)+mup)/((q(1)+mup)**2+q(2)**2)**(3.0_dp/2)&
        -mup*(q(1)-mu)/((q(1)-mu)**2+q(2)**2)**(3.0_dp/2)
s(4)=q(2)-2.0_dp*s(1)-mu*q(2)/((q(1)+mup)**2+q(2)**2)**(3.0_dp/2)-&
        mup*q(2)/((q(1)-mu)**2+q(2)**2)**(3.0_dp/2)

end function forbi


end module RK_fun
