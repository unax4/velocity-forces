module funtzio_bektorialak 
use tipoak
public::fg,fmag,fharm,forbi

contains

subroutine fg(fkarga, fmasa, ft, fq, flurra, findarra)

        use tipoak
real(kind=dp),intent(in):: fkarga, fmasa, ft
integer, intent(inout)::flurra
real(kind=dp), intent(in), dimension(:):: fq
real(kind=dp),intent(inout), dimension(:)::findarra
real(kind=dp),parameter::b=1.5_dp, g=9.8_dp

if (fq(2)<=0.0_dp .and. fq(4)<0.0_dp) then
        flurra=1
end if

findarra(1)=-b*fq(3)
findarra(2)=-b*fq(4) - g

end subroutine fg



subroutine fmag (fkarga, fmasa, ft, fq, flurra, findarra)

   use tipoak
   integer, intent(inout):: flurra
   real (kind=dp), intent(in):: fkarga, fmasa, ft
   real (kind=dp), parameter:: eremua=1.0_dp
   real (kind=dp), intent(in), dimension(:)::fq
   real (kind=dp), intent(inout), dimension(:):: findarra

   findarra(1)=(fkarga/fmasa)*fq(4)*eremua
   findarra(2)=-(fkarga/fmasa)*fq(3)*eremua

end subroutine fmag 


subroutine fharm(fkarga, fmasa, ft, fq, flurra, findarra)
use tipoak
integer, intent(inout)::flurra
real(kind=dp),intent(in):: ft, fkarga, fmasa
real(kind=dp), intent(in), dimension(:):: fq
real(kind=dp),intent(inout), dimension(:):: findarra
real(kind=dp),parameter::pi=acos(-1.0_dp)

real(kind=dp), parameter:: c1=100.0_dp, c2=0.0_dp, w1=2*pi, gamba=1.0_dp, A=1.0_dp, w=pi, psi=0.0_dp
real (kind=dp):: w0, F, nu
w0= sqrt(w1**2+gamba**2)
F= A*sqrt((w0**2-w**2)**2+(2.0_dp*gamba*w)**2)
nu= psi - atan((2.0_dp*gamba*w)/(w0**2-w**2))

!findarra(1)=-(w0**2)*fq(1) - 2.0_dp*gamba*fq(3) + F*cos(w*ft-nu)
findarra(2)=-(w0**2)*fq(2) - 2.0_dp*gamba*fq(4) + F*cos(w*ft-nu)
findarra(1)=0.0_dp

end subroutine fharm


subroutine forbi(fkarga, fmasa, ft, fq, flurra, findarra)

        use tipoak
integer, intent(inout)::flurra
real(kind=dp),intent(in):: fkarga, fmasa, ft
real (kind=dp), intent(in), dimension(:)::fq
real(kind=dp), intent(inout), dimension(:)::findarra
real(kind=dp), parameter:: pi=acos(-1.0_dp), mu=81.45_dp/82.45_dp, mup=1.0_dp/82.45_dp

findarra(1)=fq(1)+2.0_dp*fq(4)-&
        mu*(fq(1)+mup)/((fq(1)+mup)**2+fq(2)**2)**(3.0_dp/2)-mup*(fq(1)-mu)/((fq(1)-mup)**2+fq(2)**2)**(3.0_dp/2)
findarra(2)= fq(2) - 2.0_dp*fq(3) - mu*fq(2)/((fq(1)+mup)**2+fq(2)**2)**(3.0_dp/2) - mup*fq(2)/((fq(1)-mup)**2+fq(2)**2)**(3.0_dp/2)

end subroutine forbi

end module funtzio_bektorialak
