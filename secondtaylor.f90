program taylor
use tipoak
real(kind=dp)::t,ta,tb,vint,xint,yint,fint,h
real(kind=dp),dimension(2)::posi,velo,velint,for,forint
integer::i,j,k
integer,parameter::n=500


!bektore gabe
ta=0_dp
tb=10_dp
t=ta
h=(tb-ta)/n


posi=(/0.0_dp,4.0_dp/)
velo=(/10.0_dp,10.0_dp/)



do i=1,n
write(unit=*,fmt="(3f12.4)")t,posi

!x=x+h*n*v+f*h**2/2
!vint=v+f*h
!fint=force(t+h,x,vint)

!v=v+(f+fint)*h/2
!t=t+h
!f=force(t,x,v) !aire erresistentzia problementzat azken  jarri, f txikia bada indar kontserbatzaileekin konparatuz kendu eta f=fint)

call force(t,posi(1),posi(2),velo(1),velo(2),for(1),for(2))
posi=posi + h*velo + for*h**2/2

velint=velo + for*h
call force(t+h,posi(1),posi(2),velint(1),velint(2),forint(1),forint(2))

velo=velo + (for+forint)*h/2
t=t+h

if (posi(2)<0) then

exit
end if

end do

contains
subroutine force(temp,x,y,vx,vy,fx,fy)
real(kind=dp),intent(in)::temp,x,y,vx,vy
real(kind=dp),intent(out)::fx,fy
real(kind=dp),parameter::c=-0.3_dp,m=2_dp,g=9.8_dp

fx=-c*vx
fy=-c*vy - g
end subroutine force

end program taylor
