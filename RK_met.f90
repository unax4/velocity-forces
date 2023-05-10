program RK_metodoa

use tipoak
use RK_fun
use rk4

integer:: i, j, k
real(kind=dp), dimension(4):: q
real(kind=dp):: t, h
real(kind=dp), parameter:: t_amai=6.192163_dp, pi=acos(-1.0_dp)
!real(kind=dp), dimension(7), parameter:: tita=(2.0*pi/360)*(/15.0_dp,20.0_dp,25.0_dp,30.0_dp,35.0_dp,40.0_dp,45.0_dp/)
integer, parameter:: n=100000
!h=t_amai/n
h=0.0001_dp

!do k=1,7
t=0.0_dp
q=(/1.2_dp, 0.0_dp, 0.0_dp, -1.049357509830_dp/)

open (unit=20, file="rk_emaitzak.dat", status="replace", action="write")

do i=1,n
 if (t>t_amai) then
         exit
 end if
 write(unit=20, fmt="(2f16.8)"), q(1), q(2)
 call rk4_paso_dp(t,q,forbi,h)
end do
!write (unit=20, fmt="(f16.8)"), q(1)
!end do

close (unit=20)




end program RK_metodoa      
