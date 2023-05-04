    program ERM

    use tipoak
   ! use funtzioak

    real(kind=dp), dimension(4) :: q, s, q_m, q_prima, s_prima, delta_bek
    real (kind=dp):: t, h, delta
    integer:: i, j
    integer, parameter::n=15000
    real(kind=dp), dimension(2):: indarra
    real (kind=dp), parameter:: eps=0.01_dp, masa=1.0_dp, karga=1.0_dp
    q_m=1.0_dp
    t=0.0_dp
    h=0.02_dp
    delta=0.0_dp
    open (unit=12, file="emaitzak.dat", action="write", status="replace")
    q(1)=-1.0_dp !x
    q(2)=0.0_dp !y
    q(3)=0.0_dp !vx
    q(4)=1.0_dp !vy
    

    do i=1, n


    call indar_magnetikoa (karga, masa, t, q, indarra)
    s(1)=q(3) !vx
    s(2)=q(4) !vy
    s(3)=indarra(1) !ax
    s(4)=indarra(2) !ay

    q_prima=q+s*(h/2)
    call indar_magnetikoa (karga, masa, (t+(h/2)), q_prima, indarra)
    s_prima(1)=q_prima(3)
    s_prima(2)=q_prima(4)
    s_prima(3)=indarra(1)
    s_prima(4)=indarra(2)
   
    do j=1,4
    delta_bek(j)= abs(((s_prima(j)-s(j))*h)/(2*q_m(j)))
    if (delta_bek(j)>delta) then
            delta=delta_bek(j)
    end if
    end do

    if (delta<eps) then
       write (unit=12, fmt="(2f16.8)"), q(1), q(2)
       t=t+h
       q=q+s_prima*h
    end if
    h=0.9_dp*h*sqrt(eps/delta)

    end do
    close (unit=12)

    contains

    subroutine indar_magnetikoa (fkarga, fmasa, ft, fq, findarra)

       use tipoak
       real (kind=dp), intent(in):: fkarga, fmasa, ft
       real (kind=dp):: eremua
       real (kind=dp), intent(in), dimension(:)::fq
       real (kind=dp), intent(inout), dimension(:):: findarra
       eremua=1.0_dp

       findarra(1)=(fkarga/fmasa)*fq(4)*eremua
       findarra(2)=-(fkarga/fmasa)*fq(3)*eremua

    end subroutine indar_magnetikoa 


    end program ERM
