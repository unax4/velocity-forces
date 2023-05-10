    program ERM

    use tipoak
    use funtzio_bektorialak

    ! Indar motaren arabera aldatu beharreko parametroak: n, eps, t_amai, q_m, q(i)-ak

    real(kind=dp), dimension(4) :: q, s, q_m, q_prima, s_prima, delta_bek
    real (kind=dp):: t, h, delta
    integer:: i, j, lurra, k
    integer, parameter::n= 4318 !3000 harmonikoan, 1386 orbitan
    real(kind=dp), dimension(2):: indarra
    real (kind=dp), parameter:: eps=0.0001_dp, masa=1.0_dp, karga=1.0_dp, t_amai=6.192169331396_dp ! t_amai=10.0 harmonikoan
    real(kind=dp), parameter::pi=acos(-1.0_dp)
    !G2!real(kind=dp), dimension(7), parameter:: tita=(2.0_dp*pi/360)*(/15.0_dp,20.0_dp,25.0_dp,30.0_dp,35.0_dp,40.0_dp,45.0_dp/)

    q_m=1.0_dp ! 1.0 g eta mag, 10.0 harm
    h=t_amai/n
    open (unit=12, file="emaitzak.dat", action="write", status="replace")
    !G2!do k=1,7
    t=0.0_dp
    lurra=0 ! Lurra noiz ikutzen duen ikusteko
    q(1)=1.2_dp !x
    q(2)=0.0_dp !y
    q(3)=0.0_dp !vx
    q(4)=-1.0493575098300_dp !vy
    

    do 

    delta=0.0_dp
    if (lurra==1 .or. t>t_amai) then
            exit
    end if

    call forbi (karga, masa, t, q, lurra, indarra)
    s(1)=q(3) !vx
    s(2)=q(4) !vy
    s(3)=indarra(1) !ax
    s(4)=indarra(2) !ay

    q_prima=q+s*(h/2)
    call forbi (karga, masa, (t+(h/2)), q_prima, lurra, indarra)
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
print*, t
   !G2!write (unit=12, fmt="(f16.8)"), q(1)
   !G2!end do
    close (unit=12)


    end program ERM
