    program ERM

    use tipoak
    use funtzio_bektorialak
    use mcf_interpoli
!----------------------------------------------------- IZENDAPENA
    ! Indar motaren arabera aldatu beharreko parametroak: n, eps, t_amai, q_m, q(i)-ak

    real(kind=dp), dimension(4) :: q, s, q_m, q_prima, s_prima, delta_bek, ald_q
    real(kind=dp), dimension(2):: t_interp, x_interp, y_interp, vx_interp, vy_interp !Interpolaketa egiteko
    real(kind=dp):: amai_x,amai_y,amai_vx,amai_vy,errorea,amaierako_puntua ! Interpolatuko dugun puntua; kasu bakoitzean ezberdina.
    integer:: k !Zikloak zenbatzeko
    real (kind=dp):: t, h, delta
    integer:: i, j, lurra
    integer, parameter::n= 3000 !3000 harmonikoan, 1386 orbitan
    real(kind=dp), dimension(2):: indarra
    real (kind=dp), parameter:: eps=0.02_dp, masa=1.0_dp, karga=-1.0_dp, t_amai=10.0_dp ! t_amai=10.0 harmonikoan
    real (kind=dp), parameter:: mu=81.45_dp/82.45_dp, mup=1.0_dp/82.45_dp
    real (kind=dp):: ald_c
    real(kind=dp), parameter::pi=acos(-1.0_dp)
    
    k=0
    q_m=10.0_dp ! 1.0 g eta mag, 10.0 harm
   h=t_amai/n
   ! h=0.02_dp
    open (unit=12, file="erm_x_y.dat", action="write", status="replace")
    !open (unit=22, file="erm_t_x.dat", action="write", status="replace")
    !open (unit=32, file="erm_t_y.dat", action="write", status="replace")
    
    t=0.0_dp
    lurra=0 ! Lurra noiz ikutzen duen ikusteko
    q(1)=101.0_dp !x
    q(2)=0.0_dp !y
    q(3)=-100.0_dp
    q(4)=0.0_dp
    !q(3)=30.0_dp*cos((2.0_dp*pi/360)*45.0_dp) !vx
    !q(4)=30.0_dp*sin((2.0_dp*pi/360)*45.0_dp) !vy

    ald_q=q
   ! ald_c=((q(3)**2+q(4)**2-q(1)**2-q(2)**2)/2)-(mu/sqrt(q(2)**2+(q(1)+mup)**2))-(mup/sqrt(q(2)**2+(q(1)-mu)**2))
    
!----------------------------------------------------- PROGRAMA
    do 

    delta=0.0_dp
    if (lurra==1 .or. t>t_amai) then
            exit
    end if

    call fharm (karga, masa, t, q, lurra, indarra)
    s(1)=q(3) !vx
    s(2)=q(4) !vy
    s(3)=indarra(1) !ax
    s(4)=indarra(2) !ay

    q_prima=q+s*(h/2)
    call fharm (karga, masa, (t+(h/2)), q_prima, lurra, indarra)
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
       write (unit=12, fmt="(2f16.8)"),t, q(1)! q(1), q(2)
     
       k=k+1
       t_interp(1)=t
       x_interp(1)=q(1)
       y_interp(1)=q(2)
       vx_interp(1)=q(3)
       vy_interp(1)=q(4)
      
       t=t+h
       q=q+s_prima*h

       t_interp(2)=t
       x_interp(2)=q(1)
       y_interp(2)=q(2)
       vx_interp(2)=q(3)
       vy_interp(2)=q(4)
    end if

    h=0.9_dp*h*sqrt(eps/delta)
   
    end do

    call polint (t_interp, x_interp, 2, t_amai, amai_x, errorea)
    q(1)=amai_x
    call polint (t_interp, y_interp, 2, t_amai, amai_y, errorea)
    q(2)=amai_y
    call polint (t_interp, vx_interp, 2, t_amai, amai_vx, errorea)
    q(3)=amai_vx
    call polint (t_interp, vy_interp, 2, t_amai, amai_vy, errorea)
    q(4)=amai_vy
    ald_q=(q-1.0_dp)/(ald_q-1.0_dp)
   ! ald_c=(((q(3)**2+q(4)**2-q(1)**2-q(2)**2)/2)-(mu/sqrt(q(2)**2+(q(1)+mup)**2))-(mup/sqrt(q(2)**2+(q(1)-mu)**2)))-ald_c
    do i=1,4
    print*, ald_q(i)
    end do

   ! print*, ald_c
!---------------------------------------------------- INTERPOLAKETA

   ! call polint (y_interp, x_interp, 2, 0.0_dp, amaierako_puntua, errorea)
   ! print*, amaierako_puntua
  !  write (unit=12, fmt="(2f16.8)"), amaierako_puntua, 0.0_dp !Bakarrik fg eta fg2-n

    close (unit=12)
    !close (unit=22)
    !close (unit=32)
   
    end program ERM
