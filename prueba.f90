        program prueba

        use indarra
        use tipoak
        use lpa_module

        real(kind=dp), dimension(2) :: rm,rg, vm,vg, am, ag,rp,vp,ap
        real(kind=dp),dimension(1) :: rh, vh, ah
        integer :: i,j,k
        real(kind=dp) :: t, step,steptp,stepim
        
        rm(1)=0.0_dp
        rm(2)=0.0_dp
        rh=3.0_dp
        rg(1)=1.2_dp
        rg(2)=0.0_dp
        rp(1) = 0.0_dp
        rp(2) = 0.0_dp
        vm(1)=2.0_dp*acos(-1.0_dp)
        vm(2)=0.0_dp
        vh = 0.0_dp
        vg(1)=0.0_dp
        vg(2)=-1.049357509830_dp
        vp(1)=20.0_dp
        vp(2)=20.0_dp
        steptp = 0.1_dp
        stepim = 0.01_dp
        step = 0.00001_dp
        t = 0.0_dp
        open(unit=13,file="LPA_IM.dat",status="replace",action="write")
        open(unit=14,file="LPA_IH.dat",status="replace",action="write")
        open(unit=15,file="LPA_IG.dat",status="replace",action="write")
        open(unit=16,file="LPA_TP.dat",status="replace",action="write")
        
        do j=1,400
        am = IM(vm)
        call lpa(rm,vm,am,stepim)
        write (unit=13,fmt="(2es16.8)"), rm(1), rm(2)
        end do

        do k=1,1000000
        ag = IG(rg,vg)
        call lpa(rg,vg,ag,step)
        write(unit=15,fmt="(2es16.8)"), rg(1), rg(2)
        end do

        do i=1,1000
                t = t + step
                ah = IH(rh,vh,t)
                call lpa(rh,vh,ah,step)
                write(unit=14,fmt="(3es16.8)"),t,rh
        end do
        
        t=0.0_dp
        do
        if (rp(2) < 0) then
                exit
        end if
        t=t+steptp
        ap=TP(vp)
        call lpa(rp,vp,ap,steptp)
        write(unit=16,fmt="(2es16.8)"), rp(1),rp(2)
        

        end do
        
        close(unit=13)
        close(unit=14)
        close(unit=15)
        close(unit=16)

        end program prueba

