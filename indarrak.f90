        module indarra

        public :: IM, TP, IH, IG
        contains

        function IM(v)result(a)
        use tipoak
        real(kind=dp), dimension(:), intent(in) :: v
        real(kind=dp), dimension(size(v)) :: a
        real(kind=dp), parameter :: B = -2.0_dp*acos(-1.0_dp), q=-1.0_dp
        real(kind=dp) :: para
        
        para = q*B
        a(1) = -v(2)*para
        a(2) = v(1)*para
    
        end function IM

        function TP(v)result(a)
        use tipoak
        real(kind=dp), dimension(:), intent(in) :: v
        real(kind=dp), dimension(size(v)) :: a
        real(kind=dp), parameter :: g = 9.81_dp, b=2.0_dp
        
        a(1) = -b*v(1)
        a(2) = -g-b*v(2)
        
        end function TP

        function IH(x,v,t)result(a)
        use tipoak
        real(kind=dp), dimension(:), intent(in) :: x,v
        real(kind=dp), intent(in) :: t
        real(kind=dp), dimension(size(v)) :: a        
        real(kind=dp) :: w0,ganma,F,mu,w
        real(kind=dp), parameter :: pi=acos(-1.0_dp)
        ganma = 1.0_dp
        w = pi
        w0 = sqrt(4.0_dp*pi**2+1.0_dp)
        F = sqrt((3.0_dp*pi**2)**2+(2.0_dp*ganma*pi)**2)
        mu = -atan(2.0_dp*ganma*pi/(3.0_dp*pi**2))
        a = -w0**2*x-2.0_dp*ganma*v+F*cos(w*t-mu)
        end function IH

        
        function IG(r,v)result(a)
        use tipoak
        real(kind=dp), dimension(:),intent(in) :: r,v
        real(kind=dp), dimension(size(v)) :: a
        real(kind=dp), parameter :: mu=81.45_dp/82.45_dp, mup=1.0_dp/82.45_dp

      a(1)=r(1)+2.0_dp*v(2)-mu*(r(1)+mup)/(((r(1)+mup)**2+&
      r(2)**2)**1.5_dp)-mup*(r(1)-mu)/(((r(1)-mu)**2+r(2)**2)**(1.5_dp))
      a(2)=r(2)-2.0_dp*v(1)-mu*r(2)/(((r(1)+mup)**2+r(2)**2)**(1.5_dp))&
      -mup*r(2)/(((r(1)-mu)**2+r(2)**2)**(1.5_dp))

        end function IG


        end module indarra

