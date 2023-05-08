        module lpa_module
        use tipoak
        use indarra
        public :: lpa        
                
        contains

        subroutine lpa(r, v, a, step)
        use tipoak
        real(kind=dp), dimension(:), intent(in) :: a
        real(kind=dp), dimension(:), intent(inout) :: r, v
        real(kind=dp), intent(in) :: step
        
        ! Calculate the next state using the LPA formula
        v = v + a*step

        r = r + v*step

        end subroutine lpa

        end module lpa_module

