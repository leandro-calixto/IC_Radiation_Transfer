subroutine SLW1(x, intensity, n, I0, k) 
    implicit none
  
    integer, intent(in) :: n
    real, intent(in) :: x(n), I0, k
    real, intent(out) :: intensity(n) 
    integer :: j  
  
    do j = 1, n
      intensity(j) = I0 * exp(-k * x(j))
    end do
  
end subroutine SLW1
  
subroutine gerar_distancias(x, n, x0, xf)
    implicit none

    integer, intent(in) :: n
    real, intent(in) :: x0, xf
    real, intent(out) :: x(n)
    integer :: i

    do i = 1, n
        x(i) = x0 + (xf - x0) * real(i - 1) / real(n - 1)
    end do

end subroutine gerar_distancias