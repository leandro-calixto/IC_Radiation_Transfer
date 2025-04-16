program main
    implicit none
  
    integer, parameter :: n = 10
    real :: I0, k
    real, dimension(n) :: x, intensity  
    integer :: i
  
    I0 = 1.0
    k = 0.5
  
    call gerar_distancias(x, n, 0.0, 9.0)
    call SLW1(x, intensity, n, I0, k)  
  
    print *, 'x     I(x)'
    do i = 1, n
       print '(F6.2, 2X, F6.4)', x(i), intensity(i)
    end do
  
end program main