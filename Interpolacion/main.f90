program main
    use interpolacion

    implicit none
    
    real(8), dimension(0:3) :: x, y
    real(8), dimension(0:3) :: lagrange
    
    x = [-4., -2., 0., 2.]
    y = [-7.38, 0.52, 2., 14.52]
    lagrange = polinomioLagrange(x, y)
    write(*, '(4F10.5)') lagrange
end program main