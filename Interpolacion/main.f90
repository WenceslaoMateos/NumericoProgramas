program main
    use interpolacion

    implicit none
    
    real(8), dimension(0:3) :: x, y
    real(8), dimension(0:3) :: lagrange
    
    x = [-4._8, -2._8, 0._8, 2._8]
    y = [-7.38_8, 0.52_8, 2._8, 14.52_8]
    lagrange = polinomioLagrange(x, y)
    write(*, *) lagrange
end program main