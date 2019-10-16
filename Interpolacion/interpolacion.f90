module interpolacion

    implicit none
    
contains
    
    function productoPolinomios(p1, p2)
        real(8), dimension(0:), intent(in) :: p1, p2
        real(8) productoPolinomios(0:ubound(p1, 1) + ubound(p2, 1))
        integer(4) i, j

        productoPolinomios = 0.
        do i = 0, ubound(p1, 1)
            do j = 0, ubound(p2, 1)
                productoPolinomios(i + j) = productoPolinomios(i + j) + p1(i) * p2(j)
            end do
        end do
    end function productoPolinomios

    function polinomioLagrange(x, y)
        real(8), dimension(0:), intent(in) :: x, y
        real(8), dimension(0:ubound(x, 1)) :: polinomioLagrange
        real(8), dimension(0:1) :: binomio
        real(8), dimension(:), allocatable :: ant, act
        real(8) denominador
        integer(4) orden, k, i

        orden = size(x) - 1
        polinomioLagrange = 0.
        do k = 0, orden
            denominador = 1.
            act = [1.]
            do i = 0, orden
                if (i /= k) then
                    ant = act
                    deallocate(act)
                    allocate(act(0:size(ant)))
                    binomio = [-x(i), 1._8]
                    denominador = denominador * (x(k) - x(i))
                    act = productoPolinomios(ant, binomio)
                    deallocate(ant)
                end if
            end do
            polinomioLagrange = polinomioLagrange + act * y(k) / denominador
            deallocate(act)
        end do
    end function polinomioLagrange

end module interpolacion