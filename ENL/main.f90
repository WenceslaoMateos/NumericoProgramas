program main
    use ENL
    implicit none
    real(8) a, b, tol, res, x0, lambda
    
    !call plot("cos(x)*cosh(x)+1", "[4:5]")
    a = -0.3
    b = -0.1
    tol = 1e-5
    res = biseccion(a, b, func, tol)
    write(*, *) "El resultado de biseccion es: ", res
    
    x0 = -0.99
    tol = 1e-5
    res = puntoFijo(x0, func, g, tol, 200)
    write(*, *) "El resultado de punto fijo es: ", res

    a = -0.3
    b = -0.1
    x0 = -0.99
    tol = 1e-5
    !El lambda a elegir depende del intervalo a buscar raiz.
    lambda = 1./maxDeriv(a, b, dfunc, tol)
    res = puntoFijoSistematico(x0, func, lambda, tol, 200)
    write(*, *) "El resultado de punto fijo sistematico es: ", res

    x0 = -0.3
    tol = 1e-5
    res = newton(x0, func, dfunc, tol, 200)
    write(*, *) "El resultado de newton es: ", res

    !call plot("((x+1)/(x-1))-sin(3*x)", "f(x)", "[-5:5]")
    !call plot("asin((x+1.)/(x-1.))/3.", "g(x)", "[-5:5]")
    call plot("-2./(x-1)**2-3.*cos(3*x)", "Df(x)", "[-0.3:-0.1]")


    contains

    function func(x)
        real(8), intent(in) :: x
        real(8) func

        func = ((x+1)/(x-1)) - sin(3*x)
    end function func

    function dfunc(x)
        real(8), intent(in) :: x
        real(8) dfunc

        dfunc = -2./(x-1)**2 - 3. *cos(3*x)
    end function dfunc

    function g(x)
        real(8), intent(in) :: x
        real(8) g

        g = asin((x+1.)/(x-1.))/3.
    end function g

    function maxDeriv(a, b, deriv, h)
        real(8), intent(in) :: a, b, h
        real(8) maxDeriv, x, act
        procedure(funcion) :: deriv

        maxDeriv = deriv(a)
        x = a + h
        do while (x <= b)
            act = deriv(x)
            if (act > maxDeriv) then
                maxDeriv = act
            end if
            x = x + h
        end do
    end function maxDeriv

end program main