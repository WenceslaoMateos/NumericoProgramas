module ENL
    implicit none
    
    abstract interface
        function funcion(x)
            real(8), intent(in) :: x
            real(8) funcion
        end function
    end interface

contains
    
    function biseccion(a, b, f, tol)
        real(8), intent(in) :: a, b, tol
        procedure(funcion) :: f
        real(8) biseccion, izq, der
        integer(4) n, i

        izq = a
        der = b
        n = ceiling((log(abs(der - izq)/ tol))/log(2.))
        do i = 1, n
            biseccion = (izq + der)/2.
            if (f(izq) * f(biseccion) < 0.) then
                der = biseccion
            else
                izq = biseccion
            end if
        end do
    end function biseccion

    function puntoFijo(x0, f, g, tol, max_iter)
        real(8), intent(in) :: x0, tol
        integer(4), intent(in) :: max_iter
        procedure(funcion) :: f, g
        integer(4) i
        real(8) error, x, puntoFijo
        
        x = x0
        error = tol + 1
        i = 0
        do while ((error >= tol) .AND. (i <= max_iter))
            x = g(x)
            error = abs(f(x))
            i = i + 1
        end do
        puntoFijo = x
    end function puntoFijo

    function puntoFijoSistematico(x0, f, lambda, tol, max_iter)
        real(8), intent(in) :: x0, tol
        integer(4), intent(in) :: max_iter
        procedure(funcion) :: f
        integer(4) i
        real(8) error, x, puntoFijoSistematico, lambda
        
        x = x0
        error = tol + 1
        i = 0
        do while ((error >= tol) .AND. (i <= max_iter))
            x = x - lambda * f(x)
            error = abs(f(x))
            i = i + 1
        end do
        puntoFijoSistematico = x
    end function puntoFijoSistematico

    function newton(x0, f, df, tol, max_iter)
        real(8), intent(in) :: x0, tol
        integer(4), intent(in) :: max_iter
        procedure(funcion) :: f, df
        integer(4) i
        real(8) error, x, newton
        
        x = x0
        error = tol + 1
        i = 0
        do while ((error >= tol) .AND. (i <= max_iter))
            x = x - f(x)/df(x)
            error = abs(f(x))
            i = i + 1
        end do
        newton = x
    end function newton

    subroutine plot(f, nombre, rango)
        character (LEN=*), intent(in) :: f, rango, nombre
        optional rango

        open(unit=2, file="temporal.p", access='SEQUENTIAL', status='REPLACE')
        write(2, *) "set autoscale"
        write(2, *) "unset log"
        write(2, *) "unset label"
        write(2, *) "set grid"
        write(2, *) "set xtic auto"
        write(2, *) "set ytic auto"
        write(2, *) "set title 'Spline'"
        write(2, *) "set xlabel 'x'"
        write(2, *) "set ylabel 'f(x)'"
        write(2, *) "plot ", rango ," ", f ," title '", nombre ,"' with lines"
        call system('gnuplot -persist temporal.p')
        close(2, STATUS='DELETE')
    end subroutine plot
end module ENL