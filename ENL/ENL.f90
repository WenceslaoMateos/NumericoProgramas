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

    function biseccionTolY(a, b, f, tol)
        real(8), intent(in) :: a, b, tol
        procedure(funcion) :: f
        real(8) biseccionTolY, izq, der, fmedio

        izq = a
        der = b
        biseccionTolY = (izq + der)/2.
        fmedio = f(biseccionTolY)
        do while (abs(fmedio) >= tol)
            if (f(izq) * fmedio < 0.) then
                der = biseccionTolY
            else
                izq = biseccionTolY
            end if
            biseccionTolY = (izq + der)/2.
            fmedio = f(biseccionTolY)
        end do
    end function biseccionTolY

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

    function puntoFijoSistematico(x0, f, a, b, deriv, tol, max_iter)
        real(8), intent(in) :: x0, a, b, tol
        integer(4), intent(in) :: max_iter
        procedure(funcion) :: f, deriv
        integer(4) i
        real(8) error, x, puntoFijoSistematico, lambda
        
        x = x0
        error = tol + 1
        i = 0
        lambda = 1./maxDeriv(a, b, deriv, (b-a)/100)
        do while ((error >= tol) .AND. (i <= max_iter))
            x = x - lambda * f(x)
            error = abs(f(x))
            i = i + 1
        end do
        puntoFijoSistematico = x
    end function puntoFijoSistematico

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

    function newtonTolX(x0, f, df, tol, max_iter)
        real(8), intent(in) :: x0, tol
        integer(4), intent(in) :: max_iter
        procedure(funcion) :: f, df
        integer(4) i
        real(8) error, x, xant, newtonTolX
        
        x = x0
        error = tol + 1
        i = 0
        do while ((error >= tol) .AND. (i <= max_iter))
            xant = x
            x = x - f(x)/df(x)
            error = abs(x - xant)
            i = i + 1
        end do
        newtonTolX = x
    end function newtonTolX

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
        write(2, *) "set title '", nombre, "'"
        write(2, *) "set xlabel 'x'"
        write(2, *) "set ylabel 'f(x)'"
        write(2, *) "plot ", rango ," ", f ," title '", f ,"' with lines"
        call system('gnuplot -persist temporal.p')
        close(2, STATUS='DELETE')
    end subroutine plot

    subroutine plotVarias(funciones, nombre, rango)
        character(len=*), dimension(:), intent(in) :: funciones
        character (LEN=*), intent(in) :: rango, nombre
        optional rango
        integer(4) i
    
        open(unit=2, file="temporal.p", access='SEQUENTIAL', status='REPLACE')
        write(2, *) "set autoscale"
        write(2, *) "unset log"
        write(2, *) "unset label"
        write(2, *) "set grid"
        write(2, *) "set xtic auto"
        write(2, *) "set ytic auto"
        write(2, *) "set title '", nombre, "'"
        write(2, *) "set xlabel 'x'"
        write(2, *) "set ylabel 'y'"
        write(2, "(7A)", advance="NO") "plot ", rango ," ", funciones(1) ," title '", funciones(1) ,"' with lines"
        do i = 2, size(funciones)
            write(2, *) ", \"
            write(2, "(6A)", advance="NO") funciones(i) ," title '", funciones(i) ,"' with lines"
        end do
        write(2, *)
        call system('gnuplot -persist temporal.p')
        close(2, STATUS='DELETE')
    end subroutine plotVarias
    
end module ENL