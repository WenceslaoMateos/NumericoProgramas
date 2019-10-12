program principal
    use EDOs

    implicit none

    integer(4), parameter :: cant_ec = 1
    character(len=*), parameter :: archivo = "datos.dat"

    integer(4) repeticiones
    real(8) v(0:cant_ec), h, tol, xf

    repeticiones = 15
    h = 1.
    v(0) = -5.
    v(1) = 25.
    tol = 0.01
    xf = 20.
    call iterar(eulerSimple, fp, v, h, estrategia1, tol, xf, archivo)
    call plot(archivo)
    call iterar(eulerModificado, fp, v, h, estrategia1, tol, xf, archivo)
    call plot(archivo)
    call iterar(eulerMejorado, fp, v, h, estrategia1, tol, xf, archivo)
    call plot(archivo)
    call iterar(rk, fp, v, h, estrategia1, tol, xf, archivo)
    call plot(archivo)
    call iterar(rkf, fp, v, h, estrategia2, tol, xf, archivo)
    call plot(archivo)
    
contains

    function f(v)
        intent (in) :: v
        real(8) v(0:cant_ec), f
        
        f = v(0)**2     ! Cambiar la funcion aca
    end function f

    function fp(v)      ! Cambiar funcion aca
        intent (in) :: v
        real(8) v(0:), fp(0:size(v) - 1)

        fp(0) = 1       ! Variable independiente
        fp(1) = 2.*v(0) 
    end function fp

    subroutine plot(archivo)
        intent(in) :: archivo
        character (LEN=*) :: archivo

        open(unit=2, file="temporal.p", access='SEQUENTIAL', status='REPLACE')
        write(2, *) "set autoscale"
        write(2, *) "unset log"
        write(2, *) "unset label"
        write(2, *) "set xtic auto"
        write(2, *) "set ytic auto"
        write(2, *) "set title 'E2'"
        write(2, *) "set xlabel 'x'"
        write(2, *) "set ylabel 'f(x)'"
        write(2, *) "plot '", archivo, "' using 1:2 title 'metodo' with linespoints"
        call system('gnuplot -persist temporal.p')
        close(2, STATUS='DELETE')
    end subroutine plot
end program principal