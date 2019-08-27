module EDOs
    implicit none
    
contains

    function f(x)
        real(8) x, f
        
        f = x**2 !Cambiar la funcion aca
    end function f

    function fp(x, y)
        intent (in) :: x, y
        real(8) x, y, fp

        fp = 2*x !Cambiar funcion aca
    end function fp

    function eulerSimple(x, y, h)
        intent (in) :: x, y, h
        real(8) x, h, y, eulerSimple
        
        eulerSimple = y + h * fp(x, y)
    end function eulerSimple
    
end module EDOs

program principal
    use EDOs

    implicit none


    integer(4) i, repeticiones
    real(8) x, y, h

    repeticiones = 15
    x = -5
    h = 1
    y = 25
    open(2, file='datos.dat')
    do i = 1, repeticiones
        y = eulerSimple(x, y, h)
        write(2,'(2F15.10)') x, y
        x = x + h
    end do
    close(2)
    call system('gnuplot -persist valores.p')
end program principal