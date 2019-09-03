module EDOs
    implicit none
    
    integer, parameter :: cant_ec = 1
contains

    function f(v)
        intent (in) :: v
        real(8) v(0:cant_ec), f
        
        f = v(0)**2     ! Cambiar la funcion aca
    end function f

    function fp(v)      ! Cambiar funcion aca
        intent (in) :: v
        real(8), dimension(0:cant_ec) :: v, fp

        fp(0) = 1       ! Variable independiente
        fp(1) = 2.*v(0) 
    end function fp

    function eulerSimple(v, h)
        intent (in) :: v, h
        real(8), dimension(0:cant_ec) :: v, eulerSimple
        real(8) h
        
        eulerSimple = v + h*fp(v)
    end function eulerSimple

    function eulerModificado(v, h)
        intent (in) :: v, h
        real(8) h
        real(8), dimension(0:cant_ec) :: v, incr1, incr2, eulerModificado

        incr1 = h*fp(v)
        incr2 = h*fp(v + incr1)
        eulerModificado = v + (incr1 + incr2)/2.
    end function eulerModificado

    function eulerMejorado(v, h)
        intent (in) :: v, h
        real(8), dimension(0:cant_ec) :: v, incr, eulerMejorado
        real(8) h

        incr = (h/2.)*fp(v)        
        eulerMejorado = v + h*fp(v + incr)
    end function eulerMejorado

    function rk(v, h)
        intent (in) :: v, h
        real(8), dimension(0:cant_ec):: k1, k2, k3, k4, rk, v
        real(8) h

        k1 = h*fp(v)
        k2 = h*fp(v + k1/2.)
        k3 = h*fp(v + k2/2.)
        k4 = h*fp(v + k3)
        rk = v + (k1 + 2.*k2 + 2.*k3 + k4) / 6.
    end function rk

    function rkf(v, h)
        intent (in) :: v, h
        real(8), dimension(0:cant_ec) :: v, k1, k2, k3, k4, k5, k6, rkf, e
        real(8) h
            
        k1 = h*fp(v)
        k2 = h*fp(v + k1/4.0)
        k3 = h*fp(v + (3.0*k1 + 9.0*k2)/32.0)
        k4 = h*fp(v + (1932.0*k1 - 7200.0*k2 + 7296.0*k3)/2197.0)
        k5 = h*fp(v + 439.0*k1/216.0 - 8.0*k2 + 3680.0*k3/513.0 - 845.0*k4/4104.0)
        k6 = h*fp(v - 8.0*k1/27.0 + 2.0*k2 - 3544.0*k3/2565.0 + 1859.0*k4/4104.0 - 11.0*k5/40.0)
        rkf = v + (25.0*k1/216.0 + 1408.0*k3/2565.0 + 2197.0*k4/4104.0 - k5/5.0)
        e = k1/360.0 - 128.0*k3/4275.0 - 2197.0*k4/75240.0 + k5/50.0 + 2.0*k6/55.0
    end function rkf

end module EDOs

program principal
    use EDOs

    implicit none


    integer(4) i, repeticiones
    real(8) v(0:cant_ec), h

    repeticiones = 15
    h = 1
    v(0) = -5
    v(1) = 25
    open(2, file='datos.dat')
    write(2, *) v
    do i = 1, repeticiones
        v = rkf(v, h)
        write(2, *) v
    end do
    close(2)
    call system('gnuplot -persist valores.p')
end program principal