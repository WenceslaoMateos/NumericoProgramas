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

        fp = 2. * x !Cambiar funcion aca
    end function fp

    function eulerSimple(x, y, h)
        intent (in) :: x, y, h
        real(8) x, h, y, eulerSimple
        
        eulerSimple = y + h * fp(x, y)
    end function eulerSimple

    function eulerModificado(x, y, h)
        intent (in) :: x, y, h
        real(8) x, h, y, yaux, eulerModificado

        yaux = y + h * fp(x, y)        
        eulerModificado = y + h * ((fp(x, y) + fp(x + h, yaux)) / 2.)
    end function eulerModificado

    function eulerMejorado(x, y, h)
        intent (in) :: x, y, h
        real(8) x, h, y, yaux, eulerMejorado

        yaux = y + (h / 2.) * fp(x, y)        
        eulerMejorado = y + h * fp(x + (h / 2.), yaux)
    end function eulerMejorado

    function rk(x, y, h)
        intent (in) :: x, y, h
        real(8) x, h, y, k1, k2, k3, k4, rk

        k1 = h * fp(x, y)
        k2 = h * fp(x + h / 2., y + k1 / 2.)
        k3 = h * fp(x + h / 2., y + k2 / 2.)
        k4 = h * fp(x + h, y + k3)
        rk = y + (1. / 6.) * (k1 + 2. * k2 + 2. * k3 + k4)
    end function rk

    function rkf(x, y, h)
        intent (in) :: x, y, h
        real(8) x, h, y, k1, k2, k3, k4, k5, k6, rkf

        k1 = h * fp(x, y)
        k2 = h * fp(x + h/4., y + k1/4.)
        k3 = h * fp(x + 3.*h/8., y + (3.*k1 + 9.*k2)/32.)
        k4 = h * fp(x + 12.*h/13., y + (1932.*k1 - 7200.*k2 + 7296.*k3)/2197)
        k5 = h * fp(x + h, y + (439.*k1/216.) - 8.*k2 + (3680.*k3/513.) - (845.*k4/4104.))
        k6 = h * fp(x + h/2., y - 8.*k1/27. + 2.*k2 - 3544.*k3/2565. - 1859.*k4/4104. - 11.*k5/40.)
        rkf = y + 25.*k1/216. + 1408.*k3/2565. + 2197.*k4/4104. - k5/5.
    end function rkf

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
        write(2,'(2F15.10)') x, y
        y = eulerSimple(x, y, h)
        x = x + h
    end do
    close(2)
    call system('gnuplot -persist valores.p')
end program principal