module EDOs
    implicit none
    
    abstract interface 
        function funcion(v)
            real(8), intent(in) :: v(0:)
            real(8) f
        end function funcion
    end interface

    abstract interface 
        function derivada(v)
            intent(in) :: v
            real(8) v(0:), derivada(0:size(v)-1)
        end function derivada
    end interface

    abstract interface
        function metodo(v, h, fp)
            procedure(derivada) :: fp
            real(8), intent(in) :: v(0:), h
            real(8) metodo(0:size(v)-1)
        end function metodo
    end interface

    interface iterar
        module procedure iterarVeces, iterarValorFinal
    end interface
contains

    function eulerSimple(v, h, fp)
        intent (in) :: v, h
        procedure(derivada) :: fp
        real(8) v(0:), eulerSimple(0:size(v)-1), h
  
        eulerSimple = v + h*fp(v)
    end function eulerSimple

    function eulerModificado(v, h, fp)
        intent (in) :: v, h
        procedure(derivada) :: fp
        real(8) h, v(0:)
        real(8), dimension(0:size(v)-1) :: incr1, incr2, eulerModificado

        incr1 = h*fp(v)
        incr2 = h*fp(v + incr1)
        eulerModificado = v + (incr1 + incr2)/2.
    end function eulerModificado

    function eulerMejorado(v, h, fp)
        intent (in) :: v, h
        procedure(derivada) :: fp
        real(8) h, v(0:)
        real(8), dimension(0:size(v)-1) :: incr, eulerMejorado

        incr = (h/2.)*fp(v)        
        eulerMejorado = v + h*fp(v + incr)
    end function eulerMejorado

    function rk(v, h, fp)
        intent (in) :: v, h
        procedure(derivada) :: fp
        real(8) h, v(0:)
        real(8), dimension(0:size(v)-1):: k1, k2, k3, k4, rk

        k1 = h*fp(v)
        k2 = h*fp(v + k1/2.)
        k3 = h*fp(v + k2/2.)
        k4 = h*fp(v + k3)
        rk = v + (k1 + 2.*k2 + 2.*k3 + k4) / 6.
    end function rk

    function rkf(v, h, fp)
        intent (in) :: v, h
        procedure(derivada) :: fp
        real(8) v(0:), h
        real(8), dimension(0:size(v) - 1) :: k1, k2, k3, k4, k5, k6, rkf, e
            
        k1 = h*fp(v)
        k2 = h*fp(v + k1/4.0)
        k3 = h*fp(v + (3.0*k1 + 9.0*k2)/32.0)
        k4 = h*fp(v + (1932.0*k1 - 7200.0*k2 + 7296.0*k3)/2197.0)
        k5 = h*fp(v + 439.0*k1/216.0 - 8.0*k2 + 3680.0*k3/513.0 - 845.0*k4/4104.0)
        k6 = h*fp(v - 8.0*k1/27.0 + 2.0*k2 - 3544.0*k3/2565.0 + 1859.0*k4/4104.0 - 11.0*k5/40.0)
        rkf = v + (25.0*k1/216.0 + 1408.0*k3/2565.0 + 2197.0*k4/4104.0 - k5/5.0)
        e = k1/360.0 - 128.0*k3/4275.0 - 2197.0*k4/75240.0 + k5/50.0 + 2.0*k6/55.0
    end function rkf

    subroutine iterarVeces(met, fp, vinicial, h, rep)
        intent(in) :: vinicial, h, rep
        procedure(metodo) :: met
        procedure(derivada) :: fp
        integer(4) rep, i
        real(8) h, vinicial(0:), v(0:size(vinicial - 1))

        open(2, file='datos.dat')
        v = vinicial
        write(2, *) v
        do i = 1, rep
            v = met(v, h, fp)
            write(2, *) v
        end do
        close(2)
        call system('gnuplot -persist valores.p')
    end subroutine iterarVeces
    
    subroutine iterarValorFinal(met, fp, vinicial, h, xf)
        intent(in) :: vinicial, h, xf
        procedure(metodo) :: met
        procedure(derivada) :: fp
        real(8) h, xf, vinicial(0:), v(0:size(vinicial - 1))

        open(2, file='datos.dat')
        v = vinicial
        write(2, *) v
        do while (v(0) <= xf)
            v = met(v, h, fp)
            write(2, *) v
        end do
        close(2)
        call system('gnuplot -persist valores.p')
    end subroutine iterarValorFinal

end module EDOs

program principal
    use EDOs

    implicit none

    integer(4), parameter :: cant_ec = 1

    integer(4) repeticiones
    real(8) v(0:cant_ec), h, xf

    repeticiones = 15
    h = 1.
    v(0) = -5.
    v(1) = 25.
    xf = 25.
    call iterar(eulerSimple, fp, v, h, repeticiones)
    call iterar(eulerModificado, fp, v, h, repeticiones)
    call iterar(eulerMejorado, fp, v, h, repeticiones)
    call iterar(rk, fp, v, h, repeticiones)
    call iterar(rkf, fp, v, h, xf)

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

end program principal