module integrales
    implicit none
    
contains
    
    function trapecios(x, y)
        real(8), dimension(0:), intent(in) :: x, y
        real(8) h, trapecios
        integer(4) i, n

        n = ubound(x, 1)
        !n es la cantidad de segmentos, nosotros tenemos n-1 puntos
        trapecios = 0.
        do i = 0, n-1
            h = x(i+1) - x(i)
            trapecios = trapecios + h * (y(i) + y(i + 1))
        end do
        trapecios = trapecios / 2.
    end function trapecios

    !ojo que solo sirve si los puntos son equiespacidos (de 3 en 3), no se toma con rangos diferentes
    function simpsonUnTercio(x, y)
        real(8), dimension(0:), intent(in) :: x, y
        real(8) h, simpsonUnTercio
        integer(4) n

        n = ubound(x, 1)
        h = x(1) - x(0)
        simpsonUnTercio = h * (y(0) + 4 * sum(y(1: n-1: 2)) + 2 * sum(y(2: n-2: 2)) + y(n)) / 3.
    end function simpsonUnTercio

    !ojo que solo sirve si los puntos son equiespacidos (de 3 en 3), no se toma con rangos diferentes
    function simpsonTresOctavos(x, y)
        real(8), dimension(0:), intent(in) :: x, y
        real(8) h, simpsonTresOctavos
        integer(4) i, n

        n = ubound(x, 1)
        simpsonTresOctavos = 0.
        do i = 0, n-3, 3
            simpsonTresOctavos = simpsonTresOctavos + (y(i) + 3 * y(i+1) + 3 * y(i+2) + y(i+3))
        end do
        h = x(1) - x(0)
        simpsonTresOctavos = 3 * h * simpsonTresOctavos / 8.
    end function simpsonTresOctavos

    ! Necesita 2**n + 1 puntos
    function romberg(x, y)
        real(8), dimension(0:), intent(in) :: x, y
        real(8) dx, romberg
        integer(4) paso, i, cant_trap, grado, reducciones
        real(8), dimension(:), allocatable :: trap

        dx = x(1) - x(0)
        paso = ubound(x, 1)
        cant_trap = log(real(paso)) / log(2.)
        allocate(trap(0:cant_trap))
    
        ! Genera los trapecios iniciales
        do i = 0, cant_trap
            trap(i) = trapecios(x(::paso), y(::paso))
            paso = paso / 2
        end do

        ! Extrapolacion de Richardson
        reducciones = cant_trap - 1
        do grado = 2, cant_trap + 1
            do i = 0, reducciones
                trap(i) = (4**(grado-1) * trap(i+1) - trap(i)) / (4**(grado-1) - 1)
            end do
            reducciones = reducciones - 1
        end do

        romberg = trap(0)
    end function romberg

end module integrales
