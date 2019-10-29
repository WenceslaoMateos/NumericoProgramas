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

    

end module integrales