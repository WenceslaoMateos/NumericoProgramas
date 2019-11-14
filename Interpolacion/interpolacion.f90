module interpolacion
    use SELs
    use arreglos
    implicit none

    abstract interface
        function ajusteMinCuad(x, y, coeficientes)
            real(8), intent(in), dimension(0:) :: x, y, coeficientes
            real(8) ajusteMinCuad
        end function ajusteMinCuad
    end interface

    abstract interface
        function funcionDerivadaNMas1(x)
            real(8), intent(in) :: x
            real(8) funcionDerivadaNMas1
        end function funcionDerivadaNMas1
    end interface

contains

    function evaluarPolinomio(coeficientes, x)
        real(8), dimension(0:), intent(in) :: coeficientes
        real(8), intent(in) :: x
        real(8) potencia, evaluarPolinomio
        integer(4) i
        
        potencia = 1.
        evaluarPolinomio = coeficientes(0)
        do i = 1, ubound(coeficientes, 1)
            potencia = potencia * x
            evaluarPolinomio = evaluarPolinomio + coeficientes(i) * potencia
        end do
    end function evaluarPolinomio

    function evaluarSplines(a, b, c, d, x, punto)
        real(8), dimension(0:), intent(in) :: a, b, c, d, x
        real(8) punto, evaluarSplines
        integer(4) i, n, j

        n = ubound(x, 1)
        i = 0
        do while((i <= n) .and. (punto >= x(i)))
            i = i + 1
        end do
        if (i <= n) then
            j = i - 1
            evaluarSplines = a(j) + b(j) * (punto - x(j)) + c(j) * (punto - x(j))**2 + d(j) * (punto - x(j))**3
        end if
    end function evaluarSplines
    
    function productoPolinomios(p1, p2)
        real(8), dimension(0:), intent(in) :: p1, p2
        real(8) productoPolinomios(0:ubound(p1, 1) + ubound(p2, 1))
        integer(4) i, j

        productoPolinomios = 0.
        do i = 0, ubound(p1, 1)
            do j = 0, ubound(p2, 1)
                productoPolinomios(i + j) = productoPolinomios(i + j) + p1(i) * p2(j)
            end do
        end do
    end function productoPolinomios

    function polinomioLagrange(x, y)
        real(8), dimension(0:), intent(in) :: x, y
        real(8), dimension(0:ubound(x, 1)) :: polinomioLagrange, numerador
        real(8) denominador
        integer(4) orden, k, i, limite

        orden = ubound(x, 1)
        polinomioLagrange = 0.
        do k = 0, orden
            denominador = 1.
            numerador(0) = 1.
            limite = 0
            do i = 0, k - 1
                denominador = denominador * (x(k) - x(i))
                limite = limite + 1
                numerador(:limite) = productoPolinomios(numerador(:limite - 1), [-x(i), 1._8])
            end do
            do i = k + 1, orden
                denominador = denominador * (x(k) - x(i))
                limite = limite + 1
                numerador(:limite) = productoPolinomios(numerador(:limite - 1), [-x(i), 1._8])
            end do
            polinomioLagrange = polinomioLagrange + numerador * y(k) / denominador
        end do
    end function polinomioLagrange
    
    function polinomioAproximante(x, y)
        real(8), dimension(0:), intent(in) :: x, y
        real(8), dimension(0:ubound(x, 1)) :: polinomioAproximante
        real(8), dimension(0:ubound(x, 1), 1) :: aux
        real(8), dimension(size(x), size(x)) :: matriz
        real(8), dimension(size(x), 1) :: term_ind
        integer(4) columna
        
        term_ind(:, 1) = y
        matriz(:, 1) = 1.
        do columna = 2, size(x)
            matriz(:, columna) = x ** (columna - 1)
        end do
        !Tener cuidado a la hora de ejecutar el metodo, a veces es necesario pivotear
        !Usamos un metodo directo por que la matriz casi nunca va a ser diagonalmente dominante
        aux = solucionGaussJordan(matriz, term_ind)
        polinomioAproximante = aux(:, 1)      
    end function polinomioAproximante

    subroutine diferenciasDivididas(x, y, polinomioDescendente, polinomioAscendente)
        real(8), dimension(0:), intent(in) :: x, y
        real(8), dimension(0:ubound(x, 1)), intent(out) :: polinomioDescendente, polinomioAscendente
        real(8), dimension(0:ubound(x, 1)) :: diferencias, diferenciasAuxD, diferenciasAuxA, auxAsc, auxDesc
        integer(4) orden, i, diferencia, cantDif, salto, factor

        orden = ubound(x, 1)
        diferencias = y
        diferenciasAuxA(0) = y(0)
        diferenciasAuxD(0) = y(orden)
        cantDif = orden
        salto = 0
        do i = 1, orden
            cantDif = cantDif - 1
            salto = salto + 1
            do diferencia = 0, cantDif
                diferencias(diferencia) = (diferencias(diferencia + 1) - diferencias(diferencia)) &
                                          / (x(diferencia + salto) - x(diferencia))
            end do
            diferenciasAuxA(i) = diferencias(0)
            diferenciasAuxD(i) = diferencias(cantDif)
        end do
        
        polinomioDescendente = 0.
        polinomioAscendente = 0.
        polinomioAscendente(0) = diferenciasAuxA(0)
        polinomioDescendente(0) = diferenciasAuxD(0)
        auxDesc(0) = 1.
        auxAsc(0) = 1.
        do i = 1, orden
            factor = orden - i + 1
            auxAsc(:i) = productoPolinomios(auxAsc(0:i - 1), [-x(i - 1), 1._8])
            auxDesc(:i) = productoPolinomios(auxDesc(0:i - 1), [-x(factor), 1._8])
            polinomioDescendente(:i) = polinomioDescendente(:i) + diferenciasAuxD(i) * auxDesc(:i)
            polinomioAscendente(:i) = polinomioAscendente(:i) + diferenciasAuxA(i) * auxAsc(:i)
        end do
    end subroutine diferenciasDivididas
    
    function diferenciasEquiespaciado(x, y, xdesco)
        real(8), dimension(0:), intent(in) :: x, y
        real(8), intent(in) :: xdesco
        real(8) :: diferenciasEquiespaciado

        if (xdesco < (x(0) + x(ubound(x, 1))) / 2.) then
            diferenciasEquiespaciado = equiespaciadoAscendente(x, y, xdesco)
        else
            diferenciasEquiespaciado = equiespaciadoDescendente(x, y, xdesco)
        end if
    end function diferenciasEquiespaciado

    function equiespaciadoDescendente(x, y, xdesco)
        real(8), dimension(0:), intent(in) :: x, y
        real(8), intent(in) :: xdesco
        real(8), dimension(0:ubound(x, 1)) :: diferencias
        real(8) :: equiespaciadoDescendente, s, fact, acum
        integer(4) orden, i, diferencia, cantDif
        
        orden = ubound(x, 1)
        cantDif = orden
        s = (xdesco - x(orden))/(x(1) - x(0))
        acum = 1.
        diferencias = y
        equiespaciadoDescendente = y(orden)
        fact = 1.
        do i = 1, orden
            cantDif = cantDif - 1
            do diferencia = 0, cantDif
                diferencias(diferencia) = (diferencias(diferencia + 1) - diferencias(diferencia))
            end do
            fact = fact * i
            acum = acum * (s + i - 1)
            equiespaciadoDescendente = equiespaciadoDescendente + (acum * diferencias(cantDif) / fact)
        end do
    end function equiespaciadoDescendente

    function equiespaciadoAscendente(x, y, xdesco)
        real(8), dimension(0:), intent(in) :: x, y
        real(8), intent(in) :: xdesco
        real(8), dimension(0:ubound(x, 1)) :: diferencias
        real(8) :: equiespaciadoAscendente, s, fact, acum
        integer(4) orden, i, diferencia, cantDif
        
        orden = ubound(x, 1)
        cantDif = orden
        s = (xdesco - x(0))/(x(1) - x(0))
        acum = 1.
        equiespaciadoAscendente = y(0)
        fact = 1.
        diferencias = y
        do i = 1, orden
            cantDif = cantDif - 1
            do diferencia = 0, cantDif
                diferencias(diferencia) = (diferencias(diferencia + 1) - diferencias(diferencia))
            end do
            fact = fact * i
            acum = acum * (s - i + 1)
            equiespaciadoAscendente = equiespaciadoAscendente + (acum * diferencias(0) / fact)
        end do
    end function equiespaciadoAscendente

    function calculaH(x)
        real(8), dimension(0:) :: x
        real(8), dimension(0:ubound(x, 1)-1) :: calculaH
        integer(4) n

        n = ubound(x, 1)
        calculaH = x(1:n) - x(0:n-1)
    end function calculaH

    subroutine splinesCubicos(a, b, c, d, h)
        real(8), dimension(0:), intent(in) :: a, h
        real(8), dimension(0:), intent(out) :: b, c, d
        real(8), dimension(:), allocatable :: u_o, d_o, l_o
        !real(8), dimension(0:ubound(a, 1), 1) :: term_ind
        integer i, n, INFO
        
        n = ubound(a, 1)
        allocate(u_o(0:n - 1), d_o(0:n), l_o(0:n - 1))
        u_o(0) = 0
        d_o(0) = 1
        l_o(0) = h(0)
        c(0) = 0
        do i = 1, n - 1
            u_o(i) = h(i)
            d_o(i) = 2 * (h(i - 1) + h(i))
            l_o(i) = h(i)
            c(i) = 3 * ((a(i + 1) - a(i)) / h(i) - (a(i) - a(i - 1)) / h(i - 1))
            !term_ind(i, 1) = 3 * ((a(i + 1) - a(i)) / h(i) - (a(i) - a(i - 1)) / h(i - 1))
        end do
        d_o(n) = 1
        l_o(n - 1) = 0
        c(n) = 0
        !term_ind(n, 1) = 0
        !term_ind = thomas(u_o, d_o, l_o, term_ind)
        !c = term_ind(:, 1)
        call DGTSV(n, 1, l_o, d_o, u_o, c, n, INFO)

        do i = 0, n - 1
            b(i) = (a(i + 1) - a(i)) / h(i) - h(i) * (2 * c(i) + c(i + 1)) / 3
            d(i) = (c(i + 1) - c(i)) / (3 * h(i))
        end do
        
        deallocate(u_o, d_o, l_o)
    end subroutine splinesCubicos

    function minimosCuadrados(x, y, grado)
        ! Calculo de coef
        real(8), intent(in), dimension(0:) :: x, y
        integer(4), intent(in) :: grado
        real(8), dimension(0:grado) :: minimosCuadrados
        real(8), dimension(0:grado, 0:grado) :: matriz
        real(8), dimension(0:grado, 1) :: term_ind, aux
        integer(4) i
        
        matriz(0, 0) = size(x)
        term_ind(0, 1) = sum(y)
        do i = 1, grado
            matriz(0, i) = sum(x**i)
            term_ind(i, 1) = sum(y * x**i)
        end do
        do i = 1, grado
            matriz(i, 0:grado - 1) = matriz(i - 1, 1:grado)
            matriz(i, grado) = sum(x**(grado + i))
        end do
        aux = solucionGaussJordan(matriz, term_ind)
        minimosCuadrados = aux(:, 1)
    end function minimosCuadrados

    subroutine mejorMinimosCuadrados(x, y, coeficientes, criterio)
        real(8), intent(in), dimension(0:) :: x, y
        real(8), dimension(:), allocatable, intent(out) :: coeficientes
        procedure(ajusteMinCuad) :: criterio
        real(8), dimension(:), allocatable :: nuevoMinimos, mejorMinimos
        real(8) vMejor, vNuevo
        integer(4) i, noMejoro

        mejorMinimos = minimosCuadrados(x, y, 1)
        vMejor = criterio(x, y, mejorMinimos)
        nuevoMinimos = minimosCuadrados(x, y, 2)
        vNuevo = criterio(x, y, nuevoMinimos)
        if (vMejor > vNuevo) then
            noMejoro = 0
            vMejor = vNuevo
            mejorMinimos = nuevoMinimos
        else
            noMejoro = 1
        end if

        i = 3
        do while(i < ubound(x, 1) .and. noMejoro < 2)
            nuevoMinimos = minimosCuadrados(x, y, i)
            vNuevo = criterio(x, y, nuevoMinimos)
            if (vMejor > vNuevo) then
                noMejoro = 0
                vMejor = vNuevo
                mejorMinimos = nuevoMinimos
            else
                noMejoro = noMejoro + 1
            end if

            i = i + 1 
        end do
        allocate(coeficientes(0:ubound(mejorMinimos, 1) - 1))
        coeficientes = mejorMinimos
        deallocate(nuevoMinimos, mejorMinimos)
    end subroutine mejorMinimosCuadrados

    function varianza(x, y, coeficientes)
        real(8), intent(in), dimension(0:) :: x, y
        real(8), dimension(0:), intent(in) :: coeficientes
        real(8) varianza
        integer(4) i, n, M
        
        M = size(x)
        n = ubound(coeficientes, 1)
        varianza = 0.
        do i = 0, M - 1
            varianza = varianza + (evaluarPolinomio(coeficientes, x(i)) - y(i)) ** 2
        end do
        varianza = varianza / (M - n - 1)
    end function varianza

    function RMS(x, y, coeficientes)
        real(8), intent(in), dimension(0:) :: x, y, coeficientes
        real(8) RMS
        integer(4) i, M
        
        M = size(x)
        RMS = 0.
        do i = 0, M - 1
            RMS = RMS + (evaluarPolinomio(coeficientes, x(i)) - y(i)) ** 2
        end do
        RMS = sqrt(RMS / M)
    end function RMS

    subroutine mostrarPolinomio(coeficientes)
        real(8), dimension(0:) :: coeficientes
        integer(4) orden, i
        
        orden = ubound(coeficientes, 1)
        
        write(*, '(A)', advance='NO') 'F(X)='
        do i = 0, orden - 1
            write(*, '(F12.8, A, I2, A)', advance='NO') coeficientes(i), ' X^', i, ' + '
        end do 
        write(*, '(F12.8, A, I2)') coeficientes(orden), ' X^', orden         
    end subroutine mostrarPolinomio

    !func es la derivada n+1 de la funcion real que estoy intentando aproximar
    function calculaError(punto, x, derivada)
        real(8), dimension(0:), intent(in) :: x
        real(8), intent(in) :: punto
        procedure(funcionDerivadaNMas1) :: derivada
        real(8) calculaError, m
        integer(4) n, i, fact

        n = ubound(x, 1)
        m = (x(n) - x(0))/2.
        calculaError = punto - x(0)
        fact = 1
        do i = 1, n
            calculaError = calculaError * (punto - x(i))
            fact = fact * i
        end do  
        calculaError = calculaError * derivada(m) / (fact * (n+1))
    end function calculaError
    
end module interpolacion