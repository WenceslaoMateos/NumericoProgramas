module interpolacion
    use SELs
    use arreglos
    implicit none
    
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
        real(8), dimension(size(x), 1) :: term_ind, xini
        integer(4) columna
        
        term_ind(:, 1) = y
        matriz(:, 1) = 1.
        do columna = 2, size(x)
            matriz(:, columna) = x ** (columna - 1)
        end do
        xini = 0.
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
    
end module interpolacion