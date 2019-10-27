program main
    use interpolacion

    implicit none
    
    real(8), dimension(0:8) :: x, y, minCuad, lagrange, polinomioDescendente, polinomioAscendente, directo
    real(8), dimension(0:5) :: b, c, d, h
    real(8) xdesco
    
    xdesco = 5.35
    !x = [-4., -2., 0., 2.]
    !y = [-7.38, 0.52, 2., 14.52]
    !x = [0., 2., 4., 6.]
    !y = [0.25, 0.6, 0.9, 1.]
    !x = [5.1, 5.2, 5.3, 5.4, 5.5, 5.6]
    !y = [0.37798, 0.46852, 0.55437, 0.63469, 0.70867, 0.77557]
    
    !para splines
    !x = [27.70, 28., 29., 30.]
    !y = [4.10, 4.30, 4.10, 3.00]
    !x = [17.00, 20.00, 23.00, 24.00, 25.00, 27.00, 27.70]
    !y = [4.50, 7., 6.10, 5.60, 5.80, 5.20, 4.10]
    x = [1.00, 2.00, 5.00, 6.00, 7.00, 8.00, 10.00, 13.00, 17.00]
    y = [3.00, 3.70, 3.90, 4.20, 5.70, 6.60, 7.10, 6.70, 4.50]
    
    directo = polinomioAproximante(x, y)
    lagrange = polinomioLagrange(x, y)
    minCuad = mejorMinimosCuadrados(x, y)
    call diferenciasDivididas(x, y, polinomioDescendente, polinomioAscendente)
    !h = calculaH(x)
    !call splinesCubicos(y, b, c, d, h)
    write(*, *) 'Divididas Ascendente: ',evaluarPolinomio(polinomioAscendente, xdesco)
    write(*, *) 'Divididas Descendente: ',evaluarPolinomio(polinomioDescendente, xdesco)
    write(*, *) 'Equiespaciado: ',diferenciasEquiespaciado(x, y, xdesco)
    write(*, *) 'Lagrange: ',evaluarPolinomio(lagrange, xdesco)
    write(*, *) 'Directo: ',evaluarPolinomio(directo, xdesco)
    !write(*, *) 'Splines: ',evaluarSplines(y, b, c, d, x, xdesco)
    write(*, *) 'Minimos Cuadrados: ',evaluarPolinomio(minCuad, xdesco)
    call graficarPolinomio(minCuad, x, 0.01_8, 'archivo.dat')
    call graficarPolinomio(lagrange, x, 0.01_8, 'archivo.dat')

    !call graficarSplines(y, b, c, d, x, 0.01_8, 'archivo.dat')
    
contains

    function calculaH(x)
        real(8), dimension(0:) :: x
        real(8), dimension(0:ubound(x, 1)-1) :: calculaH
        integer(4) n

        n = ubound(x, 1)
        calculaH = x(1:n) - x(0:n-1)
    end function calculaH

    subroutine graficarSplines(a, b, c, d, x, salto, archivo)
        real(8), dimension(0:), intent(in) :: a, b, c, d, x
        real(8), intent(in) :: salto
        character(LEN=*), intent(in) :: archivo
        integer(4) n

        open(2, FILE=archivo)
        n = ubound(x, 1)
        xdesco = x(0)
        do while (xdesco<x(n))
            write(2, *) xdesco, evaluarSplines(a, b, c, d, x, xdesco)
            xdesco = xdesco + salto
        end do
        close(2)
        call plot(archivo)
    end subroutine graficarSplines

    subroutine graficarPolinomio(coeficientes, x, salto, archivo)
        real(8), dimension(0:), intent(in) :: coeficientes, x
        real(8), intent(in) :: salto
        character(LEN=*), intent(in) :: archivo
        integer(4) n

        open(2, FILE=archivo)
        n = ubound(x, 1)
        xdesco = x(0)
        do while (xdesco<x(n))
            write(2, *) xdesco, evaluarPolinomio(coeficientes, xdesco)
            xdesco = xdesco + salto
        end do
        close(2)
        call plot(archivo)
    end subroutine graficarPolinomio

    subroutine plot(archivo)
        intent(in) :: archivo
        character (LEN=*) :: archivo

        open(unit=2, file="temporal.p", access='SEQUENTIAL', status='REPLACE')
        write(2, *) "set autoscale"
        write(2, *) "unset log"
        write(2, *) "unset label"
        write(2, *) "set grid"
        write(2, *) "set xtic auto"
        write(2, *) "set ytic auto"
        write(2, *) "set title 'Spline'"
        write(2, *) "set xlabel 'x'"
        write(2, *) "set ylabel 'f(x)'"
        write(2, *) "plot '", archivo, "' using 1:2 title 'Spline' with lines"
        call system('gnuplot -persist temporal.p')
        close(2, STATUS='DELETE')
    end subroutine plot

end program main