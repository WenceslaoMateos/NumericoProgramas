program main
    use interpolacion

    implicit none
    
    real(8), dimension(0:3) :: x, y, lagrange, polinomioDescendente, polinomioAscendente, directo
    real(8) xdesco
    
    xdesco = 5.
    !x = [-4., -2., 0., 2.]
    !y = [-7.38, 0.52, 2., 14.52]
    x = [0., 2., 4., 6.]
    y = [0.25, 0.6, 0.9, 1.]
    !x = [5.1, 5.2, 5.3, 5.4, 5.5, 5.6]
    !y = [0.37798, 0.46852, 0.55437, 0.63469, 0.70867, 0.77557]
    directo = polinomioAproximante(x, y)
    lagrange = polinomioLagrange(x, y)
    call diferenciasDivididas(x, y, polinomioDescendente, polinomioAscendente)
    write(*, *) 'Divididas Ascendente: ',evaluarPolinomio(polinomioAscendente, xdesco)
    write(*, *) 'Divididas Descendente: ',evaluarPolinomio(polinomioDescendente, xdesco)
    write(*, *) 'Equiespaciado: ',diferenciasEquiespaciado(x, y, xdesco)
    write(*, *) 'Lagrange: ',evaluarPolinomio(lagrange, xdesco)
    write(*, *) 'Directo: ',evaluarPolinomio(directo, xdesco)

end program main