module SELs

    implicit none

contains

function matrizAmpliada(matriz, term_ind)
    real(8), dimension(:,:), intent(in) :: matriz, term_ind
    real(8) matrizAmpliada(size(matriz, dim=1), size(matriz, dim=2) + size(term_ind, dim=2))
    integer(4) col_mat, col_amp

    col_mat = size(matriz, dim=2)
    col_amp = size(matriz, dim=2) + size(term_ind, dim=2)

    matrizAmpliada(:, :col_mat) = matriz(:, :)
    matrizAmpliada(:, col_mat + 1:col_amp) = term_ind(:, :)                
end function matrizAmpliada

function gauss(matriz, term_ind)
    real(8), dimension(:, :), intent(in) :: matriz, term_ind
    real(8) gauss(size(matriz, dim=1), size(matriz, dim=2) + size(term_ind, dim=2))
    integer(4) i, j, filas, columnas
    
    gauss = matrizAmpliada(matriz, term_ind)
    filas = size(gauss, dim=1)
    columnas = size(gauss, dim=2)
    do j = 1, columnas - 1
        do i = j + 1, filas
            gauss(i, :) = gauss(i, :) - gauss(j, :) * gauss(i, j)  / gauss(j, j)
            gauss(i, j) = 0.0
        end do
    end do 
end function gauss

function gaussJordan(matriz, term_ind)
    real(8), dimension(:, :), intent(in) :: matriz, term_ind
    real(8) gaussJordan(size(matriz, dim=1), size(matriz, dim=2) + size(term_ind, dim=2))
    integer(4) i, j, filas, columnas
    
    gaussJordan = gauss(matriz, term_ind)

    filas = size(matriz, dim=1)
    columnas = size(matriz, dim=2)
    do j = 2, columnas
        do i = 1, j - 1
            gaussJordan(i, :) = gaussJordan(i, :) - gaussJordan(j, :) * gaussJordan(i, j)  / gaussJordan(j, j)
            gaussJordan(i, j) = 0.0
        end do
    end do 
        
end function gaussJordan

end module SELs

program principal
    use SELs
    use arreglos

    implicit none

    integer, parameter :: orden = 4, cant_sis = 1
    real(8) matriz(orden, orden), term_ind(orden,cant_sis), aux(orden, orden+cant_sis)


    matriz(1,1) = 2.1756
    matriz(1,2) = 4.0231
    matriz(1,3) = -2.1732
    matriz(1,4) = 5.1967
    matriz(2,1) = -4.0231
    matriz(2,2) = 6.0000
    matriz(2,3) = 0.
    matriz(2,4) = 1.1973
    matriz(3,1) = -1.0000
    matriz(3,2) = -5.2107
    matriz(3,3) = 1.1111
    matriz(3,4) = 0.
    matriz(4,1) = 6.0235
    matriz(4,2) = 7.0000
    matriz(4,3) = 0.
    matriz(4,4) = -4.1561
    term_ind(1,1) = 17.102
    term_ind(2,1) = -6.1593
    term_ind(3,1) = 3.0004
    term_ind(4,1) = 0.0000
    aux = gauss(matriz, term_ind)
    call mostrarMatriz(aux, orden, orden+cant_sis)
    !el resultado es 1, 2 ,3

contains

end program principal