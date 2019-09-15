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
                gauss(i, j+1:) = gauss(i, j+1:) - gauss(j, j+1:) * gauss(i, j)  / gauss(j, j)
                gauss(i, j) = 0.0
            end do
        end do 
    end function gauss
end module SELs

program principal
    use SELs
    use arreglos

    implicit none

    integer, parameter :: orden = 4, cant_sis = 1
    real(8) matriz(orden, orden), term_ind(orden,cant_sis), aux(orden, orden+cant_sis)
    matriz(1,1) = 2.
    matriz(1,2) = 0.
    matriz(1,3) = 0.
    matriz(1,4) = 0.
    matriz(2,1) = 1.
    matriz(2,2) = 1.5
    matriz(2,3) = 0.
    matriz(2,4) = 0.
    matriz(3,1) = 0.
    matriz(3,2) = -3.
    matriz(3,3) = 0.5
    matriz(3,4) = 0.
    matriz(4,1) = 2.
    matriz(4,2) = -2.
    matriz(4,3) = 1.
    matriz(4,4) = 1.
    term_ind(1,1) = 3.
    term_ind(2,1) = 4.5
    term_ind(3,1) = -6.6
    term_ind(4,1) = 0.8
    aux = gauss(matriz, term_ind)
    call mostrarMatriz(aux, orden, orden+cant_sis)
    !el resultado es 1, 2 ,3

contains

end program principal