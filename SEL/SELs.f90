module SELs

    implicit none

contains

function matrizAmpliada(matriz, term_ind)
    real(8), dimension(:, :), intent(in) :: matriz, term_ind
    real(8) matrizAmpliada(size(matriz, dim=1), size(matriz, dim=2) + size(term_ind, dim=2))
    integer(4) col_mat, col_amp

    col_mat = size(matriz, dim=2)
    col_amp = size(matriz, dim=2) + size(term_ind, dim=2)

    matrizAmpliada(:, :col_mat) = matriz
    matrizAmpliada(:, col_mat + 1:col_amp) = term_ind                
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

function solucionGaussJordan(matriz, term_ind)
    real(8), dimension(:, :), intent(in) :: matriz, term_ind
    real(8) solucionGaussJordan(size(matriz, dim=1), size(matriz, dim=2) + size(term_ind, dim=2))
    integer(4) i, filas

    solucionGaussJordan = gaussJordan(matriz, term_ind)
    filas = size(matriz, dim=1)
    do i = 1, filas
        solucionGaussJordan(i, :) = solucionGaussJordan(i, :) / solucionGaussJordan(i, i)
        solucionGaussJordan(i, i) = 1.
    end do
end function solucionGaussJordan

function jacobi(matriz, term_ind, xini, tol)
    real(8), dimension(:, :), intent(in) :: matriz, term_ind, xini
    real(8), intent(in) :: tol
    real(8), dimension(size(term_ind, dim = 1), size(term_ind, dim= 2 )) :: jacobi, xant
    real(8), dimension(size(term_ind, dim = 1)) :: sum
    real(8) e1, e2
    integer(4) i, j, n

    xant = xini
    n = size(jacobi, dim = 1)
    e1 = tol + 1
    e2 = tol + 1
    do while((e1 > tol) .and. (e2 > tol))
        do i = 1, n
            sum = 0
            do j = 1, i - 1
                sum = sum + matriz(i, j) * xant(j, :)
            end do
            do j = i + 1, n
                sum = sum + matriz(i, j) * xant(j, :)
            end do
            jacobi(i, :) = (term_ind(i, :) - sum) / matriz(i, i)
        end do
        e1 = maxval(abs(jacobi-xant))
        e2 = normaMatriz(residuo(matriz, jacobi, term_ind))
    end do
end function jacobi

function gaussSeidel(matriz, term_ind, xini, tol)
    real(8), dimension(:, :), intent(in) :: matriz, term_ind, xini
    real(8), intent(in) :: tol
    real(8), dimension(size(term_ind, dim = 1), size(term_ind, dim= 2 )) :: gaussSeidel, xant
    real(8), dimension(size(term_ind, dim = 1)) :: sum
    real(8) e1, e2
    integer(4) i, j, n

    xant = xini
    n = size(gaussSeidel, dim = 1)
    e1 = tol + 1
    e2 = tol + 1
    do while((e1 > tol) .and. (e2 > tol))
        do i = 1, n
            sum = 0
            do j = 1, i - 1
                sum = sum + matriz(i, j) * gaussSeidel(j, :)
            end do
            do j = i + 1, n
                sum = sum + matriz(i, j) * xant(j, :)
            end do
            gaussSeidel(i, :) = (term_ind(i, :) - sum) / matriz(i, i)
        end do
        e1 = maxval(abs(gaussSeidel-xant))
        e2 = normaMatriz(residuo(matriz, gaussSeidel, term_ind))
    end do
end function gaussSeidel

!function reduccionCrout(matriz)
!    intent(in) :: matriz
!    real(8) matriz(:, :), reduccionCrout(size(matriz, dim=1), size(matriz, dim=2))
!    integer(4) orden, i, fila, k, col
!
!    reduccionCrout = matriz
!    reduccionCrout(1, 2:) = reduccionCrout(1, 2:) / reduccionCrout(1, 1)
!    orden = size(matriz, dim=1)
!    do i = 2, orden
!        Calculo de la columna i de L
!        do fila = i + 1, orden
!            do k = 1, i - 1
!                reduccionCrout(fila, i) = reduccionCrout(fila, i) - reduccionCrout(fila, k) * reduccionCrout(k, i)
!            end do
!        end do
!
!        Calculo de la fila i de U
!        do col = i + 1, orden
!            do k = 1, i - 1
!                reduccionCrout(i, col) = reduccionCrout(i, col) - reduccionCrout(i, k) * reduccionCrout(k, col)
!            end do
!        
!            reduccionCrout(i, col) = reduccionCrout(i, col) / reduccionCrout(i, i)
!        end do
!    end do
!end function reduccionCrout
!
!function solucionCrout(matriz, term_ind)
!    real(8), dimension(:, :), intent(in) :: matriz, term_ind
!    real(8), dimension(size(term_ind, dim=1), size(term_ind, dim=2)) :: solucionCrout, c
!    real(8) aux(size(matriz, dim=1), size(matriz, dim=2))
!    integer(4) i, orden, k
!
!    orden = size(matriz, dim=1)
!    aux = reduccionCrout(matriz)
!
!    ! Calcula c
!    c(1, :) = term_ind(1, :) / aux(1, 1)
!    do i = 2, orden
!        c(i, :) = term_ind(i, :)
!        do k = 1, i - 1
!            c(i, :) = c(i, :) - aux(i, k) * c(k, :)
!        end do
!        c(i, :) = c(i, :) / aux(i, i)
!    end do
!
!    ! Calcula la Solución
!    solucionCrout(orden, :) = c(orden, :)
!    do i = orden -1, 1, -1
!        solucionCrout(i, :) = c(i, :)
!        do k = i + 1, orden
!            solucionCrout(i, :) = solucionCrout(i, :) - aux(i, k) * solucionCrout(k, :)
!        end do
!    end do
!end function

function identidad(orden)
    integer(4), intent(in) :: orden
    integer(4) i
    real(8) identidad(orden, orden)

    identidad = 0.
    do i = 1, orden
        identidad(i, i) = 1.
    end do
end function identidad

function matrizInversa(matriz) ! Llamar solo con matrices cuadradas
    intent(in) :: matriz
    integer(4) orden
    real(8), allocatable :: aux(:, :)
    real(8) matriz(:, :), matrizInversa(size(matriz, dim=1), size(matriz, dim=2))

    orden = size(matriz, dim=1)
    allocate(aux(orden, orden*2))
    aux = solucionGaussJordan(matriz, identidad(orden))
    matrizInversa = aux(:, orden + 1:)
    deallocate(aux)
end function

function residuo(mat, sol, term_ind)
    real(8), intent(in) :: mat(:,:), sol(:,:), term_ind(:,:)
    real(8) residuo(size(mat, dim=1), size(mat, dim=2))
    
    residuo = matmul(mat,sol) - term_ind
end function residuo

function normaMatriz(matriz)
    intent(in) :: matriz
    real(8) matriz(:, :), normaMatriz

    normaMatriz = maxval(sum(abs(matriz), 2))
end function normaMatriz

function condicion(matriz)
    intent(in) :: matriz
    real(8) matriz(:, :), condicion

    condicion = normaMatriz(matriz) * normaMatriz(matrizInversa(matriz))
end function condicion

end module SELs

program principal
    use SELs
    use arreglos

    implicit none

    integer, parameter :: orden = 3
    real(8) matriz(orden, orden), inversa(orden, orden)

    matriz(1, 1) = 1.
    matriz(1, 2) = 0.
    matriz(1, 3) = 0.
    matriz(2, 1) = -1.
    matriz(2, 2) = 2.
    matriz(2, 3) = 3.
    matriz(3, 1) = 0.
    matriz(3, 2) = 1.
    matriz(3, 3) = 2.
    inversa = matrizInversa(matriz)
    call mostrarMatriz(matriz)
    write(*, *)
    call mostrarMatriz(inversa)
    write(*, *)
    call mostrarMatriz(matrizInversa(inversa))
    write(*, *)
    write(*, *) "Condicion = ", condicion(matriz)
contains

end program principal