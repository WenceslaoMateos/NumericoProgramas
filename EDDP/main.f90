program principal
    use EDDP
    use arreglos
    use SELs

    implicit none
    
    integer(4), parameter :: n = 9, m = 5, orden = (n-2)*(m-2)
    real(8), dimension(1:orden, 1:orden) :: mat
    real(8), dimension(1:n-2) :: superior, inferior
    real(8), dimension(1:m-2) :: izquierda, derecha
    real(8), dimension(1:orden, 1) :: term_ind, xini, res
    real(8), dimension(n+1, m+1) :: distribucion
    integer(4) i, j

    superior = 0.
    inferior = 0.
    izquierda = 0.
    derecha = 100.
    call generaMatriz(mat, term_ind, n-2, m-2, superior, inferior, izquierda, derecha)
    
    xini = 0.
    res = gaussSeidel(mat, term_ind, xini, 0.000000000001_8)
    call mostrarMatriz(res)
    distribucion(1, 1) = 0.
    distribucion(1, m) = 100.
    distribucion(n, 1) = 0.
    distribucion(n, m) = 100.
    do i = 2, n-1
        distribucion(i, 1) = izquierda(i-1)
        distribucion(i, m) = derecha(i-1)
    end do
    do j = 2, m-1
        distribucion(1, j) = superior(j-1)
        distribucion(n, j) = inferior(j-1)
    end do
    do i = 1, n
        do j = 1, m
            distribucion(i+1, j+1) = res(i+j-1, 1)
        end do
    end do

    call grabarMatriz(distribucion,'valores.dat')
    call system('gnuplot -persist matriz.p')

end program principal