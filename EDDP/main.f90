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
    real(8), dimension(m, n) :: distribucion
    integer(4) i, offset

    superior = 0.
    inferior = 0.
    izquierda = 0.
    derecha = 100.
    call generaMatriz(mat, term_ind, n-2, m-2, superior, inferior, izquierda, derecha)
    
    xini = 0.
    res = gaussSeidel(mat, term_ind, xini, 0.000000000001_8)
    call mostrarMatriz(res)

    ! Esquinas
    distribucion(1, 1) = 0.
    distribucion(1, n) = 100.
    distribucion(m, 1) = 0.
    distribucion(m, n) = 100.

    ! Bordes
    distribucion(1, 2:n-1) = superior
    distribucion(m, 2:n-1) = inferior
    distribucion(2:m-1, 1) = izquierda
    distribucion(2:m-1, n) = derecha
    offset = 0
    do i = 2, m - 1
        distribucion(i, 2:n-1) = res(1+offset:n-2, 1)
        offset = offset + n - 2
    end do

    call mostrarMatriz(distribucion, '(9F10.2)')
    call grabarMatriz(distribucion, 'valores.dat')
    call generarDatos(distribucion, 0._8, 20._8, 0._8, 10._8, n, m, 'valores.dat')
    call system('gnuplot -persist matriz.p')

contains

    subroutine generarDatos(distribucion, x0, x1, y0, y1, n, m, archivo)
        intent(in) :: distribucion, x0, x1, y0, y1, n, m, archivo
        integer(4) n, m, i, j
        real(8) distribucion(m, n), x0, x1, y0, y1, h, k, x, y
        character(len=*) archivo

        open(unit=2, file=archivo, access='SEQUENTIAL', status='REPLACE')
        h = (x1 - x0) / (n - 1)
        k = (y1 - y0) / (m - 1)
        y = y0
        do i = 1, m
            x = x0
            do j = 1, n
                write(2, *) x, y, distribucion(i, j)
                x = x + h
            end do
            write(2, *)
            y = y + k
        end do
        close(2, status='KEEP')
    end subroutine

end program principal