program principal
    use EDDP
    use arreglos
    use SELs

    implicit none
    
    real(8) x0, x1, y0, y1
    integer(4) i, offset, n, m, orden
    real(8), dimension(:, :), allocatable :: mat, distribucion
    real(8), dimension(:), allocatable :: superior, inferior, izquierda, derecha
    real(8), dimension(:, :), allocatable :: term_ind, xini, res

    x0 = 0.
    x1 = 20.
    y0 = 0.
    y1 = 10.
    n = 140
    m = 60
    orden = (n - 1) * (m - 1)

    allocate(mat(1:orden, 1:orden), distribucion(1:m+1, 1:n+1))
    allocate(term_ind(1:orden, 1), xini(1:orden, 1), res(1:orden, 1))
    allocate(superior(1:n+1), inferior(1:n+1))
    allocate(izquierda(1:m+1), derecha(1:m+1))
    superior = 0.
    inferior = 0.
    izquierda = 0.
    derecha = 100.

    call generaMatriz2(mat, term_ind, x0, x1, y0, y1, n, m, superior, inferior, izquierda, derecha)
    ! call mostrarMatriz(mat)
    ! write (*, *)
    ! call mostrarMatriz(term_ind)
    ! write (*, *)
    
    xini = 0.
    res = gaussSeidel(mat, term_ind, xini, 0.000001_8)
    call mostrarMatriz(res)
    write (*, *)

    ! Esquinas
    distribucion(1, 1) = 0.
    distribucion(1, n+1) = 100.
    distribucion(m+1, 1) = 0.
    distribucion(m+1, n+1) = 100.

    ! Bordes
    distribucion(1, 2:n) = superior
    distribucion(m+1, 2:n) = inferior
    distribucion(2:m, 1) = izquierda
    distribucion(2:m, n+1) = derecha
    offset = 0
    do i = 2, m
        distribucion(i, 2:n) = res(1+offset:n-1, 1)
        offset = offset + n - 1
    end do

    ! call mostrarMatriz(distribucion)
    call generarDatos(distribucion, 0._8, 20._8, 0._8, 10._8, n, m, 'valores.dat')
    call system('gnuplot -persist matriz.p')

    deallocate(mat, distribucion, superior, inferior, izquierda, derecha, term_ind, xini, res)
contains

    subroutine generarDatos(distribucion, x0, x1, y0, y1, n, m, archivo)
        intent(in) :: distribucion, x0, x1, y0, y1, n, m, archivo
        integer(4) n, m, i, j
        real(8) distribucion(:, :), x0, x1, y0, y1, h, k, x, y
        character(len=*) archivo

        open(unit=2, file=archivo, access='SEQUENTIAL', status='REPLACE')
        h = (x1 - x0) / n
        k = (y1 - y0) / m
        y = y0
        do i = 1, m + 1
            x = x0
            do j = 1, n + 1
                write(2, *) x, y, distribucion(i, j)
                x = x + h
            end do
            write(2, *)
            y = y + k
        end do
        close(2, status='KEEP')
    end subroutine

end program principal