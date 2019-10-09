program principal
    use EDDP
    use arreglos
    use SELs

    implicit none
    
    real(8) x0, x1, y0, y1
    integer(4) n, m, orden
    real(8), dimension(:, :), allocatable :: mat, distribucion
    real(8), dimension(:), allocatable :: superior, inferior, izquierda, derecha
    real(8), dimension(:, :), allocatable :: term_ind, xini, res

    x0 = 0.
    x1 = 20.
    y0 = 0.
    y1 = 10.
    n = 40
    m = 20
    orden = (n - 1) * (m - 1)

    allocate(mat(1:orden, 1:orden), distribucion(1:m+1, 1:n+1))
    allocate(term_ind(1:orden, 1), xini(1:orden, 1), res(1:orden, 1))
    allocate(superior(1:n-1), inferior(1:n-1))
    allocate(izquierda(1:m-1), derecha(1:m-1))

    superior = 0.
    inferior = 0.
    izquierda = 0.
    derecha = 100.
    call generarSistema(mat, term_ind, x0, x1, y0, y1, n, m, superior, inferior, izquierda, derecha)
    ! call mostrarMatriz(mat)
    ! write (*, *)
    ! call mostrarMatriz(term_ind)
    ! write (*, *)
    
    xini = 0.
    res = gaussSeidel(mat, term_ind, xini, 0.000001_8)
    call mostrarMatriz(res)
    write (*, *)

    distribucion = generarDistribucion(n, m, res(:, 1), 0._8, 100._8, 0._8, 100._8, superior, inferior, izquierda, derecha)
    ! call mostrarMatriz(distribucion)
    call grabarDatos(distribucion, x0, x1, y0, y1, n, m, 'valores.dat')
    call system('gnuplot -persist matriz.p')

    deallocate(mat, distribucion, superior, inferior, izquierda, derecha, term_ind, xini, res)

end program principal