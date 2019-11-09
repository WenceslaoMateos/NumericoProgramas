program principal
    use EDDP
    use arreglos
    use SELs

    implicit none
    
    real(8) x0, x1, y0, y1, tol
    integer(4) n, m
    real(8), dimension(:, :), allocatable :: distribucion
    type(frontera) si, sd, ii, id
    type(frontera), dimension(:), allocatable :: superior, inferior, izquierda, derecha

    x0 = 0.
    x1 = 3.
    y0 = 0.
    y1 = 3.
    n = 3
    m = 3
    tol = 1e-5

    allocate(superior(1:n-1), inferior(1:n-1))
    allocate(izquierda(1:m-1), derecha(1:m-1))
    
    superior(1)%valor = 2.3
    superior(2)%valor = 4.4
    inferior(1)%valor = 2.3
    inferior(2)%valor = 3.4
    izquierda%valor = 1.2
    derecha(1)%valor = 6.5
    derecha(2)%valor = 5.5

    superior%tipo = DIRICHLET
    inferior%tipo = DIRICHLET
    izquierda%tipo = DIRICHLET
    derecha%tipo = DIRICHLET

    si%valor = 1.2
    si%tipo = DIRICHLET
    sd%valor = 7.5
    sd%tipo = DIRICHLET
    ii%valor = 1.2
    ii%tipo = DIRICHLET
    id%valor = 4.5
    id%tipo = DIRICHLET

    distribucion = elipticas(x0, x1, y0, y1, n, m, &
        superior, inferior, izquierda, derecha, si, sd, ii, id, &
        laplace, tol)
    call grabarDatos(distribucion, x0, x1, y0, y1, n, m, 'valores.dat')
    call plot('valores.dat')

    deallocate(distribucion, superior, inferior, izquierda, derecha)

contains

    subroutine plot(archivo)
        character(len=*), intent(in) :: archivo

        open(unit=2, file="temporal.p", access='SEQUENTIAL', status='REPLACE')
        write(2, *) 'set title "Distribución de Temperaturas en una placa rectangular"'
        write(2, *) 'set xlabel "x"'
        write(2, *) 'set ylabel "y"'
        write(2, *) 'set pm3d map'
        write(2, *) 'set nokey'
        write(2, *) 'splot "', archivo, '"'
        call system('gnuplot -persist temporal.p')
        close(2, STATUS='DELETE')
    end subroutine plot

end program principal