program principal
    use EDDP
    use arreglos
    use SELs

    implicit none

    real(8) t0, x0, xf, tf
    type(frontera) ci, cd
    integer(4) particionx, particiont
    real(8), dimension(:), allocatable :: iniciales
    character(len=20) archivo

    particionx = 5
    particiont = 4
    ci%tipo = DIRICHLET
    ci%valor = 0
    cd%tipo = DIRICHLET
    cd%valor = 10
    t0 = 0.
    x0 = 0.
    tf = 537.6
    xf = 20
    allocate(iniciales(particionx - 1))
    iniciales(:) = [2., 2., 2., 2.]
    archivo = 'resultados.dat'

    
    call implicito(iniciales, ci, cd, x0, xf, t0, tf, r, particionx, particiont, archivo)
    
    contains

    function r(dx, dt)
        real(8), intent(in) :: dx, dt
        real(8) r

        r = 1.
    end function r

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