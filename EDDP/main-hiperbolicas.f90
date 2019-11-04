program principal
    use EDDP
    use arreglos
    use SELs

    implicit none
    
    real(8) dx, dt, tf
    real(8), dimension(1:9) :: vector

    vector(:) = [0., 0.3, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.]
    dx = 0.01
    dt = dx/sqrt(0.000179*9.8/40)
    tf = dt * 20
    call hiperbolicas(vector, erre(dx, dt), tf, dt, dx, "resultados.dat")

contains

    function erre(x, y)
        real(8), intent(in) :: x, y
        real(8) erre

        erre = x / y
    end function erre

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