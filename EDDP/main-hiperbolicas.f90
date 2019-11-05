program principal
    use EDDP
    use arreglos
    use SELs

    implicit none
    
    real(8) dx, dt, tf, T, g, w
    real(8), dimension(1:9) :: vector, velocidades

    !Valores iniciales
    vector(:) = [0., 0.3, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.]
    T = 0.000179
    g = 9.8
    w = 40. 
    dx = 0.01
    dt = dx/sqrt(T*g/w)
    tf = dt * 200
    velocidades = 0.
    call hiperbolicas(vector, erre(dx, dt, T, g, w), tf, dt, dx, velocidades, "resultados.dat")

contains

    function erre(dx, dt, T, g, w)
        real(8), intent(in) :: dx, dt, T, g, w
        real(8) erre

        erre = T * g * dt**2/(w * dx**2)
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