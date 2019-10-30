program main
    use integrales
    implicit none

    real(8), dimension(0:5) :: x, y

    x = [0., 0.1, 0.2, 0.3, 0.4, 0.5]
    y = [1., 7., 4., 3., 5., 2.]

    write(*, *) 'Trapecios: ', trapecios(x, y)
    write(*, *) 'Simpson 1/3 solo al principio: ', simpsonUnTercio(x(0:4), y(0:4))
    write(*, *) 'Trapecios solo al final: ', trapecios(x(4:5), y(4:5))
    write(*, *) 'Simpson de 1/3 + Trapecios: ',trapecios(x(4:5), y(4:5)) + simpsonUnTercio(x(0:4), y(0:4))
    write(*, *) 'Simpson de 3/8 al principio: ',simpsonTresOctavos(x(0:3), y(0:3))
    write(*, *) 'Simpson de 1/3 al final: ',simpsonUnTercio(x(3:5), y(3:5))
    write(*, *) 'Simpson de 1/3 + Simpson de 3/8: ',simpsonTresOctavos(x(0:3), y(0:3)) + simpsonUnTercio(x(3:5), y(3:5))
    write(*, *) 'Romberg: ', romberg(x(0:4), y(0:4))
    write(*, *) 'Romberg + Trapecios: ', romberg(x(0:4), y(0:4)) + trapecios(x(4:5), y(4:5))
end program main