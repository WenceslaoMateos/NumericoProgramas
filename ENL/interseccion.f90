program interseccion
    use ENL

    implicit none

    call plotVarias([character(len=10) :: "-3*x+1", "x**2"], "interseccion", "[-5:5]")
    write(*, *) "Interseccion 1: ", biseccion(-4._8, -2._8, fx, 0.0001_8)
    write(*, *) "Interseccion 2: ", biseccion(0._8, 1._8, fx, 0.0001_8)

contains

function fx(x)
    real(8), intent(in) :: x
    real(8) fx

    fx = -3.*x + 1 - x**2.
end function
    
end program interseccion
