program E14
    use interpolacion

    implicit none
    
    real(8), dimension(0:17) :: x, y

    x = [1., 2., 5., 6., 7., 8., 10., 13., 17., 20., 23., 24., 25., 27., 27.7, 28., 29., 30.]
    y = [3., 3.7, 3.9, 4.2, 5.7, 6.6, 7.1, 6.7, 4.5, 7., 6.1, 5.6, 5.8, 5.2, 4.1, 4.30, 4.1, 3.]

    call grabarTraza(x(0:8), y(0:8), 0.01_8)
    call grabarTraza(x(8:14), y(8:14), 0.01_8)
    call grabarTraza(x(14:), y(14:), 0.01_8)
    call graficarSplines(x, y, 'Splines')

contains

    subroutine grabarTraza(x, y, salto)
        real(8), dimension(0:), intent(in) :: x, y
        real(8), intent(in) :: salto
        real(8), dimension(0:ubound(x, 1)) :: a, b, c, d
        logical exist
        integer(4) n
        real(8) xdesco

        a = y
        call splinesCubicos(a, b, c, d, calculaH(x))

        inquire(file="spline.dat", exist=exist)
        if (exist) then
            open(12, file="spline.dat", status="old", position="append", action="write")
        else
            open(12, file="spline.dat", status="new", action="write")
        end if

        n = ubound(x, 1)
        xdesco = x(0)
        do while (xdesco <= x(n))
            write(12, *) xdesco, evaluarSplines(a, b, c, d, x, xdesco)
            xdesco = xdesco + salto
        end do

        close(12)
    end subroutine

    subroutine graficarSplines(x, y, nombre)
        real(8), dimension(0:), intent(in) :: x, y
        character(LEN=*), intent(in) :: nombre
        integer(4) i
        real(8) xdesco

        open(4, FILE="puntos.dat")
        do i = 0, ubound(x, 1)
            write(4, *) x(i), y(i)
        end do

        open(unit=2, file="temporal.p", access='SEQUENTIAL', status='REPLACE')
        write(2, *) "set autoscale"
        write(2, *) "unset log"
        write(2, *) "unset label"
        write(2, *) "set grid"
        write(2, *) "set xtic auto"
        write(2, *) "set ytic auto"
        write(2, *) "set title '", nombre,"'"
        write(2, *) "set xlabel 'x'"
        write(2, *) "set ylabel 'f(x)'"
        write(2, *) "plot 'spline.dat' using 1:2 title '", nombre ,"' with lines, \"
        write(2, *) "'puntos.dat' using 1:2 title 'puntos' with points"
        call system('gnuplot -persist temporal.p')

        close(2, STATUS='DELETE')
        close(4, STATUS='DELETE')
    end subroutine graficarSplines

end program E14