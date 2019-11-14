program E14
    use interpolacion

    implicit none
    
    real(8), dimension(0:8) :: x, y

    x = [0.5, 2., 5., 7.5, 9., 18., 26., 28.5, 31.5]
    y = [0., 3., 4., 3., 2., 1.5, 1.66, 2.5, 0.]

    call grabarTraza(x, y, 0.01_8)
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