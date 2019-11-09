program ejercicio5
    use EDDP
    use arreglos
    use SELs

    implicit none

    real(8) t0, x0, xf, tf
    integer(4) particionx, particiont
    real(8), dimension(:), allocatable :: ini
    character(len=20) archivo
    type(ProblemaParabolicas) pp

    ! Adaptar segun problema
    x0 = 0.
    xf = 1.
    t0 = 0.
    tf = 600.
    particionx = 20
    particiont = 24000
    allocate(ini(0:particionx))
    call cargarIniciales(ini, x0, xf)
    archivo = 'resultados.dat'

    pp = formularProblemaParabolicas(x0, xf, particionx, t0, tf, particiont, &
        ini, internos, cIzquierda, cDerecha)

    call explicito(pp, archivo)
    
contains

subroutine cargarIniciales(ini, x0, xf)
    real(8), dimension(0:), intent(out) :: ini
    real(8), intent(in) :: x0, xf
    real(8) x, dx
    integer(4) i
    real(8), parameter :: PI = 3.141592653589793238462643

    dx = (xf - x0) / size(ini)
    ini(0) = 0.
    x = x0
    do i = 1, ubound(ini, 1) - 1
        x = x + dx
        ini(i) = sin(PI*x)
    end do
    ini(ubound(ini, 1)) = 0.
end subroutine cargarIniciales

! Adaptar segun discretizacion

function internos(pp, i)
    type(ProblemaParabolicas), intent(in) :: pp
    integer(4), intent(in) :: i
    real(8) internos
    real(8), parameter :: r = 0.5

    internos = r * (pp%u(i + 1) - 2*pp%u(i) + pp%u(i - 1)) + pp%u(i)
end function internos

function cIzquierda(pp, i)
    type(ProblemaParabolicas), intent(in) :: pp
    integer(4), intent(in) :: i
    real(8) cIzquierda

    cIzquierda = 0.
end function cIzquierda

function cDerecha(pp, i)
    type(ProblemaParabolicas), intent(in) :: pp
    integer(4), intent(in) :: i
    real(8) cDerecha

    cDerecha = 0.
end function cDerecha

end program ejercicio5