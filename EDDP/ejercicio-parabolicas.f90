program ejercicio
    use EDDP
    use arreglos
    use SELs
    
    real(8) t0, x0, xf, tf
    integer(4) particionx, particiont
    real(8), dimension(:), allocatable :: ini
    character(len=20) archivo
    type(ProblemaParabolicas) pp

    allocate(ini(11))
    ini = 15.
    x0 = 0.
    xf = 0.1
    t0 = 0.
    tf = 3600.
    particionx = 10
    particiont = 500
    archivo = 'resultados.dat'

    pp = formularProblemaParabolicas(x0, xf, particionx, t0, tf, particiont, &
        ini, internos, cIzquierda, cDerecha)

    call explicito2(pp, archivo)

contains

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

    cIzquierda = pp%u(i + 1)
end function cIzquierda

function cDerecha(pp, i)
    type(ProblemaParabolicas), intent(in) :: pp
    integer(4), intent(in) :: i
    real(8) cDerecha
    real(8), parameter :: h = 186., k = 37.2, Tinf = 1400. 

    cDerecha = (pp%u(i - 1) + pp%dx * h * Tinf / k) / (1 + pp%dx * h / k)
end function cDerecha

end program ejercicio