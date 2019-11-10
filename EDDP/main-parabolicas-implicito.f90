program ejercicio
    use EDDP
    use arreglos
    use SELs

    implicit none
    
    real(8) t0, x0, xf, tf, tol
    integer(4) particionx, particiont
    real(8), dimension(:), allocatable :: ini
    character(len=20) archivo
    type(ParabolicasImplicito) pi

    allocate(ini(11))
    ini = 15.
    x0 = 0.
    xf = 0.1
    t0 = 0.
    tf = 3600.
    particionx = 10
    particiont = 500
    archivo = 'resultados.dat'

    pi = formularProblemaImplicito(x0, xf, particionx, t0, tf, particiont, &
        ini, internos, cIzquierda, cDerecha)

    tol = 0.0001
    call implicito(pi, archivo, tol)

contains

subroutine internos(pi, i, d, ld, rd, term_ind)
    type(ParabolicasImplicito), intent(inout) :: pi
    integer(4), intent(in) :: i
    real(8), dimension(:) :: d, ld, rd, term_ind
    real(8) r
    real(8), parameter :: alfa = 6.94e-6

    r = alfa * pi%dt / pi%dx**2.
    d(i) = 2. + 2.*r
    ld(i) = -r
    rd(i) = -r
    term_ind(i) = r * (pi%u(i+1) + pi%u(i-1)) + (2. - 2.*r) * pi%u(i)
end subroutine internos

subroutine cIzquierda(pi, i, d, ld, rd, term_ind)
    type(ParabolicasImplicito), intent(inout) :: pi
    integer(4), intent(in) :: i
    real(8), dimension(:) :: d, ld, rd, term_ind
    real(8) r
    real(8), parameter :: alfa = 6.94e-6

    r = alfa * pi%dt / pi%dx**2.
    d(i) = 2. + 2.*r
    ld(i) = 0.
    rd(i) = -2.*r
    term_ind(i) = 2.*r*pi%u(i+1) + (2. - 2.*r) * pi%u(i)
end subroutine cIzquierda

subroutine cDerecha(pi, i, d, ld, rd, term_ind)
    type(ParabolicasImplicito), intent(inout) :: pi
    integer(4), intent(in) :: i
    real(8), dimension(:) :: d, ld, rd, term_ind
    real(8) r
    real(8), parameter :: alfa = 6.94e-6, h = 186., k = 37.2, Tinf = 1400. 

    r = alfa * pi%dt / pi%dx**2.
    d(i) = 2. + 2.*r + 2.*r*h*pi%dx/k
    ld(i) = -2.*r
    rd(i) = 0.
    term_ind(i) = 2.*r*pi%u(i-1) + (2. - 2.*r - 2.*r*h*pi%dx/k) * pi%u(i) + 4.*r*h*pi%dx*Tinf/k
end subroutine cDerecha

end program ejercicio