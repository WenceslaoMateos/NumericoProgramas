module EDDP
    use SELs
    use arreglos
    implicit none

    integer(4), parameter :: DIRICHLET = 0, NEUMANN = 1

    type frontera
        integer(4) tipo
        real(8) valor
    end type frontera

    type frontera2
        integer(4) tipo
        procedure(funcionContorno), pointer, nopass :: valor
    end type frontera2

    type ProblemaParabolicas
        real(8) t0, tf, dt, x0, xf, dx
        integer(4) particionx, particiont
        real(8), dimension(:), allocatable :: iniciales, u
        procedure(discretizacionParabolicas), pointer, nopass :: calculoInterno, calculoIzquierda, calculoDerecha
    end type ProblemaParabolicas

    abstract interface
        function discretizacionParabolicas(pp, i)
            import ProblemaParabolicas
            type(ProblemaParabolicas), intent(in) :: pp
            integer(4), intent(in) :: i
            real(8) discretizacionParabolicas
        end function
    end interface
    
    abstract interface
        function funcionContorno(x, t, algo)
            real(8), intent(in) :: x, t, algo
            real(8) funcionContorno
        end function funcionContorno
    end interface

    abstract interface
        function poisson(x, y)
            real(8), intent(in) :: x, y
            real(8) poisson
        end function poisson
    end interface

    abstract interface
        function funr(dx, dt)
            real(8), intent(in) :: dx, dt
            real(8) funr
        end function funr
    end interface

contains

!--------------------------------ELIPTICAS--------------------------------!

    function laplace(x, y)
        real(8), intent(in) :: x, y
        real(8) laplace

        laplace = 0.
    end function laplace

    function elipticas(x0, xf, y0, yf, nx, my, superior, inferior, izquierda, derecha, si, sd, ii, id, f, tol)
        integer(4), intent(in) :: nx, my
        real(8), intent(in) :: x0, xf, y0, yf, tol
        type(frontera), intent(in) :: si, sd, ii, id
        type(frontera), dimension(1:nx-1), intent(in) :: superior, inferior
        type(frontera), dimension(1:my-1), intent(in) :: izquierda, derecha
        procedure(poisson) :: f
        real(8), dimension(0:my, 0:nx) :: elipticas, term_ind, d, ud, bd, ld, rd, xini
        integer(4) desde, hasta, i, j, offset
        real(8) h, k, h2, k2, x, y

        h = (xf - x0) / nx
        k = (yf - y0) / my
        h2 = h**2.
        k2 = k**2.

        ! Nodos internos
        y = yf
        do j = 1, my - 1
            y = y - k
            x = x0
            do i = 1, nx - 1
                x = x + h
                term_ind(j, i) = f(x, y) * h2 * k2
                d(j, i) = -2. * (h2 + k2)
                ud(j, i) = h2
                bd(j, i) = h2
                ld(j, i) = k2
                rd(j, i) = k2
            end do
        end do

        ! Esquinas
        ud(0, 0) = 0.
        ld(0, 0) = 0.
        if (si%tipo == DIRICHLET) then
            term_ind(0, 0) = si%valor
            d(0, 0) = 1.
            bd(0, 0) = 0.
            rd(0, 0) = 0.
        elseif (si%tipo == NEUMANN) then
            term_ind(0, 0) = f(x0, yf) * h2 * k2 + 2.*k2*h * si%valor - 2.*h2*k * si%valor
            d(0, 0) = -2. * (h2 + k2)
            bd(0, 0) = 2.*h2
            rd(0, 0) = 2.*k2
        end if

        ud(0, nx) = 0.
        rd(0, nx) = 0.
        if (sd%tipo == DIRICHLET) then
            term_ind(0, nx) = sd%valor
            d(0, nx) = 1.
            bd(0, nx) = 0.
            ld(0, nx) = 0.
        elseif (sd%tipo == NEUMANN) then
            term_ind(0, nx) = f(xf, yf) * h2 * k2 - 2.*k2*h * sd%valor - 2.*h2*k * sd%valor
            d(0, nx) = -2. * (h2 + k2)
            bd(0, nx) = 2.*h2
            ld(0, nx) = 2.*k2
        end if

        bd(my, 0) = 0.
        ld(my, 0) = 0.
        if (ii%tipo == DIRICHLET) then
            term_ind(my, 0) = ii%valor
            d(my, 0) = 1.
            ud(my, 0) = 0.
            rd(my, 0) = 0.
        elseif (ii%tipo == NEUMANN) then
            term_ind(my, 0) = f(x0, y0) * h2 * k2 + 2.*k2*h * ii%valor + 2.*h2*k * ii%valor
            d(my, 0) = -2. * (h2 + k2)
            ud(my, 0) = 2.*h2
            rd(my, 0) = 2.*k2
        end if

        bd(my, nx) = 0.
        rd(my, nx) = 0.
        if (id%tipo == DIRICHLET) then
            term_ind(my, nx) = id%valor
            d(my, nx) = 1.
            ud(my, nx) = 0.
            ld(my, nx) = 0.
        elseif (id%tipo == NEUMANN) then
            term_ind(my, nx) = f(xf, y0) * h2 * k2 - 2.*k2*h * id%valor + 2.*h2*k * id%valor
            d(my, nx) = -2. * (h2 + k2)
            ud(my, nx) = 2.*h2
            ld(my, nx) = 2.*k2
        end if

        ! Bordes
        x = x0
        do i = 1, nx - 1
            x = x + h

            ud(0, i) = 0.
            if (superior(i)%tipo == DIRICHLET) then
                term_ind(0, i) = superior(i)%valor
                d(0, i) = 1.
                bd(0, i) = 0.
                ld(0, i) = 0.
                rd(0, i) = 0.
            elseif (superior(i)%tipo == NEUMANN) then
                term_ind(my, i) = f(x, yf) * h2 * k2 - 2.*h2*k * superior(i)%valor
                d(0, i) = -2. * (h2 + k2)
                bd(0, i) = 2.*h2
                ld(0, i) = k2
                rd(0, i) = k2
            end if

            bd(my, i) = 0.
            if (inferior(i)%tipo == DIRICHLET) then
                term_ind(0, i) = inferior(i)%valor
                d(my, i) = 1.
                ud(my, i) = 0.
                ld(my, i) = 0.
                rd(my, i) = 0.
            elseif (inferior(i)%tipo == NEUMANN) then
                term_ind(my, i) = f(x, y0) * h2 * k2 + 2.*h2*k * inferior(i)%valor
                d(my, i) = -2. * (h2 + k2)
                ud(my, i) = 2.*h2
                ld(my, i) = k2
                rd(my, i) = k2
            end if
        end do

        y = yf
        do j = 1, my - 1
            y = y - k

            ld(j, 0) = 0.
            if (izquierda(j)%tipo == DIRICHLET) then
                term_ind(j, 0) = izquierda(j)%valor
                d(j, 0) = 1.
                ud(j, 0) = 0.
                bd(j, 0) = 0.
                rd(j, 0) = 0.
            elseif (izquierda(j)%tipo == NEUMANN) then
                term_ind(j, 0) = f(x0, y) * h2 * k2 + 2.*k2*h * izquierda(j)%valor
                d(j, 0) = -2. * (h2 + k2)
                ud(j, 0) = h2
                bd(j, 0) = h2
                rd(j, 0) = 2.*k2
            end if

            rd(j, nx) = 0.
            if (derecha(j)%tipo == DIRICHLET) then
                term_ind(j, nx) = derecha(j)%valor
                d(j, nx) = 1.
                ud(j, nx) = 0.
                bd(j, nx) = 0.
                ld(j, nx) = 0.
            elseif (derecha(j)%tipo == NEUMANN) then
                term_ind(j, nx) = f(xf, y) * h2 * k2 - 2.*k2*h * derecha(j)%valor
                d(j, nx) = -2. * (h2 + k2)
                ud(j, nx) = h2
                bd(j, nx) = h2
                ld(j, nx) = 2.*k2
            end if
        end do

        xini = 0.
        elipticas = gaussSeidelMatricial(d, ud, bd, ld, rd, term_ind, xini, tol)
    end function elipticas

    function generarDistribucion(nx, my, x0, x1, y0, y1, resul, si, sd, ii, id, superior, inferior, izquierda, derecha)
        integer(4), intent(in) :: nx, my
        type(frontera), dimension(:), intent(in) :: superior, inferior, izquierda, derecha
        real(8), intent(in) :: resul(1:(nx - 1) * (my - 1)), si, sd, ii, id, x0, x1, y0, y1
        real(8) generarDistribucion(1:my+1, 1:nx+1), h, k
        integer(4) i, offset

        h = (x1 - x0) / nx
        k = (y1 - y0) / my

        ! Interno
        offset = 0
        do i = 2, my
            generarDistribucion(i, 2:nx) = resul(1+offset:nx-1)
            offset = offset + nx - 1
        end do

        ! Esquinas
        generarDistribucion(1, 1) = si
        generarDistribucion(1, nx+1) = sd
        generarDistribucion(my+1, 1) = ii
        generarDistribucion(my+1, nx+1) = id

        ! Bordes
        do i = 2, nx
            if (superior(i - 1)%tipo == DIRICHLET) then
                generarDistribucion(1, i) = superior(i - 1)%valor
            elseif (superior(i - 1)%tipo == NEUMANN) then
                generarDistribucion(1, i) = generarDistribucion(3, i) + 2. * k * superior(i - 1)%valor
            end if

            if (inferior(i - 1)%tipo == DIRICHLET) then
                generarDistribucion(my + 1, i) = inferior(i - 1)%valor
            elseif (inferior(i - 1)%tipo == NEUMANN) then
                generarDistribucion(my + 1, i) = generarDistribucion(my - 1, i) - 2. * k * inferior(i - 1)%valor
            end if
        end do

        do i = 2, my
            if (izquierda(i - 1)%tipo == DIRICHLET) then
                generarDistribucion(i, 1) = izquierda(i - 1)%valor
            elseif (izquierda(i - 1)%tipo == NEUMANN) then
                generarDistribucion(i, 1) = generarDistribucion(i, 3) - 2. * h * izquierda(i - 1)%valor
            end if

            if (derecha(i - 1)%tipo == DIRICHLET) then
                generarDistribucion(i, nx + 1) = derecha(i - 1)%valor
            else if (derecha(i - 1)%tipo == NEUMANN) then
                generarDistribucion(i, nx + 1) = generarDistribucion(i, nx - 1) + 2. * h * derecha(i - 1)%valor
            end if
        end do
    end function generarDistribucion

    subroutine grabarDatos(distribucion, x0, x1, y0, y1, nx, my, archivo)
        intent(in) :: distribucion, x0, x1, y0, y1, nx, my, archivo
        integer(4) nx, my, i, j
        real(8) distribucion(:, :), x0, x1, y0, y1, h, k, x, y
        character(len=*) archivo

        open(unit=2, file=archivo, access='SEQUENTIAL', status='REPLACE')
        h = (x1 - x0) / nx
        k = (y1 - y0) / my
        y = y1
        do i = 1, my + 1
            x = x0
            do j = 1, nx + 1
                write(2, *) x, y, distribucion(i, j)
                x = x + h
            end do
            write(2, *)
            y = y - k
        end do
        close(2, status='KEEP')
    end subroutine

!--------------------------------PARABOLICAS--------------------------------!

    function ejemploDiscretizacion(pp, i)
        type(ProblemaParabolicas), intent(in) :: pp
        integer(4), intent(in) :: i
        real(8) ejemploDiscretizacion
        real(8), parameter :: r = 0.5

        ejemploDiscretizacion = r * (pp%u(i + 1) + pp%u(i - 1)) + (1 - 2 * r) * pp%u(i)
    end function

    function formularProblemaParabolicas(x0, xf, particionx, t0, tf, particiont, &
            iniciales, calculoInterno, calculoIzquierda, calculoDerecha)
        real(8), intent(in) :: x0, xf, t0, tf, iniciales(:)
        integer(4), intent(in) :: particionx, particiont
        procedure(discretizacionParabolicas) calculoInterno, calculoIzquierda, calculoDerecha
        type(ProblemaParabolicas) pp, formularProblemaParabolicas

        pp%x0 = x0
        pp%xf = xf
        pp%particionx = particionx
        pp%dx = (xf - x0) / particionx

        pp%t0 = t0
        pp%tf = tf
        pp%particiont = particiont
        pp%dt = (tf - t0) / particiont

        pp%iniciales = iniciales

        pp%calculoInterno => calculoInterno
        pp%calculoIzquierda => calculoIzquierda
        pp%calculoDerecha => calculoDerecha

        formularProblemaParabolicas = pp
    end function

    subroutine explicito(pp, archivo)
        type(ProblemaParabolicas), intent(inout) :: pp
        character(len=*), intent(in) :: archivo
        real(8), dimension(pp%particionx + 2) :: x
        real(8), dimension(size(pp%iniciales)) :: usig
        real(8) t
        integer(4) n, i

        !escritura inicial en el archivo
        open(2, FILE=archivo)
        pp%dx = (pp%xf - pp%x0) / pp%particionx
        x(1) =  0.
        x(2) = pp%x0
        do i = 3, pp%particionx + 2
            x(i) = x(i - 1) + pp%dx
        end do
        write(2, *) x
        write(2, *)

        !core de metodo explicito
        n = size(pp%iniciales)
        pp%u = pp%iniciales
        t = pp%t0
        write(2, *) t, pp%u

        pp%dt = (pp%tf - pp%t0) / pp%particiont
        do while(t <= pp%tf)
            t = t + pp%dt
            usig(1) = pp%calculoIzquierda(pp, 1)
            do i = 2, n - 1
                usig(i) = pp%calculoInterno(pp, i)
            end do
            usig(ubound(usig, 1)) = pp%calculoDerecha(pp, ubound(usig, 1))
            pp%u = usig
            write(2, *) t, pp%u
        end do
        close(2)
    end subroutine explicito

    subroutine implicitoThomas(iniciales, ci, cd, x0, xf, t0, tf, erre, particionx, particiont, archivo)
        real(8), intent(in) :: t0, x0, xf, tf
        type(frontera), intent(in) :: ci, cd
        integer(4), intent(in) :: particionx, particiont
        real(8), dimension(particionx - 1), intent(in) :: iniciales
        real(8), dimension(particionx + 2) :: x
        character(len=*), intent(in) :: archivo
        procedure(funr) :: erre
        real(8), dimension(size(iniciales) + 2, 1) :: uant, u
        real(8), dimension(particionx - 1, 1) :: term_ind
        real(8), dimension(particionx - 1) :: u_o, d_o, l_o
        real(8) t, r, dt, dx
        integer(4) n, i
        
        !escitura inicial en el archivo
        open(2, FILE=archivo)
        dx = (xf - x0) / particionx
        x(1) =  0.
        x(2) = x0
        do i = 3, particionx + 2
            x(i) = x(i-1) + dx
        end do
        write(2, *) x
        write(2, *)

        n = size(iniciales) + 2
        u(1, 1) = ci%valor
        u(n, 1) = cd%valor
        u(2:n-1 ,1) = iniciales
        t = t0
        write(2, *) t, u(:, 1)
        
        dt = (tf - t0) / particiont
        r = erre(dx, dt)
        d_o = 2 + 2 * r
        u_o = -r
        u_o(particionx - 1) = 0.
        l_o = -r
        l_o(1) = 0.
        do while(t <= tf)
            t = t + dt
            uant = u
            term_ind = 0.
            term_ind(1, 1) = r * ci%valor 
            term_ind(n-2, 1) = r * cd%valor
            !terminos independientes
            do i = 1, n-2
                term_ind(i, 1) = term_ind(i, 1) + r * uant(i, 1) + (2 - 2 * r) * uant(i+1, 1) + r * uant(i+2, 1)
            end do
            u(2:n-1, :) = thomas(u_o, d_o, l_o, term_ind)
            write(2, *) t, u(:, 1)
        end do
        close(2)
    end subroutine implicitoThomas

    subroutine implicito(iniciales, ci, cd, x0, xf, t0, tf, erre, particionx, particiont, archivo)
        real(8), intent(in) :: t0, x0, xf, tf
        type(frontera), intent(in) :: ci, cd
        integer(4), intent(in) :: particionx, particiont
        real(8), dimension(particionx - 1), intent(in) :: iniciales
        real(8), dimension(particionx + 2) :: x
        character(len=*), intent(in) :: archivo
        procedure(funr) :: erre
        real(8), dimension(size(iniciales) + 2, 1) :: uant, u
        real(8), dimension(particionx - 1, particionx - 1) :: matriz
        real(8), dimension(particionx - 1, 1) :: term_ind
        real(8) t, r, dt, dx
        integer(4) n, i
        
        !escitura inicial en el archivo
        open(2, FILE=archivo)
        dx = (xf - x0) / particionx
        x(1) =  0.
        x(2) = x0
        do i = 3, particionx + 2
            x(i) = x(i-1) + dx
        end do
        write(2, *) x
        write(2, *)

        n = size(iniciales) + 2
        u(1, 1) = ci%valor
        u(n, 1) = cd%valor
        u(2:n-1 ,1) = iniciales
        t = t0
        write(2, *) t, u(:, 1)
        
        dt = (tf - t0) / particiont
        r = erre(dx, dt)
        matriz = 0.
        !banda central
        do i = 1, n-2
            matriz(i, i) = 2 + 2 * r
        end do
        ! banda derecha
        do i = 1, n-3
            matriz(i, i+1) = -r
        end do
        ! banda izquierda
        do i = 2, n-2
            matriz(i, i-1) = -r
        end do
        do while(t <= tf)
            t = t + dt
            uant = u
            term_ind = 0.
            term_ind(1, 1) = r * ci%valor
            term_ind(n-2, 1) = r * cd%valor
            !terminos independientes
            do i = 1, n-2
                term_ind(i, 1) = term_ind(i, 1) + r * uant(i, 1) + (2 - 2 * r) * uant(i+1, 1) + r * uant(i+2, 1)
            end do
            u(2:n-1, :) = gaussSeidel(matriz, term_ind, uant(2:n-1, :), 0.0001_8)
            write(2, *) t, u(:, 1)
        end do
        close(2)
    end subroutine implicito
    
!--------------------------------HIPERBOLICAS--------------------------------!
    !preguntar si x inicial es distinto de 0
    subroutine hiperbolicas(v, r, tf, dt, dx, velocidades, archivo)
        real(8), intent(in) :: tf, r, dx, dt
        character(len=*), intent(in) :: archivo
        real(8) v(:), x, t, longitud
        real(8), dimension(1:ubound(v, 1)) :: vAnt, vAct, velocidades
        integer(4)  n, i
        
        n = ubound(v, 1)
        open(2, FILE=archivo)
        
        !escritura de las posiciones de X
        vAnt(1) = 0
        do i = 2, n
            vAnt(i) = vAnt(i-1) + dx
        end do
        write(2, *) vAnt

        !escritura del vector inicial
        write(2, *) v

        longitud = (n-1)*dx
        vAct = v
        vAnt = v
        vAnt(1) = 0
        vAnt(n) = 0
        x = dx

        !calculo y escritura de la primera estimacion
        do i = 2, n - 1
            vAct(i) = r * (vAnt(i+1) + vAnt(i-1)) / 2. + velocidades(i) * dt + (1.-r) * vAnt(i)
        end do
        write(2, *) vAct

        t = 2*dt
        do while (t <= tf)
            do i = 2, n - 1
                v(i) = r * (vAct(i+1) + vAct(i-1)) - vAnt(i) + (2. - 2. * r) * vAct(i)
            end do
            t = t + dt
            write(2, *) v
            vAnt = vAct
            vAct = v
        end do  
    end subroutine hiperbolicas

end module EDDP