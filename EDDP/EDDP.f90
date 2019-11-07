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

    function gaussSeidelElipticas(matriz, term_ind, xini, tol, n, m)
        real(8), dimension(:, :), intent(in) :: matriz, term_ind, xini
        real(8), intent(in) :: tol
        integer(4), intent(in) :: n, m
        real(8), dimension(size(term_ind, dim=1), size(term_ind, dim=2)) :: gaussSeidelElipticas, xant
        real(8) e1
        integer(4) i, orden, cont

        gaussSeidelElipticas = xini
        orden = size(gaussSeidelElipticas, dim=1)
        e1 = tol + 1
        cont = 0
        do while(e1 > tol)
            xant = gaussSeidelElipticas
            do i = 1, orden
                gaussSeidelElipticas(i, :) = term_ind(i, :)
                ! tiene izquierda
                if (mod(i, n) /= 1) then
                    gaussSeidelElipticas(i, :) = gaussSeidelElipticas(i, :) - matriz(i, i-1) * gaussSeidelElipticas(i-1, :)
                end if
                ! tiene derecha
                if (mod(i, n) /= 0) then
                    gaussSeidelElipticas(i, :) = gaussSeidelElipticas(i, :) - matriz(i, i+1) * gaussSeidelElipticas(i+1, :)
                end if
                ! tiene arriba
                if (i > n) then
                    gaussSeidelElipticas(i, :) = gaussSeidelElipticas(i, :) - matriz(i, i-n) * gaussSeidelElipticas(i-n, :)
                end if
                ! tiene abajo
                if (i / n + 1 < m) then
                    gaussSeidelElipticas(i, :) = gaussSeidelElipticas(i, :) - matriz(i, i+n) * gaussSeidelElipticas(i+n, :)
                end if
                gaussSeidelElipticas(i, :) = gaussSeidelElipticas(i, :) / matriz(i, i)
            end do
            e1 = errorRelativo(gaussSeidelElipticas, xant, mNormaM)
            cont = cont + 1
        end do
        write(*, *) "Iteraciones: ", cont
    end function gaussSeidelElipticas

    function laplace(x, y)
        real(8), intent(in) :: x, y
        real(8) laplace

        laplace = 0.
    end function laplace
    
    subroutine generarSistema(mat, term_ind, x0, x1, y0, y1, nx, my, superior, inferior, izquierda, derecha, f)
        integer(4), intent(in) :: nx, my
        real(8), intent(in) :: x0, x1, y0, y1
        real(8), dimension(:, :), intent(out) :: mat, term_ind
        type(frontera), dimension(:), intent(in) :: superior, inferior, izquierda, derecha
        procedure(poisson) :: f
        integer(4) desde, hasta, i, n, m, offset
        real(8) h, k, x, y

        h = (x1 - x0) / nx
        k = (y1 - y0) / my
        n = nx - 1
        m = my - 1
        mat = 0.

        ! Iniciacion de terminos independientes
        x = x0
        y = y1
        do i = 1, size(term_ind, dim=1)
            if (mod(i, n) == 1) then
                x = x0 + h
                y = y - k
            else
                x = x + h
            end if
            term_ind(i, :) = f(x, y) * h**2. * k**2
        end do

        ! ---------------GENERACION DE TERMINOS INDEPENDIENTES--------------- !
        ! Condiciones superiores
        desde = 1
        hasta = n
        do i = desde, hasta
            if (superior(i)%tipo == DIRICHLET) then
                term_ind(i, :) = term_ind(i, :) - superior(i)%valor * h**2.
            elseif (superior(i)%tipo == NEUMANN) then
                term_ind(i, :) = term_ind(i, :) - 2. * k * superior(i)%valor
                mat(i, i + n) = h**2.
            end if
        end do

        ! Condiciones inferiores
        desde = n * (m - 1) + 1
        hasta = n * m
        offset = 0
        do i = desde, hasta
            offset = offset + 1
            if (inferior(offset)%tipo == DIRICHLET) then
                term_ind(i, :) = term_ind(i, :) - inferior(offset)%valor * h**2.
            elseif (inferior(offset)%tipo == NEUMANN) then
                term_ind(i, :) = term_ind(i, :) + 2. * k * inferior(offset)%valor
                mat(i, i - n) = h**2.
            end if
        end do

        ! Condiciones izquierda
        desde = 1
        hasta = n * (m - 1) + 1
        offset = 0
        do i = desde, hasta, n
            offset = offset + 1
            if (izquierda(offset)%tipo == DIRICHLET) then
                term_ind(i, :) = term_ind(i, :) - izquierda(offset)%valor * k**2.
            elseif (izquierda(offset)%tipo == NEUMANN) then
                term_ind(i, :) = term_ind(i, :) + 2. * h * izquierda(offset)%valor
                mat(i, i + 1) = k**2.
            end if
        end do
        
        ! Condiciones derecha
        desde = n
        hasta = n * m
        offset = 0
        do i = desde, hasta, n
            offset = offset + 1
            if (derecha(offset)%tipo == DIRICHLET) then
                term_ind(i, :) = term_ind(i, :) - derecha(offset)%valor * k**2.
            elseif (derecha(offset)%tipo == NEUMANN) then
                term_ind(i, :) = term_ind(i, :) - 2. * h * derecha(offset)%valor
                mat(i, i - 1) = k**2.
            end if
        end do

        ! ----------------------GENERACION DE DIAGONAL---------------------- !
        do i = 1, n * m
            mat(i, i) = -2. * (h**2. + k**2.)
        end do

        ! ------------------------GENERACION DE UNOS------------------------ !
        ! Banda superiores
        do i = n + 1, n * m
            mat(i, i - n) = mat(i, i - n) + h**2.
        end do

        ! Banda inferiores
        do i = 1, n * (m - 1)
            mat(i, i + n) = mat(i, i + n) + h**2.
        end do

        ! Banda izquierda
        do i = 2, n * m
            if (mod(i, n) /= 1) then
                mat(i, i - 1) = mat(i, i - 1) + k**2.
            end if
        end do

        ! Banda derecha
        do i = 1, n * m - 1
            if (mod(i, n) /= 0) then
                mat(i, i + 1) = mat(i, i + 1) + k**2.
            end if
        end do
    end subroutine generarSistema

    function elipticas(x0, x1, y0, y1, nx, my, superior, inferior, izquierda, derecha, f, tol)
        integer(4), intent(in) :: nx, my
        real(8), intent(in) :: x0, x1, y0, y1, tol
        real(8), dimension((nx - 1) * (my - 1)) :: elipticas
        type(frontera), dimension(:), intent(in) :: superior, inferior, izquierda, derecha
        procedure(poisson) :: f
        real(8), dimension(:), allocatable :: term_ind, d, ud, bd, ld, rd, xini
        integer(4) desde, hasta, i, n, m, offset, orden
        real(8) h, k, x, y

        h = (x1 - x0) / nx
        k = (y1 - y0) / my
        n = nx - 1
        m = my - 1
        orden = n * m
        allocate(term_ind(orden), d(orden), ud(orden), bd(orden), ld(orden), rd(orden), xini(orden))
        ud = 0.
        bd = 0.
        ld = 0.
        rd = 0.

        ! Iniciacion de terminos independientes
        x = x0
        y = y1
        do i = 1, size(term_ind, dim=1)
            if (mod(i, n) == 1) then
                x = x0 + h
                y = y - k
            else
                x = x + h
            end if
            term_ind(i) = f(x, y) * h**2. * k**2
        end do

        ! ---------------GENERACION DE TERMINOS INDEPENDIENTES--------------- !
        ! Condiciones superiores
        desde = 1
        hasta = n
        do i = desde, hasta
            if (superior(i)%tipo == DIRICHLET) then
                term_ind(i) = term_ind(i) - superior(i)%valor * h**2.
            elseif (superior(i)%tipo == NEUMANN) then
                term_ind(i) = term_ind(i) - 2. * k * superior(i)%valor
                ud(i) = h**2.
            end if
        end do

        ! Condiciones inferiores
        desde = n * (m - 1) + 1
        hasta = n * m
        offset = 0
        do i = desde, hasta
            offset = offset + 1
            if (inferior(offset)%tipo == DIRICHLET) then
                term_ind(i) = term_ind(i) - inferior(offset)%valor * h**2.
            elseif (inferior(offset)%tipo == NEUMANN) then
                term_ind(i) = term_ind(i) + 2. * k * inferior(offset)%valor
                bd(i) = h**2.
            end if
        end do

        ! Condiciones izquierda
        desde = 1
        hasta = n * (m - 1) + 1
        offset = 0
        do i = desde, hasta, n
            offset = offset + 1
            if (izquierda(offset)%tipo == DIRICHLET) then
                term_ind(i) = term_ind(i) - izquierda(offset)%valor * k**2.
            elseif (izquierda(offset)%tipo == NEUMANN) then
                term_ind(i) = term_ind(i) + 2. * h * izquierda(offset)%valor
                ld(i) = k**2.
            end if
        end do
        
        ! Condiciones derecha
        desde = n
        hasta = n * m
        offset = 0
        do i = desde, hasta, n
            offset = offset + 1
            if (derecha(offset)%tipo == DIRICHLET) then
                term_ind(i) = term_ind(i) - derecha(offset)%valor * k**2.
            elseif (derecha(offset)%tipo == NEUMANN) then
                term_ind(i) = term_ind(i) - 2. * h * derecha(offset)%valor
                rd(i) = k**2.
            end if
        end do

        ! ----------------------GENERACION DE DIAGONAL---------------------- !
        d = -2. * (h**2. + k**2.)

        ! Banda superiores
        ud(2:orden) = ud(2:orden) + h**2.

        ! Banda inferiores
        bd(1:n*(m-1)) = bd(1:n*(m-1)) + h**2.

        ! Banda izquierda
        do i = 2,orden
            if (mod(i, n) /= 1) then
                ld(i) = ld(i) + k**2.
            end if
        end do

        ! Banda derecha
        do i = 1, orden - 1
            if (mod(i, n) /= 0) then
                rd(i) = rd(i) + k**2.
            end if
        end do

        xini = 0.
        elipticas = gaussSeidel2D(d, ud, bd, ld, rd, term_ind, n, xini, tol)
        deallocate(term_ind, d, ud, bd, ld, rd, xini)
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