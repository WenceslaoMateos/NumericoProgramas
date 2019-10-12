module EDDP
    implicit none

    integer(4), parameter :: DIRICHLET = 0, NEUMANN = 1

    type condicion
        integer(4) tipo
        real(8) valor
    end type condicion

    abstract interface
        function poisson(x, y)
            real(8), intent(in) :: x, y
            real(8) poisson
        end function poisson
    end interface
    
contains

    function laplace(x, y)
        real(8), intent(in) :: x, y
        real(8) laplace

        laplace = 0.
    end function laplace
    
    subroutine generarSistema(mat, term_ind, x0, x1, y0, y1, nx, my, superior, inferior, izquierda, derecha, f)
        integer(4), intent(in) :: nx, my
        real(8), intent(in) :: x0, x1, y0, y1
        real(8), dimension(:, :), intent(out) :: mat, term_ind
        real(8), dimension(:), intent(in) :: superior, inferior, izquierda, derecha
        procedure(poisson) :: f
        integer(4) desde, hasta, i, n, m
        real(8) h, k, x, y

        h = (x1 - x0) / nx
        k = (y1 - y0) / my
        n = nx - 1
        m = my - 1
        mat = 0.

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
        term_ind(desde:hasta, 1) = term_ind(desde:hasta, 1) - superior(:) * h**2.

        ! Condiciones inferiores
        desde = n * (m - 1) + 1
        hasta = n * m
        term_ind(desde:hasta, 1) = term_ind(desde:hasta, 1) - inferior(:) * h**2.

        ! Condiciones izquierda
        desde = 1
        hasta = n * (m - 1) + 1
        do i = desde, hasta, n
            term_ind(i, 1) = term_ind(i, 1) - izquierda(mod(i, n)+1) * k**2.
        end do
        
        ! Condiciones derecha
        desde = n
        hasta = n * m     
        do i = desde, hasta, n
            term_ind(i, 1) = term_ind(i, 1) - derecha(i/n) * k**2.
        end do

        ! ----------------------GENERACION DE DIAGONAL---------------------- !
        do i = 1, n * m
            mat(i, i) = -2. * (h**2. + k**2.)
        end do

        ! ------------------------GENERACION DE UNOS------------------------ !
        ! Banda superiores
        do i = n + 1, n * m
            mat(i, i - n) = h**2.
        end do

        ! Banda inferiores
        do i = 1, n * (m - 1)
            mat(i, i + n) = h**2.
        end do

        ! Banda izquierda
        do i = 2, n * m
            if (mod(i, n) /= 1) then
                mat(i, i - 1) = k**2.
            end if
        end do

        ! Banda derecha
        do i = 1, n * m - 1
            if (mod(i, n) /= 0) then
                mat(i, i + 1) = k**2.
            end if
        end do
    end subroutine generarSistema

    subroutine generarSistema2(mat, term_ind, x0, x1, y0, y1, nx, my, superior, inferior, izquierda, derecha, f)
        integer(4), intent(in) :: nx, my
        real(8), intent(in) :: x0, x1, y0, y1
        real(8), dimension(:, :), intent(out) :: mat, term_ind
        type(condicion), dimension(:), intent(in) :: superior, inferior, izquierda, derecha
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
    end subroutine generarSistema2

    function generarDistribucion(nx, my, resul, si, sd, ii, id, superior, inferior, izquierda, derecha)
        integer(4), intent(in) :: nx, my
        real(8), dimension(:), intent(in) :: superior, inferior, izquierda, derecha
        real(8), intent(in) :: resul(1:(nx - 1) * (my - 1)), si, sd, ii, id
        real(8) generarDistribucion(1:my+1, 1:nx+1)
        integer(4) i, offset

        ! Esquinas
        generarDistribucion(1, 1) = si
        generarDistribucion(1, nx+1) = sd
        generarDistribucion(my+1, 1) = ii
        generarDistribucion(my+1, nx+1) = id

        ! Bordes
        generarDistribucion(1, 2:nx) = superior
        generarDistribucion(my+1, 2:nx) = inferior
        generarDistribucion(2:my, 1) = izquierda
        generarDistribucion(2:my, nx+1) = derecha

        ! Interno
        offset = 0
        do i = 2, my
            generarDistribucion(i, 2:nx) = resul(1+offset:nx-1)
            offset = offset + nx - 1
        end do
    end function generarDistribucion

    function generarDistribucion2(nx, my, resul, si, sd, ii, id, superior, inferior, izquierda, derecha)
        integer(4), intent(in) :: nx, my
        type(condicion), dimension(:), intent(in) :: superior, inferior, izquierda, derecha
        real(8), intent(in) :: resul(1:(nx - 1) * (my - 1)), si, sd, ii, id
        real(8) generarDistribucion2(1:my+1, 1:nx+1)
        integer(4) i, offset

        ! Interno
        offset = 0
        do i = 2, my
            generarDistribucion2(i, 2:nx) = resul(1+offset:nx-1)
            offset = offset + nx - 1
        end do

        ! Esquinas
        generarDistribucion2(1, 1) = si
        generarDistribucion2(1, nx+1) = sd
        generarDistribucion2(my+1, 1) = ii
        generarDistribucion2(my+1, nx+1) = id

        ! Bordes
        do i = 2, nx
            if (superior(i - 1)%tipo == DIRICHLET) then
                generarDistribucion2(1, i) = superior(i - 1)%valor
            elseif (superior(i - 1)%tipo == NEUMANN) then
                generarDistribucion2(1, i) = generarDistribucion2(3, i)
            end if

            if (inferior(i - 1)%tipo == DIRICHLET) then
                generarDistribucion2(my + 1, i) = inferior(i - 1)%valor
            elseif (inferior(i - 1)%tipo == NEUMANN) then
                generarDistribucion2(my + 1, i) = generarDistribucion2(my - 1, i)
            end if
        end do

        do i = 2, my
            if (izquierda(i - 1)%tipo == DIRICHLET) then
                generarDistribucion2(i, 1) = izquierda(i - 1)%valor
            elseif (izquierda(i - 1)%tipo == NEUMANN) then
                generarDistribucion2(i, 1) = generarDistribucion2(i, 3)
            end if

            if (derecha(i - 1)%tipo == DIRICHLET) then
                generarDistribucion2(i, nx + 1) = derecha(i - 1)%valor
            else if (derecha(i - 1)%tipo == NEUMANN) then
                generarDistribucion2(i, nx + 1) = generarDistribucion2(i, nx - 1)
            end if
        end do
    end function generarDistribucion2

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

end module EDDP