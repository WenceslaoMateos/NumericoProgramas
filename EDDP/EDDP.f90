module EDDP
    implicit none
    
    contains
    
    subroutine generarSistema(mat, term_ind, x0, x1, y0, y1, nx, my, superior, inferior, izquierda, derecha)
        integer(4), intent(in) :: nx, my
        real(8), intent(in) :: x0, x1, y0, y1
        real(8), dimension(:,:), intent(out) :: mat, term_ind
        real(8), dimension(:), intent(in) :: superior, inferior, izquierda, derecha
        integer(4) filas, columnas, desde, hasta, i, n, m
        real(8) h, k

        h = (x1 - x0) / nx
        k = (y1 - y0) / my
        n = nx - 1
        m = my - 1
        mat = 0.
        term_ind = 0. * h**2. * k**2
        filas = size(mat, dim=1)
        columnas = size(mat, dim=2)

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

    subroutine grabarDatos(distribucion, x0, x1, y0, y1, nx, my, archivo)
        intent(in) :: distribucion, x0, x1, y0, y1, nx, my, archivo
        integer(4) nx, my, i, j
        real(8) distribucion(:, :), x0, x1, y0, y1, h, k, x, y
        character(len=*) archivo

        open(unit=2, file=archivo, access='SEQUENTIAL', status='REPLACE')
        h = (x1 - x0) / nx
        k = (y1 - y0) / my
        y = y0
        do i = 1, my + 1
            x = x0
            do j = 1, nx + 1
                write(2, *) x, y, distribucion(i, j)
                x = x + h
            end do
            write(2, *)
            y = y + k
        end do
        close(2, status='KEEP')
    end subroutine

end module EDDP