module EDDP
    implicit none
    
    contains

    subroutine generaMatriz(mat, term_ind, n, m, superior, inferior, izquierda, derecha)
        real(8), dimension(:,:), intent(out) :: mat, term_ind
        real(8), dimension(:), intent(in) :: superior, inferior, izquierda, derecha
        integer(4), intent(in) :: n, m
        integer(4) filas, columnas, desde, hasta, i

        mat = 0.
        term_ind = 0.
        filas = size(mat, dim=1)
        columnas = size(mat, dim=2)

        ! ---------------GENERACION DE TERMINOS INDEPENDIENTES--------------- !
        ! Condiciones superiores
        desde = 1
        hasta = n
        term_ind(desde:hasta, 1) = term_ind(desde:hasta, 1) + superior(:)

        ! Condiciones inferiores
        desde = n * (m - 1) + 1
        hasta = n * m
        term_ind(desde:hasta, 1) = term_ind(desde:hasta, 1) + inferior(:)

        ! Condiciones izquierda
        desde = 1
        hasta = n * (m - 1) + 1
        do i = desde, hasta, n
            term_ind(i, 1) = term_ind(i, 1) + izquierda(mod(i, n)+1)
        end do
        
        ! Condiciones derecha
        desde = n
        hasta = n*m     
        do i = desde, hasta, n
            term_ind(i, 1) = term_ind(i, 1) + derecha(i/n)
        end do

        ! ----------------------GENERACION DE DIAGONAL---------------------- !
        do i = 1, n * m
            mat(i, i) = 4.
        end do

        ! ------------------------GENERACION DE UNOS------------------------ !
        ! Banda superiores
        do i = n + 1, n * m
            mat(i, i - n) = -1.
        end do

        ! Banda inferiores
        do i = 1, n * (m - 1)
            mat(i, i + n) = -1.
        end do

        ! Banda izquierda
        do i = 2, n * m
            if (mod(i, n) /= 1) then
                mat(i, i - 1) = -1.
            end if
        end do

        ! Banda derecha
        do i = 1, n * m - 1
            if (mod(i, n) /= 0) then
                mat(i, i + 1) = -1.
            end if
        end do

    end subroutine generaMatriz
    
    subroutine generaMatriz2(mat, term_ind, x0, x1, y0, y1, nx, my, superior, inferior, izquierda, derecha)
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
    end subroutine

end module EDDP