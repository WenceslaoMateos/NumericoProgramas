MODULE arreglos

IMPLICIT NONE

CONTAINS

SUBROUTINE cargarMatriz(A, N, M)
INTENT (IN) :: N, M
INTENT (OUT) :: A
INTEGER(4) N, M, i, j
REAL(8) A(N, M)

    DO i = 1, N
        DO j = 1, M
            WRITE(*, '(A, I2, A, I2, A)', ADVANCE='NO') 'Ingrese el elemento (', i, ', ', j, ') de la matriz: '
            READ(*, *) A(i, j)
        END DO
    END DO
END SUBROUTINE

SUBROUTINE mostrarMatriz(A, N, M)
INTENT (IN) :: N, M, A
INTEGER(4) N, M, i, j
REAL(8) A(N, M)

    DO i = 1, N
        DO j = 1, M
            WRITE(*, '(F10.5, A)', ADVANCE='NO') A(i, j), ' '
        END DO
        WRITE (*, *)
    END DO
END SUBROUTINE

SUBROUTINE leerMatriz(A, N, M, archivo)
INTENT (IN) :: N, M, archivo
INTENT (OUT) :: A
INTEGER(4) N, M, i, j
REAL(8) A(N, M)
CHARACTER (LEN=*) :: archivo

    OPEN(UNIT=2, FILE=archivo, ACCESS='SEQUENTIAL')
    DO i = 1, N
        DO j = 1, M
            READ(2, '(F10.5)', ADVANCE='NO') A(i, j)
        END DO
        READ(2, *)
    END DO
    CLOSE(2)
END SUBROUTINE

SUBROUTINE grabarMatriz(A, N, M, archivo)
INTENT (IN) :: N, M, A, archivo
INTEGER(4) N, M, i, j
REAL(8) A(N, M)
CHARACTER (LEN=*) :: archivo

    OPEN(UNIT=2, FILE=archivo, ACCESS='SEQUENTIAL', STATUS='REPLACE')
    DO i = 1, N
        DO j = 1, M
            WRITE(2, '(F10.5, A)', ADVANCE='NO') A(i, j), ' '
        END DO
        WRITE (2, *)
    END DO
    CLOSE(2, STATUS='KEEP')
END SUBROUTINE

SUBROUTINE cargarVector(V, N)
INTENT (IN) :: N
INTENT (OUT) :: V
REAL(8) V(N)
INTEGER(4) N, i
    
    DO i = 1, N
        WRITE(*, '(A, I3, A)', ADVANCE='NO') 'Ingrese el elemento ', i, ' del vector: '
        READ(*, *) V(i)
    END DO
END SUBROUTINE

SUBROUTINE mostrarVector(V, N)
INTENT (IN) :: N, V
REAL(8) V(N)
INTEGER(4) N, i
    
    DO i = 1, N
        WRITE(*, '(F10.5)') V(i)
    END DO
END SUBROUTINE

SUBROUTINE leerVector(V, N, archivo)
INTENT (IN) :: N, archivo
INTENT (OUT) :: V
REAL(8) V(N)
INTEGER(4) N, i
CHARACTER (LEN=*) :: archivo
    
    OPEN(UNIT=2, FILE=archivo, ACCESS='SEQUENTIAL')
    DO i = 1, N
        READ(2, '(F10.5)') V(i)
    END DO
    CLOSE(2)
END SUBROUTINE

SUBROUTINE grabarVector(V, N, archivo)
INTENT (IN) :: N, V, archivo
REAL(8) V(N)
INTEGER(4) N, i
CHARACTER (LEN=*) :: archivo
    
    OPEN(UNIT=2, FILE=archivo, ACCESS='SEQUENTIAL', STATUS='REPLACE')
    DO i = 1, N
        WRITE(2, '(F10.5)') V(i)
    END DO
    CLOSE(2, STATUS='KEEP')
END SUBROUTINE

END MODULE
