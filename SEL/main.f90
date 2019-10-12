program principal
    use SELs
    use arreglos

    implicit none

    integer, parameter :: orden = 4
    real(8) matriz(orden, orden), term_ind(orden, 1), xini(orden, 1), matriz2(3, 3)

    xini = 0
    call leerMatriz(matriz, "matriz_ejemplo1.txt")
    call leerMatriz(term_ind, "independientes_ejemplo1.txt")

    call mostrarMatriz(matriz)
    write(*, *)
    call mostrarMatriz(term_ind)
    write(*, *)

    write(*, *) "Norma M ", mNormaM(matriz)
    write(*, *) "Norma L ", mNormaL(matriz)
    write(*, *) "Norma F ", mNormaF(matriz)
    write(*, *)

    call pivotear(matriz, term_ind)
    call mostrarMatriz(matriz)
    write(*, *)
    call mostrarMatriz(term_ind)
    write(*, *)

    write(*, *) "Solucion Gauss:"
    call mostrarMatriz(solucionGauss(matriz, term_ind))
    write(*, *)

    write(*, *) "Solucion Gauss iterativo:"
    call mostrarMatriz(refinamientoIter(matriz, term_ind, 0.000000000001_8, solucionGauss, mNormaM))
    write(*, *)

    write(*, *) "Jacobi:"
    call mostrarMatriz(jacobi(matriz, term_ind, xini, 0.000000000001_8))
    write(*, *)

    write(*, *) "Gauss Seidel:"
    call mostrarMatriz(gaussSeidel(matriz, term_ind, xini, 0.000000000001_8))
    write(*, *)

    matriz2(1, 1) = 1.
    matriz2(1, 2) = 0.
    matriz2(1, 3) = 0.
    matriz2(2, 1) = 2.
    matriz2(2, 2) = 1.
    matriz2(2, 3) = 0.
    matriz2(3, 1) = -1.
    matriz2(3, 2) = -2.
    matriz2(3, 3) = -1.
    write(*, *) "Matriz original:"
    call mostrarMatriz(matriz2)
    write(*, *)
    write(*, *) "Matriz inversa:"
    call mostrarMatriz(matrizInversa(matriz2))
    write(*, *)
    write(*, *) "Chequeo matriz inversa:"
    call mostrarMatriz(matmul(matriz2, matrizInversa(matriz2)))
contains

end program principal
