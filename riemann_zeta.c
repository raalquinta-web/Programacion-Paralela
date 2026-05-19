#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

 double riemann_zeta(double s, int k) {
    double result = 0.0;
    for (int i = 1; i < k; i++) {
        for (int j = 1; j < k; j++) {
            result += (2*(i&1)-1) / pow(i+j, s);
        }
    }
    return result*pow(2,s);
 }

int main() {
    unsigned int N=1000;

    double *X = (double*)malloc(N*sizeof(double));

    printf("Calculando la funcion zeta de Riemann en paralelo...\n");


    #pragma omp parallel for shared(X, N) schedule(guided)
    for (unsigned int k = 0; k < N; k++) {
        X[k] = riemann_zeta(2.0, k);
    }

    // Sabemos que Zeta(2) converge a pi^2 / 6 (aprox 1.644934)
    printf("El valor aproximado para k = %d es: %f\n", N - 1, X[N - 1]);

    // Liberamos la memoria
    free(X);

    return 0;
}

// Apartado (i)

// No existen dependencias de datos entre las iteraciones porque cada iteración calcula 
// un valor independiente y lo escribe en una posición  de memoria única (X[k]). Las variables
// compartidas son el puntero del arreglo X y el límite superior N. El iterador del bucle k 
// es privado para cada hilo.

// Apartado (ii)

// Dado que la complejidad temporal de la función evaluada es O(k^2), la cantidad de trabajo 
// computacional por iteración no es constante, sino que crece cuadráticamente. 
// Las iteraciones iniciales (con valores pequeños de k) finalizan casi instantáneamente, 
// mientras que las últimas iteraciones (con valores grandes de k) acaparan mucho tiempo de CPU.
// Si el trabajo se divide ingenuamente en partes iguales, se producirá un desbalanceo severo:
// los hilos a cargo de las primeras iteraciones terminarán rápido y quedarán inactivos (idle), 
// esperando a que los hilos con las iteraciones más pesadas terminen su trabajo.

// Apartado (iii)

// Static scheduling no es adecuado para este caso debido a la naturaleza no uniforme del 
// trabajo en cada iteración. Con static scheduling, cada hilo recibe un bloque fijo de 
// iteraciones, lo que puede resultar en un desbalanceo significativo. 

// Dynamic scheduling es una mejor opción, ya que asigna iteraciones a los hilos de manera 
// dinámica a medida que terminan su trabajo. Esto permite que los hilos que terminan 
// rápidamente puedan tomar nuevas iteraciones, reduciendo el tiempo de inactividad y 
// mejorando la eficiencia general del programa.

// Guided scheduling también es una buena opción, ya que asigna iteraciones en bloques que
// disminuyen de tamaño a medida que se asignan, lo que puede ayudar a equilibrar la carga 
// de trabajo a medida que se procesan iteraciones más pesadas. Sin embargo, 
// dynamic scheduling es más flexible y puede adaptarse mejor a la variabilidad extrema 
// en el tiempo de ejecución de las iteraciones. 





