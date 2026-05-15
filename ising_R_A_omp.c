#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


#define L 32

// Funcion para llenar la red con spines aleatorios (-1 o 1)
void inicializar_red(int red[L][L]) {
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            // rand() % 2 da 0 o 1
            // Si rand() % 2 es 0, entonces (rand() % 2) * 2 - 1 es -1
            // Si rand() % 2 es 1, entonces (rand() % 2) * 2 - 1 es 1   
            red[i][j] = (rand() % 2) * 2 - 1; // Esto genera 0 o 1, luego lo convierte a -1 o 1
        }
    }
}

// Funcion para imprimir la red en pentalla y ver que funciona
//void imprimir_red(int red[L][L]) {
//    for (int i = 0; i < L; i++) {
//        for (int j = 0; j < L; j++) {
//            // Imprimimos un + si es 1 y un - si es -1 para una mejor visualización
//            if (red[i][j] == 1) {
//                printf("+ ");
//            } else {
//                printf("- ");
//            }
//        }
//        printf("\n"); // salto de linea al terminar cada fila
//    }
//   printf("\n"); // salto de linea adicional para separar la red de otras salidas
//}

// Funcion para calcular la energía promedio por sitio
double calcular_energia(int red[L][L]) {

    double energia_total = 0;
    #pragma omp parallel for collapse(2) reduction(+:energia_total)
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            // Miramos solo a la derecha y abajo (con efecto Pac-Man) para no contar doble
            int derecha = red[i][(j + 1) % L];
            int abajo   = red[(i + 1) % L][j];
            
            // Si son iguales, da 1. Si son distintos, da -1. 
            // Como la fórmula tiene un signo negativo al inicio (H = -J * ...), multiplicamos por -1.0
            energia_total += -1.0 * red[i][j] * (derecha + abajo);
        }
    }
    
    // Dividimos por el número total de sitios para tener la energía promedio
    return energia_total / (L * L);
}

// Funcion para calcular la magnetizacion promedio por sitio
double calcular_magnetizacion(int red[L][L]) {
    int suma_spines = 0;

    #pragma omp parallel for collapse(2) reduction(+:suma_spines)
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            suma_spines += red[i][j];
        }
    }
    // La magnetización promedio por sitio es la suma total de spines dividida por el número total de sitios
    return (double)suma_spines / (L * L);
}

// Funcion que calcula cuanta energia cuesta voltear el spin en (i,j)
int calcular_delta_E(int red[L][L], int i, int j) {
    // El truco (i - 1 + L) % L asegura que si i es 0, el vecino de arriba sea L - 1
    int arriba = red[(i - 1 + L) % L][j];
    int abajo = red[(i + 1) % L][j];
    int izquierda = red[i][(j - 1 + L) % L];
    int derecha = red[i][(j + 1) % L];

    int suma_vecinos = arriba + abajo + izquierda + derecha;

    //J = 1 en este modelo
    int delta_E = 2 * red[i][j] * suma_vecinos; // El factor de 2 viene de la fórmula ΔE = 2 * J * s_ij * sum_neighbors
    return delta_E;
}

// Funcion que da un solo paso del algoritmo de Metrópolis
void paso_metropolis(int red[L][L], double T) {
    // 1. Elegir un spin al azar
    int i = rand() % L;
    int j = rand() % L;

    // 2. Calcular la diferencia de energía
    int dE = calcular_delta_E(red, i, j);

    // 3. Reglas de Metrópolis
    if (dE < 0) {
        red[i][j] = -red[i][j]; // Aceptamos el cambio siempre
    } else {
        // Calculamos la probabilidad de aceptarlo aunque cueste energía
        double p = exp(-dE / T);
        
        // Generamos un decimal aleatorio entre 0.0 y 1.0
        double r = (double)rand() / RAND_MAX; 
        
        if (r < p) {
            red[i][j] = -red[i][j]; // Lo aceptamos por suerte termodinámica
        }
    }
}

// Funcion que calcula la correlación espacial para una distancia R
double calcular_correlacion(int red[L][L], int R, double mag_promedio) {
    double suma_correlacion = 0;
    
    #pragma omp parallel for collapse(2) reduction(+:suma_correlacion)
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            int spin_actual = red[i][j];
            
            // Tomamos el spin a una distancia R hacia la derecha (con efecto Pac-Man)
            int spin_vecino = red[i][(j + R) % L]; 

            suma_correlacion += (spin_actual * spin_vecino);
        }
    }

    // Promedio de (S1 * S2) en toda la red
    double promedio_S1_S2 = suma_correlacion / (L * L);

    // Aplicamos la fórmula del PDF: C(R) = <S1*S2> - <S>^2 
    // Nota: mag_promedio es <S>, por lo que <S>^2 es mag_promedio * mag_promedio
    return promedio_S1_S2 - (mag_promedio * mag_promedio);
}

int main() {
    srand(time(NULL));
    int red_spines[L][L];
    int N_total = L * L;
   
    printf("Iniciando simulaciones del Modelo de Ising...\n");
    printf("Por favor espera, esto tomara un par de minutos al hacer cálculos precisos.\n\n");
   

    // Ejercicio 1: Tiempo de Equilibración

    printf("1/3 - Calculando tiempo de equilibracion (T=2.0)...\n");
    
    // Abrimos un archivo para escribir los datos ("w" significa write)
    FILE *archivo_eq = fopen("datos_eq.txt", "w");
    fprintf(archivo_eq, "Iteracion\tMagnetizacion\tEnergia\n");   
    
    
    inicializar_red(red_spines);
    double T = 2.0;
    
    // Aquí medimos iteración por iteración para ver cómo cae al inicio
    
    for (int k = 0; k <= 500000; k++) {
        if (k % 100 == 0) { // Guardamos una muestra cada 100 pasos
            fprintf(archivo_eq, "%d\t%f\t%f\n", k, calcular_magnetizacion(red_spines), calcular_energia(red_spines));
        }
        paso_metropolis(red_spines, T);
    }
    fclose(archivo_eq);


    // Ejercicio 2: Transición de Fase (Barrido de Temperaturas)

    printf("2/3 - Realizando barrido de Temperaturas (T=1.0 a T=4.0)...\n");
    
    FILE *archivo_fase = fopen("datos_fase.txt", "w");
    fprintf(archivo_fase, "T\t<M>\t<E>\tCv\tChi\n");
    
    // Barrido de temperaturas desde 1.0 hasta 4.0

    for (double T_fase = 1.0; T_fase <= 4.0; T_fase += 0.1) {
        inicializar_red(red_spines);

        // Termalizacion 
        int termalizacion = L*L;
        
        for (int k = 0; k < termalizacion; k++) {
            paso_metropolis(red_spines, T_fase);
        }

        // Medicion
        int medicion = 5000;
        double sum_M = 0, sum_M2 = 0, sum_E = 0, sum_E2 = 0;
        for (int m = 0; m < medicion; m++) {
            for (int k = 0; k < N_total; k++) {
                paso_metropolis(red_spines, T_fase);
            }
            // Tomamos valor absoluto de M para evitar cancelaciones
            double mag = fabs(calcular_magnetizacion(red_spines));
            double ene = calcular_energia(red_spines);

            sum_M += mag; sum_M2 += (mag * mag);
            sum_E += ene; sum_E2 += (ene * ene);
        }

        // Promedios y fluctuaciones
        double avg_M = sum_M / medicion;
        double avg_M2 = sum_M2 / medicion;
        double avg_E = sum_E / medicion;
        double avg_E2 = sum_E2 / medicion;

        double var_M = avg_M2 - (avg_M * avg_M);
        double var_E = avg_E2 - (avg_E * avg_E);

        double Cv = (var_E * N_total) / (T_fase * T_fase);
        double Chi = (var_M * N_total) / T_fase;

        fprintf(archivo_fase, "%f\t%f\t%f\t%f\t%f\n", T_fase, avg_M, avg_E, Cv, Chi);
    }
    fclose(archivo_fase); 

    // Ejercicio 3: Longitud de Correlación

    printf("3/3 - Calculando longitud de correlacion (T=2.3 y T=3.0)...\n");
    
    FILE *archivo_corr = fopen("datos_corr.txt", "w");
    fprintf(archivo_corr, "R\tC(R)_Tc_2.3\tC(R)_Hot_3.0\n");

    double sum_C_Tc[L/2 + 1] = {0};
    double sum_C_Hot[L/2 + 1] = {0};
    int mcs_corr_medicion = 5000;

    // Prueba para T = 2.3 (Cerca del punto crítico) 
    double T_c1 = 2.3;
    inicializar_red(red_spines);

    for(int k=0; k < 2000 * N_total; k++) paso_metropolis(red_spines, T_c1); 

    for(int m=0; m < mcs_corr_medicion; m++) {
        for(int k=0; k < N_total; k++) paso_metropolis(red_spines, T_c1);
        double mag = calcular_magnetizacion(red_spines);
        for(int R=1; R<=L/2; R++) sum_C_Tc[R] += calcular_correlacion(red_spines, R, mag);
    }

    // Prueba para T = 3.0
    double T_c2 = 3.0;
    inicializar_red(red_spines);

    for(int k=0; k < 2000 * N_total; k++) paso_metropolis(red_spines, T_c2); 

    for(int m=0; m < mcs_corr_medicion; m++) {
        for(int k=0; k < N_total; k++) paso_metropolis(red_spines, T_c2);
        double mag = calcular_magnetizacion(red_spines);
        for(int R=1; R<=L/2; R++) sum_C_Hot[R] += calcular_correlacion(red_spines, R, mag);
    }

    // Guardamos los datos de ambas temperaturas en el archivo
    for(int R=1; R<=L/2; R++) {
        fprintf(archivo_corr, "%d\t%f\t%f\n", R, sum_C_Tc[R]/mcs_corr_medicion, sum_C_Hot[R]/mcs_corr_medicion);
    }
    fclose(archivo_corr);

    printf("\n¡Simulaciones terminadas con exito!\n");
    printf("Revisa los archivos 'datos_eq.txt', 'datos_fase.txt' y 'datos_corr.txt'\n");

    return 0;
}