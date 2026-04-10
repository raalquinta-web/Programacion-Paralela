#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define N 100
#define dt 200 //
//Ecuacion de onda usando diferencia central
 
void onda(double u_ant[N], double u_act[N], double u_sig[N], double c){
	
	double c2 = c*c;
	//Aqui calculamos los puntos corregidos con u_sig[]
	for (int i=1; i<N-1; i++){
		 
		 u_sig[i] = 2.0 * u_act[i] - u_ant[i] + c2*(u_act[i+1] - 2.0*u_act[i] + u_act[i-1]);
	}
	
	//Condiciones
	u_sig[0] = 0.0;
	u_sig[N-1] = 0.0;
}


int main(){
	double u_ant[N] = {0};
	double u_act[N] = {0};
	double u_sig[N] = {0};
	double cte = 1.0;
	
	//Una onda que empieza en el centro	
	u_act[N/2]=1.0;
	u_ant[N/2]=1.0;
	//Aqui abrimos un archivo csv con el nombre que querramos para guardar los datos
	FILE *archivo = fopen("datos_onda.csv", "w");
	if (archivo ==NULL){
		printf("Error al crear el archivo. \n");
		return 1;
	}
	
	// Aqui comienza el bucle espacial de la simulacion
	for (int t=0; t<dt; t++){
		//aqui calculamos u_sig[]
		onda(u_ant,u_act,u_sig,cte);
		//aqui se guardan los datos actuales 
		for (int i=0; i<N; i++){
			fprintf(archivo, "%f", u_act[i]);
			if (i<N-1){
				fprintf(archivo, ",");
			}
		}
		fprintf(archivo, "\n"); //salto de linea para la sig iteracion
		
		//Aqui actualizamos los arreglos ya que el tiempo avanza
		for (int i=0; i<N; i++){
			u_ant[i] = u_act[i];
			u_act[i] = u_sig[i];
		}
	
	}
	
	fclose(archivo);
	printf("Simulacion completa con pulso en el medio. Datos guardados en 'datos_onda.csv'\n");	
	return 0;
}
