#include <stdlib.h>
#include <complex>
#include<iostream>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <cstdlib>
        

using namespace std;

//Se declaran las funciones 
int contar_lineas(FILE *in);

float* lagrange(float *tim, float *x_in, int num_dat, float espacio);

float *fouriertransr( int num_dat, float* x_lag);

float *fouriertransi( int num_dat, float* x_lag);

//Funcion que cuenta el numero de lineas del archivo 
int contar_lineas(FILE *in){

  	int cuenta_datos;
  	int a;
	  	cuenta_datos = 0;
  	do{
    	a = fgetc(in);

    	if(a=='\n'){

      	cuenta_datos++;
    }

  }while(a!=EOF);
	
  fclose(in);
  
  return cuenta_datos;
}

//Se realiza la funcion de lagrange y se toma como puntero

float* lagrange(float *tim, float *x_in, int num_dat, float espacio){
	
	float *tim_corr;
	
	float *x_corr;
	
	float tim_nue;
	
	float lag;
	
	float lagran;


	tim_corr = (float*)std::malloc(num_dat*sizeof(float));
	
	x_corr = (float*)std::malloc(num_dat*sizeof(float));
	
	

	for (int i = 0; i < num_dat; ++i){

		tim_corr[i] = i*espacio;
	}


	for (int j = 0; j < num_dat; ++j){
	  
		tim_nue = tim_corr[j];
		
		for (int k = 0; k < num_dat; ++k){
		  
			lagran = 1.0;
			
			for (int l = 0; l < num_dat; ++l){
			  
				lagran = lagran * ((tim_nue- tim[l])/(tim[k]-tim[l]));
			}

			lag += lagran*x_in[k];
		}
		x_corr[j]= lag;
		
	}

	return x_corr;
}

//Se realiza la funcion de la transfromada de fourier para datos reales, se usa coseno

float *fouriertransr( int num_dat, float* x_lag){
  
	float pi = 3.14159;
	
	float *xn_real;

	xn_real = (float*)malloc(num_dat*sizeof(float));
	

	for (int i = 0; i < num_dat; ++i){
	  
		xn_real[i]=0.0;
		
		for (int j = 0; j < num_dat; ++j){
		  
			xn_real[j] += x_lag[j]*cos((-2*pi*i*j)/num_dat);			
		}

	}

	return xn_real;
}

// Se realiza la funcion de la transformada para datos imaginarios se usa seno 

float *fouriertransi( int num_dat, float* x_lag){
  
	float pi = 3.14159;
	
	float *xn_imag;
	

	xn_imag = (float*)std::malloc(num_dat*sizeof(float));

	for (int i = 0; i < num_dat; ++i){
	  
		xn_imag[i]=0.0;
		
		for (int j = 0; j < num_dat; ++j){
		  
			xn_imag[j] += x_lag[j]*cos((-2*pi*i*j)/num_dat);			
		}

	}

	return xn_imag;
}

//Se realiza la funcion main en donde tiene parametro de entrada el archivo de datos que se va a usar 

int main(int argc, char *argv[]){
  

  //Se le asigna un nombre al archivo y se pide que se abra en modo lectura
	char *nomb_archivo = argv[1];
	
	FILE *in = fopen (nomb_archivo, "r");
	
	int num_dat = contar_lineas(in);

	// Se declaran todas las variables 
	float *tim;
	
	float *x_in;
	
	int i;
	
	float espacio;
	
	float *x_lag;
	
	float *xf_re;
	
	float *xf_im;
	
	float espaciof;
	
	float *frecuencia;
	
	FILE *salida;

	

	tim = (float*)std::malloc(num_dat*sizeof(float));
	
	x_in = (float*)std::malloc(num_dat*sizeof(float));
	
	x_lag = (float*)std::malloc(num_dat*sizeof(float));
	
	xf_re = (float*)std::malloc(num_dat*sizeof(float));
	
	xf_im = (float*)std::malloc(num_dat*sizeof(float));
	
	frecuencia = (float*)std::malloc(num_dat*sizeof(float));
	
	//Se recorre el archivo en busca de las columnas uno y dos que corresponden al tiempo y la posicion
	
	for(i=0;i<num_dat;i++){

    fscanf(in, "%f %f\n", &tim[i], &x_in[i]);

	}

	//Se busca tener intervalos de t  conocidos 
	
	espacio = (tim[num_dat-1]-tim[0])/(num_dat-1);

	// Se utiliza el polinomio de lagrange para encontrar valores de x relacionados con los intervalos de t conocidos
	
	x_lag = lagrange(tim, x_in, num_dat, espacio);
	
	xf_re = fouriertransr(num_dat, x_lag);
	
	xf_im = fouriertransi(num_dat, x_lag);

	espaciof = ((1/espacio)/2)/((num_dat/2)-1);

	
	// Se escribe el archivo de salida
	
	//salida = fopen("transformada.txt", "w");

	// Se hallan las frecuencias
	
	for (int j = 0; j < num_dat; ++j){
		frecuencia[j]= j*espaciof;
		
		cout << &frecuencia << " " << &xf_re << " " << &xf_im << endl;
	
			}
	         
	
	return 0;

}
	


