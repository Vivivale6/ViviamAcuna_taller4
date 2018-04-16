import numpy as np
import matplotlib.pyplot as plt

#Se ingresa la imagen, esta debe estar en la misma carpeta del archivo
Im1 = plt.imread("viviam.png")


#En este paso se toma la decision de cortar la imagen debido a que si son muchos pixeles, al momento de hacer el recorrido se toma basante tiempo. 
#En caso de escoger una imagen pequena este se puede omitir
Im1 = Im1[1:21,1:21,:]


m,n,I = np.shape(Im1)

#El linspace tiene los mismos puntos que el tamano de la gausiana
t = np.linspace(-10, 10, 20)

#Funcion que retorna la gausiana
def gauss(t):
    return 1.0 - np.exp(0.1*t*t)

#Funcion que retorna la gaussiana inversa
def gauss_inver(t):
    return np.exp(-0.1*t*t)

#Funcion que integra    
def integrate(f,a, b, N):
    h = (b-a)/(N-1)
    suma = 0.0

    for i in range(N-1):
        x = a + i * h 

        suma = np.sum(0.5 * (f(t) + f(t+h)) * h)
    return suma
#Funcion del nucleo o kernel
def fundamental1(gauss_nbaja):
    gauss_nubaja= gauss_nbaja + 100000
    normalizado_baja = gausiana_bajas/gauss_nubaja

    return normalizado_baja[:, np.newaxis] * normalizado_baja[np.newaxis, :]

def fundamental2(gauss_nalta):
    normalizado_alta = gausiana_altas/gauss_nalta

    return normalizado_alta[:, np.newaxis] * normalizado_alta[np.newaxis, :]


#Funcion que realiza la transformada de fourier
def transformada2d(A):
    forma = np.shape(A)
    ancho = forma[1] 
    largo = forma[0] 
    
    
    Real = np.zeros((largo,ancho))
    
    for k in range(largo):
        
        for l in range(ancho):
            
             for m in range(largo):
                    
                for n in range(ancho):

                    #Se toma solo el coseno porque solo se requieren los valores reales
                    Real[k,l] = Real[k,l] + A[m,n]*np.cos(-2.0*np.pi*(m*k/largo+n*l/ancho))
                    
    return Real    
#Funcion que realiza la transfromada inversa de fourier
def transformadaInv2d(A):
    forma = np.shape(A)
    ancho = forma[1] 
    largo = forma[0] 
    
    
    Real = np.zeros((largo,ancho))
    
    for k in range(largo):
        
        for l in range(ancho):
            
             for m in range(largo):
                    
                for n in range(ancho):

                    #Se toma solo el coseno porque solo se requieren los valores reales
                    Real[k,l] = Real[k,l] + A[m,n]*np.cos(2.0*np.pi*(m*k/largo+n*l/ancho))
                    
    return Real    
#Gausiana para frecuencias altas y bajas
gausiana_bajas = gauss(t)
gausiana_altas = gauss_inver(t)    

#Normalizacion para frecuencias bajas y altas
gauss_nbaja= integrate(gauss, 0, 10, 10) 
gauss_nalta= integrate(gauss_inver, 0, 10, 10) 

#Se hallan los valores de la matriz nucleo o  kernel

nucleo_nbajas= fundamental1(gauss_nbaja) 
nucleo_naltas= fundamental2(gauss_nalta)

#Se calculan los valores de la transformada para la imagen u para el nucleo o kernel
imagen_nuevabajas=Im1[:,:,0]
nucleo_trbajas= transformada2d(nucleo_nbajas)
imagen_trbajas= transformada2d(imagen_nuevabajas)

imagen_nuevaaltas=Im1[:,:,0]
nucleo_traltas= transformada2d(nucleo_naltas)
imagen_traltas= transformada2d(imagen_nuevaaltas)

#Se multiplican los valores de las transfromadas del nucleo y la imagen 
mult_transbajas = nucleo_trbajas*imagen_trbajas
mult_transaltas = nucleo_traltas*imagen_traltas

#Se realiza la transformada inversa de la multiplicacion anterior
imagen_bajas = transformadaInv2d(mult_transbajas)
imagen_altas = transformadaInv2d(mult_transaltas)

#Se guardan las imagenes de las frecuencias bajas y altas
plt.imshow(imagen_bajas)
plt.savefig("bajas.png")

plt.imshow(imagen_altas)
plt.savefig("altas.png")

#El filtro de la gaussiana elimina los altos y pasa los bajos
#El filtro de la gaussiana menos uno elimina los bajos y pasa los altos


