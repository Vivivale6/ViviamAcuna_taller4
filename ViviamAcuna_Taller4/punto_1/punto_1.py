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

gausiana = gauss(t)


#Funcion que integra para normalizar los datos
def integrate(f,a, b, N):
    h = (b-a)/(N-1)
    contador = 0.0

    for i in range(N-1):

        x = a + i * h 
        contador = np.sum(0.5 * (f(t) + f(t+h)) * h)

    return contador
#Se normalizan los datos    
gauss_n= integrate(gauss, 0, 10, 10) 


#Funcion que retorna el nucleo o kernel
def nucleo (gauss_n):
    
    gauss_nu= gauss_n + 100000
    normalizado = gausiana/gauss_nu
    return normalizado[:, np.newaxis] * normalizado[np.newaxis, :]

nucleo_n= nucleo(gauss_n) 


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

                    #En este caso se usa coseno porque solo se tienen en cuenta los valores reales y no los imaginarios
                    Real[k,l] = Real[k,l] + A[m,n]*np.cos(-2.0*np.pi*(m*k/largo+n*l/ancho))
                    
    return Real

#Se tienen en cuenta las dos primeras componentes de la imagen    
imagen_nueva=Im1[:,:,0]

#Se hallan las transformadas del nucleo y de la imagen 
nucleo_tr= transformada2d(nucleo_n)
imagen_tr= transformada2d(imagen_nueva)

#Se multiplican las transformadas del kernel y de la imagen 
mult_trans = nucleo_tr*imagen_tr


#Funcion que realiza la transformada de fourier inversa
def transformadaInv2d(A):
    forma = np.shape(A)
    ancho = forma[1] 
    largo = forma[0] 
    
    
    Real = np.zeros((largo,ancho))
    
    for k in range(largo):
        
        for l in range(ancho):
            
             for m in range(largo):
                    
                for n in range(ancho):


                    #Se toma coseno porque solo se tienen en cuenta los valores reales
                    Real[k,l] = Real[k,l] + A[m,n]*np.cos(2.0*np.pi*(m*k/largo+n*l/ancho))
                    
    return Real
#Se realiza la transformada inversa a la multiplicacion de las transformadas de la imagen y el kernel
imagen_suavizada = transformadaInv2d(mult_trans)

#Se guarda la imagen
plt.imshow(imagen_suavizada)
plt.savefig("suave1.png")




