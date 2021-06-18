# TP4-TSSE
El código de prueba se realizó sobre una biblioteca propia de FFT.

El objetivo de los test unitarios diseñados es probar distintas propiedades que debe tener una FFT, y servirá para validar nuevas implementaciones en la biblioteca al probar nuevos algoritmos.

Las pruebas elegidas para realizar son:
 *  Al llamar a la FFT con un vector de 0's, debe devolver un vector de 0's
 *  El error cometido por la función de raiz cuadrada debe ser menor al 1/1000 del valor esperado
 *  Al ingresar un vector simétrico en la componente real, la parte imaginaria de la FFT es nula
 *  La FFT de una senial continua tiene que valer 0 en todas las componentes que no sean la primera
 *  El balance de energías en el tiempo y frecuencia tiene que cumplir la identidad de parceval
 *  Dado un arreglo, se debe poder obtener la posicion de maxima energía del mismo
 
