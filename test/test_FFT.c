/**
 *  Al llamar a la FFT con un vector de 0's, debe devolver un vector de 0's
 *  El error cometido por la función de raiz cuadrada debe ser menor al 1/1000 del valor esperado
 *  Al ingresar un vector simétrico en la componente real, la parte imaginaria de la FFT es nula
 *  La FFT de una senial continua tiene que valer 0 en todas las componentes que no sean la primera
 *  El balance de energías en el tiempo y frecuencia tiene que cumplir la identidad de parceval
 *  Dado un arreglo, se debe poder obtener la posicion de maxima energía del mismo
 */
#include "unity.h"
#include "FFT.h"

#define TamanioVectores 2048
#define testNumber 527.0
#define sqrt527 22.9564806

float real[TamanioVectores];
float imaginario[TamanioVectores];
float zeros[TamanioVectores];
float salida[TamanioVectores/2];
const int posicion = 73;


void InicializarVector(float *vector, float valor){
    for (uint16_t i=0; i<TamanioVectores; i++){
        vector[i] = valor;
    }
}

void setUp(void){
    InicializarVector(zeros, 0);
    InicializarVector(real, 0);
    InicializarVector(imaginario, 0);
}


void tearDown(void){
}

// Al llamar a la FFT con un vector de 0's, debe devolver un vector de 0's
void test_FFTDeCeros(void){
    FFT(real, imaginario);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(zeros, real, TamanioVectores);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(zeros, imaginario, TamanioVectores);
}

// El error cometido por la función de raiz cuadrada debe ser menor al 1/1000 del valor esperado
void test_ErrorRaizCuadrada(){
    float optimizada = squareRootOptimized(527);
    TEST_ASSERT_FLOAT_WITHIN(sqrt527/1000, sqrt527, optimizada);
}

// Al ingresar un vector simétrico en la componente real, la parte imaginaria de la FFT es nula
void test_VectorRealSimetrico(){
    InicializarVector(real, 4);
    FFT(real, imaginario);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(zeros, imaginario, TamanioVectores);
}

// La FFT de una senial continua tiene que valer 0 en todas las componentes que no sean la primera
void test_VectorContinuo(){
    InicializarVector(real, 4);
    FFT(real, imaginario);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(zeros, real+1, TamanioVectores-1);
}

// El balance de energías en el tiempo y frecuencia tiene que cumplir la relación de parceval
// Para cumplir esto, es necesario hacer el escalado de los valores obtenidos de la FFT con las funciones de la biblioteca
void test_ConservacionEnergia(){
    float sumaEntrada = 0.0;
    float sumaSalida = 0.0;
    InicializarVector(real, 4);
    // Calculo la energía RMS en la entrada
    for (uint16_t i = 0; i<TamanioVectores; i++){
        sumaEntrada += real[i] * real[i] / 2;
    }
    sumaEntrada /= TamanioVectores;
    
    FFT(real, imaginario);
    FFTScaleHalf(real, imaginario);
    FFTModuleHalf(real, imaginario, salida);
    for (uint16_t i = 0; i<TamanioVectores/2; i++){
        sumaSalida += salida[i] * salida[i];
    }
    TEST_ASSERT_FLOAT_WITHIN(sumaSalida/1000, sumaSalida, sumaEntrada);    
}

// Dado un arreglo, se debe poder obtener la posicion de maxima energía del mismo
void test_FrecuenciaMaxima(void){
    float resultado = 0;
    real[posicion] = testNumber;
    FFTMaxFrecuency(real, &resultado);
    TEST_ASSERT_EQUAL_FLOAT(posicion, resultado);
}
