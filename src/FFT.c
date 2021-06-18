#include "FFT.h"

// Takes the real and imaginary vectors from the original domain, and gives back the FFT of the input, split into real and imaginary parts.
// Does not compute the module of the output, so output information is in the form of complex vectors instead of module and fase.
// A scaling is done taking into account the input size, since the algorithm itself scales the input as it progresses.
void FFT(float real[MaxValues], float imag[MaxValues]) {
	uint16_t i = 0;
	uint16_t numOfProblems = 1;
	uint16_t problemSize = MaxValues;
	uint16_t halfSize = 0;
	uint16_t valueLast, valueFirst, twiddleFactor;
	float temporalReal, temporalImaginary, auxReal, auxImaginary;
	uint16_t j, k;
	
	// Radix2 fft algorithm
	while (problemSize > 1){
		halfSize = problemSize/2;
		for (k=0; k<=(numOfProblems - 1); k++){
			valueFirst = k * problemSize;
			valueLast = valueFirst + halfSize - 1;
			twiddleFactor = 0;
			for (j=valueFirst; j<=valueLast; j++){
				temporalReal = real[j];
				temporalImaginary = imag[j];
				real[j] = temporalReal + real[j + halfSize];
				imag[j] = temporalImaginary + imag[j + halfSize];
				auxReal = real[j + halfSize];
				auxImaginary = imag[j + halfSize];
				real[j + halfSize] = (temporalReal - auxReal) * cos_table[twiddleFactor] - (temporalImaginary - auxImaginary) * sin_table[twiddleFactor];
				imag[j + halfSize] = (temporalReal - auxReal) * sin_table[twiddleFactor] + (temporalImaginary - auxImaginary) * cos_table[twiddleFactor];
				twiddleFactor = twiddleFactor + numOfProblems;
			}
		}
		numOfProblems = 2 * numOfProblems;
		problemSize = halfSize;
	}

	// Reorder values so they make sense.
	// The output of the algorithm previously used has the results in bit reversed order (i.e. 100 result is on position 001).
	for(i=0; i<MaxValues; i++){
		j = reverseBits(i);
		if (j>i){	// Only swap the values which weren't already swapped
			temporalReal = real[i];
			real[i] = real[j];
			real[j] = temporalReal;
			temporalImaginary = imag[i];
			imag[i] = imag[j];
			imag[j] = temporalImaginary;
		}
		if (i==j){
			real[i] = real[i];
			imag[i] = imag[i];
		}
	}
}

// Computes the Inverse FFT of the input vectors, to return to the time domain.
// The output is swaped, so the real vector is actually the imaginary part, and vice versa.
void IFFT(float real[MaxValues], float imag[MaxValues]){
	FFT(imag, real);
}

// Scales the real and imaginary vectors taking into account the input vector size
void FFTScale(float real[MaxValues], float imag[MaxValues]){
	real[0] = real[0] / ScalingFactor;
	imag[0] = imag[0] / ScalingFactor;
	for (uint16_t i=1; i<MaxValues; i++){
		real[i] = (real[i] * 2) / ScalingFactor;
		imag[i] = (imag[i] * 2) / ScalingFactor;
	}
}

// Scales the real and imaginary vectors taking into account the input vector size, and a minimum value to consider it valid.
void FFTScaleTruncated(float real[MaxValues], float imag[MaxValues], float offset, float underfrequencyValues){
	real[0] = real[0] / ScalingFactor;
	imag[0] = imag[0] / ScalingFactor;
	for (uint16_t i=1; i<MaxValues; i++){
		real[i] = (real[i] * 2) / ScalingFactor;
		imag[i] = (imag[i] * 2) / ScalingFactor;
	}
	for (uint16_t i=0; i<MaxValues; i++){
		if (i < underfrequencyValues){
			real[i] = 0;
			imag[i] = 0;
		} else {
			if (absuoluteValue(real[i]) < offset){
				real[i] = 0;
			} else if (real[i] > 0){
				real[i] = real[i] - offset;
			} else{
				real[i] = real[i] + offset;
			}
			if (absuoluteValue(imag[i]) < offset){
				imag[i] = 0;
			} else if (imag[i] > 0){
				imag[i] = imag[i] - offset;
			} else{
				imag[i] = imag[i] + offset;
			}
		}
	}
}

// Scales the real and imaginary vectors taking into account the input vector size
void FFTScaleHalf(float real[MaxValues], float imag[MaxValues]){
	real[0] = real[0] / ScalingFactor;
	imag[0] = imag[0] / ScalingFactor;
	for (uint16_t i=1; i<MaxValues/2; i++){
		real[i] = (real[i] * 2) / ScalingFactor;
		imag[i] = (imag[i] * 2) / ScalingFactor;
	}
}

// Scales the real and imaginary vectors. 
// Since the vectors were scaled by 1/sqrt(2) when doing the FFT, they need to be resized by multiplying them by sqrt(2)
void IFFTScale(float real[MaxValues], float imag[MaxValues]){
	real[0] = real[0];
	imag[0] = imag[0];
	for (uint16_t i=1; i<MaxValues; i++){
		real[i] = real[i];
		imag[i] = imag[i];
	}
}

// Receives a real and imaginary vector, and computes over the output vector the square root of the sum of input vector squares.
void FFTModule(float real[MaxValues], float imag[MaxValues], float *output){
	uint16_t i=0;
	for (i=0; i<MaxValues; i++){
		output[i] = squareRootOptimized((real[i]*real[i])+(imag[i]*imag[i]));
	}
}

// Receives a real and imaginary vector, and computes over the output vector the square root of the sum of input vector squares (up to half the range).
void FFTModuleHalf(float real[MaxValues], float imag[MaxValues], float *output){
	uint16_t i=0;
	for (i=0; i<MaxFFTValues; i++){
		output[i] = squareRootOptimized((real[i]*real[i])+(imag[i]*imag[i]));
	}
}

// Computes the maximum frequency from the vector input, and loads it into the peakFrequency pointer.
void FFTMaxFrecuency(float *vectorInput, float *peakFrequency){
	uint16_t i, j;
	j = 1;
	for (i=1; i<(MaxFFTValues); i++){
		if (vectorInput[i] > vectorInput[j]) j=i;
	}
	*peakFrequency = j*1.0;
}

// Computes the reverse bit to bit from the input number (for example, 1000b transforms into 0001b).
uint16_t reverseBits(uint16_t x) {
	uint16_t result = 0;
	for (int i = 0; i < 11; i++, x >>= 1)
		result = (result << 1) | (x & 0b1);
	return result;
}

float absuoluteValue(float input){ return (input>0?input:-1*input);}

// Computes the square root of a function, optimized for speed.
float squareRootOptimized (float input){
	int32_t exponent;
	int8_t sign;
	float aux,accumulator,xPower;
	
	if(!(input>0.0)) return 0.0;

	// Extracts the mantissa and exponent from the original float
	exponent = ((*((uint32_t *)&input)) & 0x7F800000)>>23;
	exponent = exponent - 127;
	if (exponent<0) sign = -1; else sign = 1;
	aux = ((*((uint32_t *)&input)) & 0x007FFFFF); // Take the mantissa bits from the original float
	aux = aux * 1.192092895507812e-07; // Divide it by 2**23 to get the value of mantissa-1
	
	accumulator =  1.0 + 0.49959804148061*aux;
	xPower = aux*aux;
	accumulator = accumulator - 0.12047308243453*xPower;
	xPower = xPower * aux;
	accumulator = accumulator + 0.04585425015501*xPower;
	xPower = xPower * aux;
	accumulator = accumulator - 0.01076564682800*xPower;

	if (exponent & 0x00000001) accumulator = accumulator * ROOT2; // An odd input exponent means an extra sqrt(2) in the output

	if (sign>0){
		exponent = exponent/2;	
		return (accumulator * (1<<exponent));
	}else{
		exponent = (int32_t)((exponent+0.1)/2)-1;
		return (accumulator / (float)(1<<(-1*exponent)));
	}
}