
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <bitset>
#include <time.h>

using namespace std;

double pi = M_PI;

typedef complex<double> doublecomplex;

int reverseBits(unsigned int numberToReverse, int limit);
void DFT(doublecomplex *timeDomain, doublecomplex *frequencyDomain, int size);
void inv_FFT(doublecomplex *timeDomain, doublecomplex *frequencyDomain, int size); // TO DO
void FFT(doublecomplex *timeDomain, doublecomplex *frequencyDomain, int size);
void inv_FFT(doublecomplex *timeDomain, doublecomplex *frequencyDomain, int size); // TO DO
void recFFT(doublecomplex *timeDomain, doublecomplex *frequencyDomain, int size);
void inv_recFFT(doublecomplex *timeDomain, doublecomplex *frequencyDomain, int size); // TO DO
doublecomplex twiddle(int power, int N);
doublecomplex inv_twiddle(int power, int N);

doublecomplex twiddle(int power, int N){
	complex<double> W = complex<double>(cos(-2.0*pi/N),sin(-2.0*pi/N));
	return pow(W,power);
}
doublecomplex inv_twiddle(int power, int N){
	complex<double> W = complex<double>(cos(2.0*pi/N),sin(2.0*pi/N));
	return pow(W,power);
}

// Radix-2 Fast Fourier Transform of a Signal (1D)
// Non-recursive

void FFT(doublecomplex *timeDomain, doublecomplex *frequencyDomain, int size){
	int oddOffset, fftSize = 2, evenIndex, oddIndex, i, j;
	doublecomplex e, o, W_e, W_o;
	for (i=0; i < size; i++){
		*(frequencyDomain+i) = *(timeDomain+reverseBits(i,size));
	}
	while(fftSize <= size){ // will be power of 2
		oddOffset = fftSize/2;
		for (i=0; i<size; i+=fftSize){
			for(j=0;j<fftSize/2;j++){
				evenIndex = i+j;
				oddIndex  = evenIndex+oddOffset;
				e = *(frequencyDomain+evenIndex);
				o = *(frequencyDomain+oddIndex);
				W_e = twiddle( (evenIndex) % fftSize,fftSize);
				W_o = twiddle( (oddIndex)  % fftSize,fftSize);
				*(frequencyDomain+evenIndex) = e + W_e*o;				
				*(frequencyDomain+oddIndex)  = e + W_o*o;
			}
		}
		fftSize = fftSize*2;
	}	
}


// Radix-2 Fast Fourier Transform of a Signal 1D
// Recursive

void recFFT(doublecomplex *timeDomain, doublecomplex *frequencyDomain, int size){

	doublecomplex *even, *odd, *Feven, *Fodd, W;
	even  = (doublecomplex*)malloc(sizeof(doublecomplex)*(size/2));
	Feven = (doublecomplex*)malloc(sizeof(doublecomplex)*(size/2));
	odd   = (doublecomplex*)malloc(sizeof(doublecomplex)*(size/2));
	Fodd  = (doublecomplex*)malloc(sizeof(doublecomplex)*(size/2));

	if ( size == 1 ){
		*frequencyDomain = *timeDomain;
	}else{
		int limit = size/2, i;
		for (i=0; i<limit; i++){
			*(even+i) = *(timeDomain+2*i);
			*(odd+i) = *(timeDomain+2*i+1);
		}
		recFFT(even,Feven,size/2);	
		recFFT(odd,Fodd,size/2);
		for (i=0; i<size;i++){
			W = twiddle(i,size);	
			*(frequencyDomain+i) = *(Feven+i%(size/2))+W*(*(Fodd+i%(size/2)));
		}
	}
	free(even);
	free(odd);
	free(Feven);
	free(Fodd);
}

void DFT(doublecomplex *timeDomain, doublecomplex *frequencyDomain, int size){
	// twiddle Factor
	complex<double> W = complex<double>(cos(2.0*pi/size),sin(2.0*pi/size));
	// Loop through and calc each value in the frequency domain
	for (int freqIndex=0; freqIndex<size; freqIndex++){
		for (int timeIndex=0; timeIndex<size; timeIndex++){
			*(frequencyDomain+freqIndex) += (*(timeDomain+timeIndex)) * pow(W,freqIndex*timeIndex);
		}
	}
}

// Bit reversal Algorithm for integer indices
int reverseBits(unsigned int numberToReverse, int limit){
	unsigned int result=0, numberOfBits=0, temp;
	for (int j=0; j<sizeof(limit)*8; j++){
		if ( limit-1 >= pow(2,j) ){
			numberOfBits++;
		}
		else{
			break;
		}
	}
	for (int i = 0; i < numberOfBits; i++){
		temp = (numberToReverse & (1 << i));
		if(temp){
			result |= (1 << ((numberOfBits - 1) - i));
		}
  	}
	return result;
}

