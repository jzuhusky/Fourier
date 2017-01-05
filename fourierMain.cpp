
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

#include "fourier.h"

using namespace std;

int main(int argc, char** argv){

	if ( argc == 1 ){
		return 0;
	}

	ifstream timeSeriesSignal;
	ofstream dftout, fftout, recfftout, benchmarkout; 
	ofstream dftbench, fftbench, recfftbench;
	
	timeSeriesSignal.open(argv[argc-1]); // Input file is the last thing specified

	int c, dftflag=0, fftflag=0, recfftflag=0, timeResults=0, FS=1;
	char *arg, *id, *outfile;
	while ( (c=getopt(argc,argv,"m:if:")) != -1 ){
                switch(c){
		// Do Benchmarking->Print Time Stats
                case 'm':
			if (*optarg == 'd'){
				dftflag = 1;
			} 
			if (*optarg == 'f'){
				fftflag = 1; 
			} 
			if (*optarg == 'r'){
				recfftflag = 1; 
			} 
                        break;
		case 'i': 
		// Give this program executing an ID for outputfiles
			id = optarg;
			break;
		case 'f':
		//Set the sampling frequency->For graphical output purposes
			FS = atoi(optarg);
			cout << "FS: " << FS << endl;	
			break;
                default:
			fftflag=1;
                        break;
                }
        }

	doublecomplex *timeDomain, *frequencyDomain, *frequencyDomainFFT, *recfrequencyDomainFFT;

	int powersOfTwo = 0, numDataPts = 0, lastAllocationSize=1;
	timeDomain = (doublecomplex*)malloc(sizeof(doublecomplex)*pow(2,powersOfTwo));
	double currentVal;
	while (timeSeriesSignal >> currentVal ){
		numDataPts++;
		if ( numDataPts > lastAllocationSize ){
			powersOfTwo++;
			timeDomain = (doublecomplex*)realloc(timeDomain, sizeof(doublecomplex)*pow(2,powersOfTwo));
			lastAllocationSize = pow(2,powersOfTwo);
		}
		*(timeDomain+numDataPts-1) = doublecomplex(currentVal,0.0);	
	} 
	time_t begin, end;
	if (dftflag==1){
		dftout.open("DFT.dat");
		frequencyDomain = (doublecomplex*)malloc(sizeof(doublecomplex)*lastAllocationSize);
		cout << "Begin DFT " << endl;
		begin = time(NULL);
		DFT(timeDomain, frequencyDomain, lastAllocationSize);
		end = time(NULL);
		if(timeResults==1){
			dftbench << "DFT: " << lastAllocationSize << " " <<  double(end-begin) << " sec" << endl;
		}
		for ( int i=0; i<lastAllocationSize; i++){
			dftout << (FS*double(i)/lastAllocationSize)  << " "<< abs(real(*(frequencyDomain+i))) << endl;
		}
		free(frequencyDomain);
		dftout.close();
	}
	if (fftflag==1){
		fftout.open("FFT.dat");
		frequencyDomainFFT    = (doublecomplex*)malloc(sizeof(doublecomplex)*lastAllocationSize);
		cout << "Begin FFT " << endl;
		begin = time(NULL);
		FFT(timeDomain, frequencyDomainFFT, lastAllocationSize);
		end = time(NULL);
		if(timeResults==1){
			fftbench << "FFT: " << lastAllocationSize << " " << double(end-begin) << " sec" << endl;
		}
		for ( int i=0; i<lastAllocationSize; i++){
			fftout << (FS*double(i)/lastAllocationSize)  << " "<< abs(real(*(frequencyDomainFFT+i))) << endl;
		}
		
		fftout.close();
		free(frequencyDomainFFT);
	}
	if (recfftflag==1){
		recfftout.open("recFFT.dat");
		recfrequencyDomainFFT    = (doublecomplex*)malloc(sizeof(doublecomplex)*lastAllocationSize);
		cout << "Begin recFFT " << endl;
		begin = time(NULL);
		FFT(timeDomain, frequencyDomainFFT, lastAllocationSize);
		end = time(NULL);
		if(timeResults==1){
			recfftbench  << "recFFT: "<< lastAllocationSize << " "  << double(end-begin) << " sec" << endl;
		}
		for ( int i=0; i<lastAllocationSize; i++){
			recfftout << (FS*double(i)/lastAllocationSize)  << " "
						<< abs(real(*(recfrequencyDomainFFT+i)))/lastAllocationSize << endl;
		}
		recfftout.close();
		free(recfrequencyDomainFFT);
	}
	free(timeDomain);
	if(timeResults==1){
		benchmarkout.close();
		dftbench.close();
		fftbench.close();
		recfftbench.close();
	}
	
	return 0;
	
}
