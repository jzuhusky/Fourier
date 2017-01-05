# Fourier
This is a side project I like to work on. Implementing my own DFT/FFT code &amp; some other cool stuff

The DFT & FFT routines are inside of the fourier.h file. 

The main program should allow you to run the code from the Unix-command line:

./a.out -m <dfr> -f <samplingFrequency> <inputFile (signal)>

But any code which calls the routines should work. 

Notes:
1. This Main code automatically pads any signal with zeros to make the signal length a power of 2,
which is necessary for the FFT algos to work. 
2. -m = MODE. d -> DFT (slow), r -> recursive FFT (much faster), f -> iterative FFT (Fastest) 
3. -f -> Sampling frequency, this is purely easy plotting. This program will print out 
data in the form "<freq Index> <value> \n". Using GNUplot is an easy way to plot your signals 


