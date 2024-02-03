%% Jeffrey Wong | ECE-210 | HW #7

clear
close all
clc

%% Problem 1

t = 0:0.001:2; % Columns correspond to times (0-2s)
f = 1:50000; % Rows correspond to signals of each frequency
[T,F] = meshgrid(t,f);
X = sin(2.*pi.*T.*F); % Generating our sine waves
x = sum(X); % Summing them together seems to cancel them out?

figure
plot(t,x)
title('Input Signal')

%% Problem 2

Hd2 = Q2Butterworth; % Label variables based on Q# for convenience
output2 = filter(Hd2,x);

Fs = 100000;
N= 2^15;
S2 = fft(output2,N);
S2 = fftshift(abs(S2))/N;
F2 = Fs.*(-N/2:N/2-1)/N;

figure
plot(F2,S2)
title('Fourier Transform of Butterworth-Filtered Signal')
xlabel('Frequency (Hz)')
ylabel('Magnitude')

% Our filtered signal has an extremely fuzzy passband in the fourier domain-
% everything cancelled out almost but not quite due to rounding I think?

%% Problem 3

Hd3 = Q3ChebyshevI; % Label variables based on Q# for convenience
output3 = filter(Hd3,x);

S3 = fft(output3,N);
S3 = fftshift(abs(S3))/N;
F3 = Fs.*(-N/2:N/2-1)/N;

figure
plot(F3,S3)
title('Fourier Transform of Chebyshev I-Filtered Signal')
xlabel('Frequency (Hz)')
ylabel('Magnitude')

%% Problem 4

Hd4 = Q4ChebyshevII; 
output4 = filter(Hd4,x);

S4 = fft(output4,N);
S4 = fftshift(abs(S4))/N;
F4 = Fs.*(-N/2:N/2-1)/N;

figure
plot(F4,S4)
title('Fourier Transform of Chebyshev II-Filtered Signal')
xlabel('Frequency (Hz)')
ylabel('Magnitude')

%% Problem 5

Hd5 = Q5Elliptic; % Label variables based on Q# for convenience
output5 = filter(Hd5,x);

S5 = fft(output5,N);
S5 = fftshift(abs(S5))/N;
F5 = Fs.*(-N/2:N/2-1)/N;

figure
plot(F5,S5)
title('Fourier Transform of Elliptic-Filtered Signal')
xlabel('Frequency (Hz)')
ylabel('Magnitude')