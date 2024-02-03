%% Jeffrey Wong | ECE-211 | HW #5

clear
close all
clc

%% Problem 1

% Note: MATLAB outputs the form denoted by [z,p,k] (zp form) or [b,a]
% (transfer function form) by calling a different version of the cheby1
% function based on the output (a version of function overloading)

% Part a

fcrit = [0.5, 0.6]; % Nyq. Bandwith is 20 kHz, normalized values are 10k/20k and 12k/20k

[z,p,k] = cheby1(3,2,fcrit);

zplane(z,p)

[b,a] = zp2tf(z,p,k);

% Part b

f = linspace(0,2e4,10000);

h = freqz(b,a,pi*f/2e4);

% Part c
mag = 20*log10(abs(h));
pha = unwrap(angle(h)*180/pi);

figure
subplot(2,1,1)
plot(f/1e3,mag);
xlabel('Frequency (kHz)');
ylabel('Gain (db)');
title('Magnitude Response');
grid on;
axis([0 20 -40 2]);

subplot(2,1,2)
plot(f/1e3,pha);
xlabel('Frequency (kHz)');
ylabel('Phase (deg)');
title('Phase Response');
grid on;

% Part d

% transitionPoint finds where a vector crosses from positive to negative value or vice versa.
% Calling transitionPoint on mag+30 allows us to find where magnitude crosses the -30dB threshold
crossings = transitionPoint(mag+30);

% Finds the highest frequency from DC to 10 kHz where the gain crosses -30 dB
% The indices correspond to values from 0 to 10 kHz spaced by ~2 Hz, so we
% have to multiply by 2 to get the actual frequency
f1 = find(crossings(1:5000),1,"last")*2; 

% Finds the lowest frequency from 12 kHz to 20 kHz where the gain crosses -30 dB
% The indices correspond to values from 12 kHz to 20 kHz spaced by ~2 Hz, so we
% have to multiply by 2 then add 12k to get the actual frequency
f2 = 12000+find(crossings(6001:end),1,"first")*2; 

fprintf("The gain for the filter falls below -30 dB at about %f kHz and %f kHz.\n", f1/1e3, f2/1e3);


%% Problem 4

Kmax = 4;
K0 = 4-2*sqrt(2);
K1 = K0;
K2 = 0.5*K0 + 0.5*Kmax;
K3 = 0.2*K0 + 0.8*Kmax;

RC = 1e-5; % RC constant is 10nF*1kOhm = 10 us = 10^-5 s

% Part a

flin = linspace(0,1e6,1e4);

h1lin = freqs([K1/RC, 0], [1, (4-K1)/RC, 2/(RC)^2], 2*pi*flin); % Use 2*pi to convert from freq. to angular freq.
h2lin = freqs([K2/RC, 0], [1, (4-K2)/RC, 2/(RC)^2], 2*pi*flin);
h3lin = freqs([K3/RC, 0], [1, (4-K3)/RC, 2/(RC)^2], 2*pi*flin);

mag1lin = 20*log10(abs(h1lin));
pha1lin = unwrap(angle(h1lin)*180/pi);
mag2lin = 20*log10(abs(h2lin));
pha2lin = unwrap(angle(h2lin)*180/pi);
mag3lin = 20*log10(abs(h3lin));
pha3lin = unwrap(angle(h3lin)*180/pi);

figure
subplot(2,1,1)
plot(flin/1e3,mag1lin);
xlabel('Frequency (kHz)');
ylabel('Gain (db)');
title('Linear Plot of Magnitude Response for K = K1');
axis([0 1000 mag1lin(end) max(mag1lin)+3]); % mag1(end) gets the magnitude at 1MHz, and max(mag1) gives the maximum gain
grid on;

subplot(2,1,2)
plot(flin/1e3,pha1lin);
xlabel('Frequency (kHz)');
ylabel('Phase (deg)');
title('Linear Plot of Phase Response for K = K1');
grid on;

figure
subplot(2,1,1)
plot(flin/1e3,mag2lin);
xlabel('Frequency (kHz)');
ylabel('Gain (db)');
title('Linear Plot of Magnitude Response for K = K2');
axis([0 1000 mag2lin(end) max(mag2lin)+3]);
grid on;

subplot(2,1,2)
plot(flin/1e3,pha2lin);
xlabel('Frequency (kHz)');
ylabel('Phase (deg)');
title('Linear Plot of Phase Response for K = K2');
grid on;

figure
subplot(2,1,1)
plot(flin/1e3,mag3lin);
xlabel('Frequency (kHz)');
ylabel('Gain (db)');
title('Linear Plot of Magnitude Response for K = K3');
axis([0 1000 mag3lin(end) max(mag3lin)+3]);
grid on;

subplot(2,1,2)
plot(flin/1e3,pha3lin);
xlabel('Frequency (kHz)');
ylabel('Phase (deg)');
title('Linear Plot of Phase Response for K = K3');
grid on;

% Part b

flog = logspace(3,6,1e4);

h1log = freqs([K1/RC, 0], [1, (4-K1)/RC, 2/(RC)^2], 2*pi*flog);
h2log = freqs([K2/RC, 0], [1, (4-K2)/RC, 2/(RC)^2], 2*pi*flog);
h3log = freqs([K3/RC, 0], [1, (4-K3)/RC, 2/(RC)^2], 2*pi*flog);

mag1log = 20*log10(abs(h1log));
pha1log = unwrap(angle(h1log)*180/pi);
mag2log = 20*log10(abs(h2log));
pha2log = unwrap(angle(h2log)*180/pi);
mag3log = 20*log10(abs(h3log));
pha3log = unwrap(angle(h3log)*180/pi);

figure
subplot(2,1,1)
semilogx(flog/1e3,mag1log);
xlabel('Frequency (kHz)');
ylabel('Gain (db)');
title('Bode Plot of Magnitude Response for K = K1');
axis([1 1000 mag1log(end) max(mag1log)+3]); % mag1(end) gets the magnitude at 1MHz, and max(mag1) gives the maximum gain
grid on;

subplot(2,1,2)
semilogx(flog/1e3,pha1log);
xlabel('Frequency (kHz)');
ylabel('Phase (deg)');
title('Bode Plot of Phase Response for K = K1');
grid on;

figure
subplot(2,1,1)
semilogx(flog/1e3,mag2log);
xlabel('Frequency (kHz)');
ylabel('Gain (db)');
title('Bode Plot of Magnitude Response for K = K2');
axis([1 1000 mag2log(end) max(mag2log)+3]);
grid on;

subplot(2,1,2)
semilogx(flog/1e3,pha2log);
xlabel('Frequency (kHz)');
ylabel('Phase (deg)');
title('Bode Plot of Phase Response for K = K2');
grid on;

figure
subplot(2,1,1)
semilogx(flog/1e3,mag3log);
xlabel('Frequency (kHz)');
ylabel('Gain (db)');
title('Bode Plot of Magnitude Response for K = K3');
axis([1 1000 mag3log(end) max(mag3log)+3]);
grid on;

subplot(2,1,2)
semilogx(flog/1e3,pha3log);
xlabel('Frequency (kHz)');
ylabel('Phase (deg)');
title('Bode Plot of Phase Response for K = K3');
grid on;

% Part c

wn = sqrt(2)/RC; % Natural frequency does not depend on K
q1 = sqrt(2)/(4-K1);
q2 = sqrt(2)/(4-K2);
q3 = sqrt(2)/(4-K3);

fprintf("The values of wn and Q for K = K1 are %f and %f, respectively\n", wn, q1);
fprintf("The values of wn and Q for K = K2 are %f and %f, respectively\n", wn, q2);
fprintf("The values of wn and Q for K = K3 are %f and %f, respectively\n", wn, q3);

% Part d

[P1,I1] = max(mag1log); % We want the index where the maximum value occurs for below and the peak to find the 3dB points
[P2,I2] = max(mag2log);
[P3,I3] = max(mag3log);

% Values in mag1log span exponentially from 10^3 to 10^6, meaning they increase by a factor of 10 after 1/3 of the matrix
% Use the equation 1000*10^(Index-1)/3333 to convert from index to frequency to find peak freq, flo, and fhi
f01 = 1000*10^((I1-1)/3333); 
f02 = 1000*10^((I2-1)/3333); 
f03 = 1000*10^((I3-1)/3333); 

crossings1 = transitionPoint(mag1log-P1+3); % Adapting code from Q1
crossings2 = transitionPoint(mag2log-P2+3);
crossings3 = transitionPoint(mag3log-P3+3);

loclo1 = find(crossings1(1:I1),1,"last");
lochi1 = find(crossings1(I1+1:end),1,"first");
loclo2 = find(crossings2(1:I2),1,"last");
lochi2 = find(crossings2(I2+1:end),1,"first");
loclo3 = find(crossings3(1:I3),1,"last");
lochi3 = find(crossings3(I3+1:end),1,"first");

flo1 = 1000*10^((loclo1-1)/3333);
fhi1 = 1000*10^((lochi1+I1)/3333);
flo2 = 1000*10^((loclo2-1)/3333);
fhi2 = 1000*10^((lochi2+I2)/3333); 
flo3 = 1000*10^((loclo3-1)/3333);
fhi3 = 1000*10^((lochi3+I3)/3333); 

fprintf("When K = K1, f0 = %f, flo = %f, and fhi = %f\n", f01, flo1, fhi1);
fprintf("When K = K2, f0 = %f, flo = %f, and fhi = %f\n", f02, flo2, fhi2);
fprintf("When K = K3, f0 = %f, flo = %f, and fhi = %f\n", f03, flo3, fhi3);

% Part e

% Estimate based off of calculated peak frequencies
wEst1 = 2*pi*f01;
wEst2 = 2*pi*f02;
wEst3 = 2*pi*f03;
% The estimated frequency is about 1.4139e5 rad/s, very close to wn =
% sqrt(2)e5- likely difference is due to rounding error with the logspace

% Part f

% Bandwidths
delf1 = fhi1 - flo1;
delf2 = fhi2 - flo2;
delf3 = fhi3 - flo3;
% Center Frequencies
fctr1 = sqrt(fhi1*flo1);
fctr2 = sqrt(fhi2*flo2);
fctr3 = sqrt(fhi3*flo3);

% Part g

qEst1 = fctr1/delf1; % Expected 0.5, estimate 0.5009- close enough
qEst2 = fctr2/delf2; % Expected 1.0, estimate 1.0016- close enough
qEst3 = fctr3/delf3; % Expected 2.5, estimate 2.5009- close enough