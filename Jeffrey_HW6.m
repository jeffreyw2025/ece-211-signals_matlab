%% Jeffrey Wong | ECE-210 | HW #6

clear
close all
clc

%% Problem 1

% Part a
b = [3/7, 2/3, 1/2];
a = [1/2, 1/3, 2];
[z,p,k] = tf2zp(b,a);

% Part b
figure
zplane(z,p)

% Part c
[h,t] = impz(b,a,50);
n = 0:49;

figure
stem(n,h,'bo')
title("Impulse Response Plot")

% Part d
X = (-0.5).^n;
Y1 = filter(b,a,X);

figure
subplot(2,1,1)
stem(n,X)
title("X")

subplot(2,1,2)
stem(n,Y1)
title("Filter via Function on X")

% Part e
Y2 = conv(X,h);
Y2 = Y2(1:50); % We only want to take the first 50 values to match the filter

figure
subplot(2,1,1)
stem(n,X);
title("X")

subplot(2,1,2)
stem(n,Y2)
title("Filter via Convolution with X")


%% Problem 2

% The difference equation is given by y[n] =
% y[n-1]+y[n-2]+delta[n-1]+delta[n] (Due to ICs)

fib = zeros(1,100); % Preallocation for a bit of speed
fib(1) = 1;
fib(2) = 1;
for(i = 3:100)
    fib(i) = fib(i-1) + fib(i-2);
end

figure
semilogy(0:99, fib, 'b.')
title("Fibonacci Plot")


% The difference equation translates to Y(z) = z^-1Y(z) + z^-2Y(z) + z^-1 + 1, or
% Y(z)/X(z) = H(z) = z^-1+1/(1-z^-1-z^-2)
% The transfer function of the system is given by (z^2+z)/(z^2 - z - 1)
% The zeros would then be at 0 and -1, and the poles are (1+-sqrt(5))/2

b2 = [1 1 0];
a2 = [1 -1 -1];
[z2,p2,k2] = tf2zp(b2,a2);

figure
zplane(z2,p2)

% The system is causal, as the output is only dependent on prior output
% terms, as indicated by the differece equations. Since our system is
% causal, our system cannot be stable. First of all, because there is a
% pole outside the unit circle, our system cannot be both stable and
% causal, and we already know our system is causal. Additionally, if we
% take the L1 norm of the impulse response, we find that it is infinite as
% our terms are monotonically increasing for n > 1.

