%% Jeffrey Wong | ECE-211 | HW #8

clear
close all
clc

%% Problem 1- Heavy Tail Distributions

% Part a

N = 10^6;
nu = 5;
alpha = 0.544; 

norm = randn(1,N);
tdist = sqrt((nu-2)/nu)*trnd(nu,1,N); % Need to multiply by sqrt((nu-2)/nu) to account for variance not being 1
x = pi*rand(1,N);
cauchy = alpha*tan(x);

gaussianwithinone = mean(abs(norm) < 1); % The fraction of Gaussian values with absolute value less than 1.
tdistwithinone = mean(abs(tdist) < 1); % The fraction of t distribution values with absolute value less than 1.
cauchywithinone = mean(abs(cauchy) < 1); % The fraction of Cauchy values with absolute value less than 1.

disp("The fraction of values in the normal distribution with absolute value less than 1 is " + gaussianwithinone);
disp("The fraction of values in the t distribution with absolute value less than 1 is " + tdistwithinone);
disp("The fraction of values in the Cauchy distribution with absolute value less than 1 is " + cauchywithinone);

figure
hold on
histogram(norm)
xline(-1,"--")
xline(1,"--")
title('Gaussian Distribution')

figure
hold on
histogram(tdist)
xline(-1,"--")
xline(1,"--")
title('Students t-Distribution')

figure
hold on
histogram(cauchy)
xline(-1,"--")
xline(1,"--")
title('Cauchy Distribution')

% Cauchy variables can result in extreme outliers which leads to the
% extremely wide x-axis we see in the graph of the Cauchy distribution


% Part b

norm2 = reshape(norm, 10, 10^5);
tdist2 = reshape(tdist, 10, 10^5);
cauchy2 = reshape(cauchy, 10, 10^5);

normsegmentmeans = mean(norm2, 2); % Returns the mean of each of the 10 segments for the normal distribtuion
tidstsegmentmeans = mean(tdist2, 2);
cauchysegmentmeans = mean(cauchy2, 2);

% It does not make sense to say the Cauchy distribution has mean 0- the
% segment means are not near zero even with 100k values 
% and go as far as +-5 or more in some cases!

%% Problem 2- ARMA and AR Models

% Part a

% 1) ARMA(2,2)
% 2) Innovations filter
% 3) v[n] = x[n] - 1.6x[n-1] + 0.81x[n-2] - 0.4v[n-1] - 0.2v[n-2]
% 4) x[n] - 1.6x[n-1] + 0.81x[n-2] = v[n] + 0.4[n-1] + 0.2[n-2]
% X(z)(1 - 1.6z^-1 + 0.81 z^-2) = V(z)(1 + 0.4z^-1 + 0.2z^-2)
% H(z) = X(z)/V(z) = (1 + 0.4z^-1 + 0.2z^-2)/(1 - 1.6z^-1 + 0.81 z^-2)
% 5) S_x(w) = sigma_v^2 |B(w)|^2/|A(w)|^2
% = 2[(e^j2w + 0.4e^jw + 0.2)(e^-j2w + 0.4e^-jw + 0.2)]/[(e^j2w - 1.6e^jw + 0.81)(e^-2jw - 1.6e^-jw + 0.81)] 

b = [1, 0.4, 0.2];
a = [1, -1.6, 0.81];
[z,p,k] = tf2zpk(b,a);

figure
zplane(z,p)
title('Pole-zero plot of Innovations Filter')

% Part b

N = 10^4;

v = sqrt(2)*randn(1,N); % White noise

x = filter(b,a,v); % Applies innovations filter to get x

r_x = zeros(1,7);

for(m = 0:6)
    r_x(m+1) = dot(x(1:N-m),x(1+m:N))/(N-m);
end

r_neg = r_x(7:-1:2); % Takes the values at m = 1 to 6 and reverses it to get the values from -6 to 1
r_x = [r_neg, r_x]; % Combines the negative and nonnegative parts to get a vector from -6 to +6

figure
stem(-6:6,r_x)
title('r_x(m) for ARMA model')

% Note the indices of r_x are 7 greater than the lag represented (so the
% first row would contain r_x(0) through r_x(6), the 3rd row r_x(-2)
% through r_x(4), etc.

R = toeplitz(r_x(7:13)); % This can be done due to symmetry of r(m) = r(-m)
disp(R)

eigval = eig(R);
isNegative = sum(eigval<=0); % Sum should be zero if every eigenvalue is positive
if(isNegative)
    disp("R is not positive definite- there is at least one negative eigenvalue")
else
    disp("R is positive definite")
end

x_rev = (x(end:-1:1)).'; % Reverses the x vector and transposes it so that latest values are on top of column, as desired

R_est = zeros(7,7,9994); % Want n + 6 to max out at 10^4
% Manually generate the Toeplitz matrix
for(i = 1:N-6)
    x_rev = (x(i+6:-1:i)).'; % Reverses the x vector and transposes it so that latest values are on top of column, as desired
    R_est(:,:,i) = x_rev*x_rev.';
end
R_est = mean(R_est, 3);
R_diff = R - R_est; % MATLAB does not seem to support spectral norm, L2 norm is around 0.01 or less

% Part c

[s_est,w] = pwelch(x, hamming(512),256,512);

figure
plot(w,s_est)
title('Power Spectral Density')
xlabel("Normalized Digital Radian Frequency")

% The peak frequency, w0, is approximately 0.47 (0.45 0.49 0.49)

% The pole angle is arctan(0.412/0.8) = 0.476, which is close to our
% estimated peak frequency

% Part d

[a0,varv] = aryule(x,4);
disp("varv = " + varv)

% varv ranges between 1.9 and 2.1- somewhat close to our actual variance of 2

x0 = filter(1,a0,v);

figure
hold on
stem(1:100,x(1:100))
stem(1:100,x0(1:100))
title('x vs x0')
xlabel("n")

r_x0 = zeros(1,7);

for(m = 0:6)
    r_x0(m+1) = dot(x0(1:N-m),x0(1+m:N))/(N-m);
end

r_neg = r_x0(7:-1:2); 
r_x0 = [r_neg, r_x0]; 

r_diff = abs(r_x - r_x0); % Returns the difference between the two correlation vectors
figure
stem(-6:6,r_diff)
title('Difference between r_x(m) and r_x0(m)') % Difference is not necessarily smaller for -3 <= m <= 3
