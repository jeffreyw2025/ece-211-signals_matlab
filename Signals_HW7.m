%% Jeffrey Wong | ECE-210 | HW # 7

clear
close all
clc

%% Problem 1b + 2- Function Call and Analysis

[A1,S1] = generate_signal(100, 200, 20, [0, -2, -4], 10);

R1 = (A1*A1.')/200;

[U1,sval1,V1] = svd(A1);
sval1 = diag(sval1);

[eigvec1,eigval01] = eig(R1);
[eigval1,idx] = sort(diag(eigval01),'descend');
eigvec1 = eigvec1(:,idx);

figure
stem(1:length(sval1),sval1)
title('Singular Values (M = 100, N = 200)')
singratio1 = sval1(3)/sval1(4);
disp("The ratio of the third singular value to the fourth when N=200 is " + singratio1);

figure
stem(1:length(eigval1),eigval1)
title('Eigenvalues (M = 100, N = 200)')
eigratio1 = eigval1(3)/eigval1(4);
disp("The ratio of the third eigenvalue to the fourth when N=200 is " + eigratio1);

UL1 = U1(:,1:3);

PN1 = eye(length(UL1)) - UL1 * UL1.';

Rinv = inv(R1);

SMUSIC_Source1 = 1./(S1.' * PN1 * S1);
% We only want the diagonal entries because we only want to consider signals multiplied with themselves- ignore "cross-terms"
SMUSIC_Source1 = diag(SMUSIC_Source1);

SMVDR_Source1 = 1./(S1.' * Rinv * S1);
SMVDR_Source1 = diag(SMVDR_Source1);

Srand1 = zeros(100,20);

for(i = 1:20)
    Srand1(:,i) = generate_sourcevector(100,20);
end

SMUSIC_Rand1 = 1./(Srand1.' * PN1 * Srand1);
SMUSIC_Rand1 = diag(SMUSIC_Rand1);

SMVDR_Rand1 = 1./(Srand1.' * Rinv * Srand1);
SMVDR_Rand1 = diag(SMVDR_Rand1);

% Distribution of random test vector spectra
max_rand_MUSIC1 = max(SMUSIC_Rand1);
mean_rand_MUSIC1 = sum(SMUSIC_Rand1)/length(SMUSIC_Rand1);
median_rand_MUSIC1 = median(SMUSIC_Rand1);

specratio_MUSIC1 = SMUSIC_Source1(3)/max_rand_MUSIC1; % Indicates how large the values are for the given vs random source vectors

max_rand_MVDR1 = max(SMVDR_Rand1);
mean_rand_MVDR1 = sum(SMVDR_Rand1)/length(SMVDR_Rand1);
median_rand_MVDR1 = median(SMVDR_Rand1);

specratio_MVDR1 = SMVDR_Source1(3)/max_rand_MVDR1;

% MUSIC seems to be more effective at identifying the correct source
% vectors- The ratio of the smallest spectrum for the actual source vectors
% to the maximum of the spectrum of the randomly generated vectors is
% consistently higher using MUSIC, and it gives higher values for the
% spectra in general

%% Problem 3- Try Again

[A2,S2] = generate_signal(100, 50, 20, [0, -2, -4], 10);

[U2,sval2,V2] = svd(A2);
sval2 = diag(sval2);

figure
stem(1:length(sval2),sval2)
title('Singular Values (M = 100, N = 50)')
singratio2 = sval2(3)/sval2(4);
disp("The ratio of the third singular value to the fourth when N=50 is " + singratio2);

UL2 = U2(:,1:3);

PN2 = eye(length(UL2)) - UL2 * UL2.';

SMUSIC_Source2 = 1./(S2.' * PN2 * S2);
SMUSIC_Source2 = diag(SMUSIC_Source2);

Srand2 = zeros(100,20);

for(i = 1:20)
    Srand2(:,i) = generate_sourcevector(100,20);
end

SMUSIC_Rand2 = 1./(Srand2.' * PN1 * Srand2);
SMUSIC_Rand2 = diag(SMUSIC_Rand2);

% Distribution of random test vector spectra
max_rand_MUSIC2 = max(SMUSIC_Rand2);
mean_rand_MUSIC2 = sum(SMUSIC_Rand2)/length(SMUSIC_Rand2);
median_rand_MUSIC2 = median(SMUSIC_Rand2);

specratio_MUSIC2 = SMUSIC_Source2(3)/max_rand_MUSIC2;

% MUSIC does work when N is reduced to 50, but it performs substantially
% worse, as the spectra values for the actual source vectors are
% substantially lower while the random vectors have similar spectras to
% before of about 1, meaning the actual source vectors are less distingushable than random ones.

%% Problem 4- One Last Thing

STS = S1.'*S1;
disp(STS)

% The entries STS_ij of STS denote the number of components Si and Sj have
% in common divided by the total number of nonzero components, which acts
% as a measure of how correlated the source vectors are.

%% Problem 1a- Signal Generation Functions

function l = generate_sourcevector(M,K)
    components = randperm(M,K);
    [ind,comp] = meshgrid(components,1:M);
    l = sum(ind == comp,2).*(1/sqrt(K));
end

function [A,S] = generate_signal(M, N, K, PdB, PndB)
    L = length(PdB);
    S = zeros(M,L);
    for(i = 1:L) % Can't call randperm to generate all the source vectors at once, no vectorization trick for that
        S(:,i) = generate_sourcevector(M,K);
    end
    
    B = zeros(L,N);
    for(j = 1:L) % Different variances for each source vector means a for loop is needed
        B(j,:) = 10^(PdB(j)/20)*randn(N,1); % Variance given by 10^(PndB/10), sqrt can be incorporated into exponent.
    end

    V = 10^(PndB/20)*randn(M,N); % Variance given by 10^(PndB/10), sqrt can be incorporated into exponent.
    A = S * B + V/sqrt(M); 
end