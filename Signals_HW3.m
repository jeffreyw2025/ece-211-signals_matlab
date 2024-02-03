%% Jeffrey Wong | ECE-211 | HW #3

clear
close all
clc

%% Problem 3

% Part b
x = 1:7;

subplot(2,1,1)
ds = downsamp(x, 2);
x1 = 0:size(ds,2)-1; % size tells us the numver of elements in ds, 0:size-1 contains size elements
p1 = stem(x1, ds, 'ko');
title('Subplot 1: downsamp(x)')

subplot(2,1,2)
us = upsamp(x, 2);
x2 = 0:size(us,2)-1; % us has length 13
p2 = stem(x2, us, 'bo');
title('Subplot 2: upsamp(x)')

% Part c
h = [3, 2, 1, 4];

g = upsamp(h, 2);

c1 = downsamp(conv(g,x),2);
c2 = conv(h,downsamp(x,2));

% Part d
diffd = max(abs(c1-c2)); % Should be zero if the two vectors match, as it means the values are all the same- They do match!
if diffd == 0
    fprintf("The vectors (v2)(g*x) and h*((v2)x) match\n")
else
    fprintf("The vectors (v2)(g*x) and h*((v2)x) do not match\n")
end
   

% Part e
e1 = upsamp(conv(h,x),2);
e2 = conv(g,upsamp(x,2));
diffe = max(abs(e1-e2)); % Should be zero if the two vectors match- They do match!
if diffe == 0
    fprintf("The vectors (^2)(h*x) and g*((^2)x) match\n")
else
    fprintf("The vectors (^2)(h*x) and g*((^2)x) do not match\n")
end

% Function Definitions (Part a)
function y = downsamp(x, M)
    y = x(1:M:end);
end

function z = upsamp(x, M)
    temp = zeros(M,size(x,2)); % The matrix we create will be "unwrapped" later- 1st row holds values in x, rest are the intermediate zeroes
    temp(1,:) = x; % Changes first row of matrix to x, the original signal
    z = temp(1:end-M+1); % "unwraps" the matrix we have created, column then row. Also truncates end zeroes.
end
