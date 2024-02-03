%% Jeffrey Wong | ECE-210 | HW #4

clear
close all
clc

%% Problem 1

% Part d- Testing
samp0 = [2, 3, 4; 0, 2, 2; 0, 1, 7]; % Test on a simple real-valued case, for personal use
samp1 = [1 + 1j, 1 + 2j; 2 + 1j, 2 + 2j; 3 + 1j, 3 + 2j]; % Generates 2 complex length 3 vectors (Case 1 for testing)
samp2 = [1 + 1j, 1 + 2j, 1 + 3j; 2 + 1j, 2 + 2j, 2 + 3j; 3 + 1j, 3 + 2j, 3 + 3j]; % Generates 3 complex length 3 vectors (Case 2 for testing)

test0 = gramSchmidt(samp0);
check0 = isOrthonormal(test0);
est0 = orthoProj([1; 2; 3], test0);
err0 = norm([1; 2; 3] - est0);

test1 = rand(size(samp1), "like", samp1); % The "like" keyword generates an object of the form of samp1 (See above)
test1 = gramSchmidt(test1);
check1 = isOrthonormal(test1);
est1 = orthoProj([3+1j; 1-2j; 4], test1);
err1 = norm([3+1j; 1-2j; 4] - est1); % Error seems to be hard to pin down for complex-valued functions- not sure why

test2 = rand(size(samp2), "like", samp2);
test2 = gramSchmidt(test2);
check2 = isOrthonormal(test2);
est2 = orthoProj([3+1j; 1-2j; 4], test2);
err2 = norm([3+1j; 1-2j; 4] - est2); % Error seems to be hard to pin down for complex-valued functions- not sure why
% Intuition suggests err2 should tend to be less than err1 because the
% third vector allows for a more complete basis and thus a more accurate
% estimation- This doesn't seem to hold that well

% Part e- Plotting
x = 0:0.001:2;
mu = 0:0.5:2;
[xs,mus] = ndgrid(x,mu);
gaussians = exp(-(xs-mus).^2)/sqrt(2*pi);

figure
hold on
plot(x, sin(pi*x),'DisplayName',"sin(x)")
for(i = 1:length(mu))
    legendString = ['\mu = ',num2str(mu(i))];
    plot(x, gaussians(:,i),'DisplayName',legendString)
end
legend
grid on
title("sin(x) and Gaussians")
xlabel("x")
ylabel("f(x)")

basis = gramSchmidt(gaussians);

sin_est = orthoProj(transpose(sin(x*pi)),basis);

subplot(2,1,1)
hold on
plot(x, sin(x*pi),'DisplayName',"sin(x)")
plot(x, sin_est,'DisplayName',"Estimate")
legend
grid on
title("Sinusoid and Estimate")
xlabel("x")
ylabel("f(x)")

subplot(2,1,2)
hold on
for(i = 1:length(mu))
    legendString = ['\phi ', num2str(i)];
    plot(x, basis(:,i),'DisplayName',legendString)
end
legend
grid on
title("Basis Functions")
xlabel("x")
ylabel("f(x)")

% Function Definitions

% Part a- Gram-Schmidt
function phis = gramSchmidt(X)
    phis = X;
    for i = 1:size(X,2)
        for j = 1:i-1
            phis(:,i) = phis(:,i) - (dot(phis(:,j), X(:,i)) * phis(:,j)); % Represents inner product times phi vectors
        end
        phis(:,i) = phis(:,i)/norm(phis(:,i)); % Normalization step
    end
end

% Part b- Determining orthonormality
function ortho = isOrthonormal(X)
    ortho = 1;
    for i = 1:size(X,2)
        for j = i:size(X,2) % If we check i and j we don't need to check j and i
            % Inner products must be zero for vectors that are different,
            % and 1 for vectors with themselves
            if(i == j & abs(dot(X(:,i), X(:,j))-1) > 10^-5)
                ortho = 0; % Will automatically return 0 (false) if any deviancy occurs
                return
            elseif(i ~= j & abs(dot(X(:,i), X(:,j))) > 10^-5)
                ortho = 0;
                return
            end
        end
    end
end

% Part c- Projection!

% Note: Error seems to be hard to pin down for complex-valued functions- not sure why
function estimate = orthoProj(v, X)
    estimate = zeros(length(v),1);
    for i = 1:size(X,2) % Finds the component of v along each vector in X
        estimate = estimate + (dot(v, X(:,i)) * X(:,i));
    end
end