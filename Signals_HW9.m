%% Jeffrey Wong | ECE-211 | HW #9

clear
close all
clc

%% Problem 1

% Euler Method

% x_de[n+1] = x_de[n] + dt(Ax_de[n] + Bu_de[n]) = (Adt + I)x_de[n] + dtBu_de[n]
% y_de[n] = Cx_de[n] + Du_de[n] % Release the bogus!

% A_de = (I + dt*A)
% B_de = (dt*B)
% C_de = C % How we compute our state transitions should not affect how we create our outputs
% D_de = D

% Midpoint Method

% x_dm[n+1] = x_dm[n] + dt(A(x_dm[n] + dt/2(Ax_dm[n]+B[I_m, 0_m; 0_m, 0_m]u_dm[n]) + B[0_m, 0_m; 0_m, I_m]u_dm[n])
% = x_dm + dt*Ax_dm[n] + dt^2/2*A^2x_dm[n] + dt^2/2*AB[I_m, 0_m; 0_m, 0_mu_dm[n] + dt*B[0_m, 0_m; 0_m, I_m]u_dm[n])
% y_dm[n] = Cx_dm[n] + Du_dm[n] 

% A_dm = (I + dt*A + dt^2/2*A^2)
% B_dm = (dt^2/2*AB[I_m, 0_m; 0_m, 0_m] + dt*B[0_m, 0_m; 0_m, I_m])
% C_dm = C
% D_dm = D

%% Problem 2

% The formulas for A_de and A_dm resemble the "MacLaurin series" for e^At at Adt up
% to the first-order and second-order terms, respectively- makes sense as an approximation.

%% Problem 3

% Part a

E = [3, 1; 2, 1];
A_0 = [-0.3, 0.4; -0.4, -0.3];
eigvalA_0 = eig(A_0);
A = E*A_0*inv(E);
eigvalA = eig(A);
sz = length(A); % Used for future reference

% Eigenvalues are -0.3 +- 0.4j for both A_0 and A

% Part b

s = sym('s');
F = inv(s*eye(sz)-A);
f = ilaplace(F);
trmatrix = matlabFunction(f); % e^At represents a state transition matrix

% Part c

nmax = 100; % We run the matrix until n = 100 (Corresponding to t = 10 s)
fs = 10; % Sample frequency

A_d = trmatrix(1/fs); % A_d = e^Adt = trmatrix(dt) = trmatrix(1/fs)

xd0 = [2; 1];
xdn = zeros(sz, nmax+1);

for(n = 0:100)
    xdn(:,n+1) = A_d^n*xd0;
end

% Part d

A_de = eye(sz) + A/fs;
xden = zeros(sz, nmax+1);

for(n = 0:100)
   xden(:,n+1) = A_de^n*xd0;
end

A_dm = eye(sz) + A/fs + (A^2)/(2*fs^2);
xdmn = zeros(sz, nmax+1);

for(n = 0:100)
   xdmn(:,n+1) = A_dm^n*xd0;
end

% Part e

figure
legend
hold on
plot3(xdn(1,:),xdn(2,:),0:0.1:10, 'DisplayName', "xd")
plot3(xden(1,:),xden(2,:),0:0.1:10, 'DisplayName', "xde")
plot3(xdmn(1,:),xdmn(2,:),0:0.1:10, 'DisplayName', "xdm")
xlabel('x_{d1}(t)'),ylabel('x_{d2}(t)'),zlabel('t')
title('State Space Trajectories')

% Part f

errde = max(abs(xdn - xden)); % errde and errdm will be used for comparison in part g
errdm = max(abs(xdn - xdmn));

maxerrde = max(errde);
maxerrdm = max(errdm);

% Part g

% The midpoint method seems to be substantially better than Euler's
% method, particularly over the 2.5 - 4.1 second range. x_de does not seem
% to curve inward as sharply as x_d or x_dm, as x_dm can "turn" twice per 
% dt as opposed to once, which causes x_dm to travel 
% further out until it turns around at 3.3 seconds (the same time the other
% curves turn around).