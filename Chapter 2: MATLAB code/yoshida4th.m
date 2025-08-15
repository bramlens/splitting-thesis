function approx = yoshida4th(data)
% This function computes Yoshida's 4th order integrator's approximations
% for a given ODE. It takes the following input:
% -> data: a data struct containing the relevant data of the problem;
% and it produces the following output values:
% -> approx: a dim x N+1 matrix containing the approximations for each
%            variable, row-wise.

% Two vectors containing the q- and p-approximations
q = zeros(0.5*data.dim, data.N+1);
p = zeros(0.5*data.dim, data.N+1);

% Initial data
q(:,1) = data.initial(1:0.5*data.dim);
p(:,1) = data.initial(0.5*data.dim+1:data.dim);

% 4th order coefficients
w1 = -nthroot(2,3) / (2 - nthroot(2,3));
w2 = 1 / (2 - nthroot(2,3));
c  = [0.5*w2, 0.5*(w1 + w2), 0.5*(w1 + w2), 0.5*w2];
d  = [w2, w1, w2];

% Integration
for iter = 1:data.N
    
    % First step: half for position, full for momentum
    Q1 = q(:, iter) + c(1) * data.gradT(p(:, iter)) * data.h;
    P1 = p(:, iter) - d(1) * data.gradV(Q1) * data.h;
    
    % Second step
    Q2 = Q1 + c(2) * data.gradT(P1) * data.h;
    P2 = P1 - d(2) * data.gradV(Q2) * data.h;
    
    % Third step
    Q3 = Q2 + c(3) * data.gradT(P2) * data.h;
    P3 = P2 - d(3) * data.gradV(Q3) * data.h;
    
    q(:, iter + 1) = Q3 + c(4) * data.gradT(P3) * data.h;
    p(:, iter + 1) = P3;
end

% Gathering the approximations
approx(1:0.5*data.dim, :)            = q;
approx(0.5*data.dim + 1:data.dim, :) = p;

end

