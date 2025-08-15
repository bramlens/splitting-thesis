function [Q, P, R] = solverMOC(data, int, T)
% This function numerically solves a linear, first order PDE using the
% method of characteristics; it is specifically implemented to solve the
% Liouville equation applied to a Hamiltonian system. It takes the input:
% -> data: a data struct containing all relevant information about the PDE,
%          underlying Hamiltonian system, etc;
% -> int: a string expressing which integrator is to be used to solve the
%         characteristic equations;
% -> T: a time at which we want to compute the MoC solution
%       (standardisation: -2*pi <= T <= 2*pi)
% and it produces the following output:
% -> Q, P, R: matrices representing the solution's values throughout time
%             (row-wise).

% Matrices that will contain the approximations
Q = sparse(length(data.sRange), length(data.tRange));
P = sparse(length(data.sRange), length(data.tRange));
R = sparse(length(data.sRange), length(data.tRange));

% Extracting variables from the data struct for legibility
s = data.sRange;
t = data.tRange; dt = t(2) - t(1);
data.h = dt;

% Method of characteristics
for iter = 1:length(data.sRange) % Picking a point on the initial curve
    
    % Computing the initial values
    s  = data.sRange(iter);
    q0 = data.q0(s);
    p0 = data.p0(s);
    data.initial = [q0; p0];
    
    % Solving the resulting system of ODEs by the chosen integrator for
    % chosen s on the initial curve
    if int == "fe"
        sol = forwardEuler(data);
    elseif int == "se"
        sol = symplecticEuler(data);
    elseif int == "sv"
        sol = stormerVerlet(data);
    elseif int == "y4"
        sol = yoshida4th(data);
    end
    
    % Saving the approximations in the matrices
    Q(iter,:) = sol(1,:);
    P(iter,:) = sol(2,:);
end

% Computation of the solution at selected time T
tPlot = floor(T/dt);
for t_iter = floor((length(data.tRange)/2)+1):length(data.tRange)
    for s_iter = 1:length(data.sRange)
        R(t_iter, s_iter) = data.rho0(Q(s_iter, t_iter - tPlot), P(s_iter, t_iter - tPlot));
    end
end

end