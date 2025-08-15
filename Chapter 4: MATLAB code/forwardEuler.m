function approx = forwardEuler(data)
% This function computes the forward (explicit) Euler approximations for a
% given ODE. It takes the following input:
% -> data: a data struct containing the relevant data of the problem;
% and it produces the following output values:
% -> approx: a dim x N+1 matrix containing the approximations for each
%            variable, row-wise.

% Extracting the function on the ODE's RHS from the data struct to lighten
% notations
f = data.f;

% Matrix containing the approximations, column-wise
% -> rows 1 until dim/2: approximations for the q-variable
% -> rows dim/2 + 1 until dim: approximations for the p-variable
approx      = sparse(data.dim, data.N+1);
approx(:,1) = data.initial;

% Implementing the integrator
for iteration = 2:data.N+1
    approx(:,iteration) = approx(:,iteration-1) + data.h*data.f(approx(:,iteration-1));
end

% For consistency: putting the rows in the right order, if needed
if data.dim == 4
    q2 = approx(3,:);
    p1 = approx(2,:);
    approx(2,:) = q2;
    approx(3,:) = p1;
end

end
