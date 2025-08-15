function approx = symplecticEuler(data)
% This function computes the symplectic Euler (VT) approximations for a
% given ODE. It takes the following input:
% -> data: a data struct containing the relevant data of the problem;
% and it produces the following output values:
% -> approx: a dim x N+1 matrix containing the approximations for each
%            variable, row-wise.

% Matrix containing the approximations, column-wise
% -> rows 1 until dim/2: approximations for the q-variable
% -> rows dim/2 + 1 until dim: approximations for the p-variable
approx      = sparse(data.dim, data.N + 1);
approx(:,1) = data.initial;

% Implementing the integrator
for iter = 2:data.N + 1
    approx(1:(0.5*data.dim), iter)            = approx(1:(0.5*data.dim), iter-1) + data.h * data.gradT(approx((0.5*data.dim + 1):data.dim, iter-1));
    approx((0.5*data.dim + 1):data.dim, iter) = approx((0.5*data.dim + 1):data.dim, iter-1) - data.h * data.gradV(approx(1:(0.5*data.dim), iter));
end

end

