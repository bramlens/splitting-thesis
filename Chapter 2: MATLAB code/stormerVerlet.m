function approx = stormerVerlet(data)
% This function computes the Stormer-Verlet (VTV) approximations for a
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
for iter = 1:data.N
    
    % Half step
    p_half = approx((0.5*data.dim + 1):data.dim, iter) - 0.5*data.h*data.gradV(approx(1:0.5*data.dim, iter));
    
    % Updating position q
    approx(1:0.5*data.dim, iter + 1) = approx(1:0.5*data.dim, iter) + data.h*data.gradT(p_half);
    
    % Updating momenta p
    approx((0.5*data.dim+1):data.dim, iter + 1) = p_half - 0.5*data.h*data.gradV(approx(1:0.5*data.dim, iter + 1));
    
end

end

