function energy = energyCalc(approx, data)
% This function computes the energy associated to each iteration of a
% numerical integrator. It takes the following input values:
% -> approx: a dim x N + 1 matrix containing the approximations (row-wise);
% -> data:   a data struct containing the problem's relevant data values;
% and it produces the following output:
% -> energy: a vector of length data.N + 1 containing the desired energies
%            per iteration.

% Initialising the output vector
energy = sparse(data.N + 1, 1);

% Implementation of the output vector
for index = 1:data.N + 1
    energy(index) = data.H(approx(1:0.5*data.dim, index), approx(0.5*data.dim+1:data.dim,index));
end
end