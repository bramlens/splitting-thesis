function y = f_rhs(state, data)
% This function computes the right-hand side of the ODE associated to the
% three body problem. It takes the following input values:
% -> state: an 18x1 vector of positions (1:9) and velocities (10:18);
% -> data:  a data struct containing all the data values relevant to the
%           3BP;
% and it produces the following output values:
% -> y:     an 18x1 vector representing the right-hand side of the ODE
%           associated to the 3BP.

% For convenience: abbreviating needed data variables' names
m     = data.m;
accel = data.accel;

% Implementing the desired 18x1 vector
y = zeros(18,1);

y(1:3)   = state(10:12);
y(4:6)   = state(13:15);
y(7:9)   = state(16:18);
y(10:12) = accel(state(1:3), state(4:6), state(7:9), m(1), m(2), m(3));
y(13:15) = accel(state(4:6), state(1:3), state(7:9), m(2), m(1), m(3));
y(16:18) = accel(state(7:9), state(1:3), state(4:6), m(3), m(1), m(2));

end

