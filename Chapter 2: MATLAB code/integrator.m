function approx = integrator(data, int)
% This function computes the approximations to a given ODE for a chosen
% numerical integrator. It takes the following input values:
% -> data: a data struct containing the important data values, see the
%          function "readData.m" for more information;
% -> int : a string expressing the integrator to be considered;
% and it produces the following output values:
% -> approx: a dim x N + 1 matrix containing the approximations to the
%            problem's solution.

if int == "fe"
    approx = forwardEuler(data);
elseif int == "se"
    approx = symplecticEuler(data);
elseif int == "sv"
    approx = stormerVerlet(data);
elseif int == "y4"
    approx = yoshida4th(data);
end

end

