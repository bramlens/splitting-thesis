function sol = exact(problem)
% This function computes the exact solution to a selected problem. It takes
% the following input values:
% -> problem: a string that determines which problem is to be considered,
%             its possible values are "osc", "kep", and "rtb";
% and it produces the following output values:
% -> sol: a vector containing the values of the problem's exact solution.

% Collecting the relevant data
data = readData(problem);
sol  = sparse(1e6, 2);

if problem == "osc"
    
    % Time interval (domain) of the problem
    time_interval = linspace(0,2*pi, 1e6);
    
    % Computing the solution's x- and y-coordinates
    sol(:,1) = data.initial(1)*cos(time_interval) + data.initial(2)*sin(time_interval);
    sol(:,2) = data.initial(2)*cos(time_interval) - data.initial(1)*sin(time_interval);
    
elseif problem == "kep"
    
    % Initial energy H0
    H0  = data.H(data.initial(1:2), data.initial(3:4));
    
    % Initial angular momentum L0
    L0  = data.initial(1)*data.initial(4) - data.initial(2)*data.initial(3);
    
    % Eccentricity of the orbit
    ecc = sqrt(1 + 2*H0*L0^2);
    
    % Computing the exact orbit in polar coordinates
    angle  = linspace(0,2*pi, 1e6);
    radius = L0^2 ./ (1 + ecc * cos(angle));
    
    % Transformation of the orbit into cartesian coordinates
    [x, y] = pol2cart(angle, radius);
    
    % Implementing a sparse matrix that will contain the solution's
    % x-coordinates (column 1) and y-coordinates (column 2)
    sol      = sparse(1e6, 2);
    sol(:,1) = x;
    sol(:,2) = y;
end

end

