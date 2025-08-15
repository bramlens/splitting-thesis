%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% CHAPTER 2: NUMERICAL EXAMPLES (ORDER) %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% BRAM LENS - SPLITTING METHODS %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;

%% Implementation of the exact solution
exact = @(q0, p0, time) [q0*cos(time) + p0*sin(time); -q0*sin(time) + p0*cos(time)];

%% Data struct
data.initial = [2; 3/2];
data.h       = 0.01;
data.dim     = 2;
data.N       = 2000;
data.H       = @(q,p) 0.5*q.^2 + 0.5*p.^2;
data.T       = @(p) 0.5*p.^2;
data.gradT   = @(p) p;
data.V       = @(q) 0.5*q.^2;
data.gradV   = @(q) q;
data.f       = @(x) [x(2); -x(1)];

stepsizes    = logspace(-3, -1, 10);
time_final   = 20; % final time needs to be fixed for good comparisons

%% Gathering the integrators
ints  = {@forwardEuler, @symplecticEuler, @stormerVerlet, @yoshida4th};
label = ["forward Euler", "Symplectic Euler (VT)", "Stormer-Verlet (VTV)", "Yoshida's 4th"];

order = zeros(1, length(label));

%% Estimating the order of the integrators

for iter = 1:length(ints)
    
    % Global error vector
    e_global = zeros(size(stepsizes));
    
    for rep = 1:length(stepsizes)
        
        % Varying step sizes
        data.h = stepsizes(rep);
        
        % Final time computation
        data.N = floor(time_final / data.h);
        Tfin = data.N * data.h;
        
        % Integration
        approx   = ints{iter}(data);
        
        % Numerical resp. exact solution at the final time
        sol_num  = approx(:, end);
        sol_ex   = exact(data.initial(1), data.initial(2), Tfin);
        
        % Computing the global error
        e_global(rep) = norm(sol_num - sol_ex);
    end
    
    % Estimating the order of the integrator
    coef  = polyfit(log(stepsizes), log(e_global), 1);
    order(iter) = coef(1);
    
    subplot(2,2,iter)
    loglog(stepsizes, e_global, "*-", "LineWidth", 1.5, "Color", ...
        [0, 0.4470, 0.7410])
     grid on
    title(sprintf("Estimated order of %s: %.3f", label(iter), order(iter)))
    xlabel("log_{10}(h) with h: step size")
    ylabel("log_{10}(global error)") 
end