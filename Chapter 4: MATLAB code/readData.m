function data = readData(problem)

if problem == "osc"
    data.initial  = [2; 1.5];
    data.center   = [1; 1];
    data.sigma    = [0.5; 0.5];
    data.a        = @(q,p) p;
    data.b        = @(q,p) -q;
    data.c        = @(q,p,rho) 0;
    data.q0       = @(s) s;
    data.p0       = @(s) 0;
    data.rho0     = @(q,p) exp(-(q - data.center(1)).^2 / ...
                    (data.sigma(1)^2) - (p-data.center(2)).^2/...
                    (data.sigma(2)^2))/(pi*data.sigma(1)*data.sigma(2));
    data.exact    = @(q,p,t) data.rho0(q*cos(t)-p*sin(t),q*sin(t)+p*cos(t));
    data.sRange   = linspace(0,5,50);
    data.tRange   = linspace(-2*pi,2*pi,90);
    data.dim      = 2;
    data.N        = length(data.tRange) - 1;
    data.h        = data.tRange(2) - data.tRange(1);
    data.f        = @(S) [data.a(S(1),S(2)); data.b(S(1),S(2))];
    data.F        = @(q,p) Fbox(q,p);
    data.gradT    = @(p) p;
    data.gradV    = @(q) q;
elseif problem == "pen"
    data.initial = [1.5, 1];
    data.h       = 0.01;
    data.dim     = 2;
    data.H       = @(q,p) 0.5*p.^2 + 1 - cos(q);
    data.T       = @(p) 0.5*p.^2;
    data.gradT   = @(p) p;
    data.gradV   = @(q) sin(q);
    data.f       = @(x) [x(2); -sin(x(1))];
    data.center  = [1; 0.4];
    data.sigma   = [0.5; 0.5];
    data.a       = @(q,p) p;
    data.b       = @(q,p) -sin(q);
    data.c       = @(q,p,rho) 0;
    data.q0       = @(s)  s;
    data.p0       = @(s) 0;
    data.sRange   = linspace(-4,4,200);
    data.tRange   = linspace(-pi,pi,200);
    data.rho0     = @(q,p) exp(-(q - data.center(1)).^2 / ...
                    (data.sigma(1)^2) - (p-data.center(2)).^2/...
                    (data.sigma(2)^2))/(pi*data.sigma(1)*data.sigma(2));
    data.exact    = @(q,p,t) data.rho0(q*cos(t)-p*sin(t),q*sin(t)+p*cos(t));
    data.N        = length(data.tRange) - 1;
    
elseif problem == "kep"
    data.initial = [1-1/2; 0; 0; sqrt((1+1/2)/(1-1/2))];
    data.center  = [data.initial(1); data.initial(2); data.initial(3); data.initial(4)];
    data.sigma   = [1; 1; 1; 1];
    data.rho0    = @(qx, qy, px, py) exp(-(qx-data.center(1)).^2/(data.sigma(1)^2) - (qy-data.center(2)).^2/(data.sigma(2))^2 - (px-data.center(3)).^2/(data.sigma(3))^2 - (py-data.center(4)).^2/(data.sigma(4))^2)/(pi*prod(data.sigma))* (1/(pi^2 * prod(data.sigma)));
    data.q0       = @(s)  s;
    data.p0       = @(s) 0;
    data.H       = @(q,p) 0.5*norm(p,2)^2 - 1/norm(q,2);
    data.T       = @(p) 0.5*norm(p,2);                                 
    data.gradT   = @(p) [p(1); p(2)];                               
    data.V       = @(q) -1/norm(q,2);                               
    data.gradV   = @(q) [q(1); q(2)] / norm(q,2)^3;
    data.f       = @(x) [x(2); -x(1)/(x(1)^2 + x(3)^2)^1.5; x(4); ...
                        -x(3)/(x(1)^2 + x(3)^2)^1.5];
    data.h       = 0.01;                                                   
    data.dim     = 4;                                                      
    data.N       = 10000;
    data.sRange  = linspace(-3,3,500);
    data.tRange  = linspace(-2*pi,2*pi,400);
    
end

end