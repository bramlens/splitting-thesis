function data = readData(problem)
% This function collects the relevant data for a specific problem and
% collects it in a single data struct. It takes the following input values:
% -> problem: a string expressing which problem is to be considered, the
%             uncoupled harmonic oscillator "osc", the Kepler problem
%             "kep", or the three body problem "3bp";
% and it produces the following output values:
% -> data: a data struct containing the relevant data for the problem to be
% considered.

if problem == "osc"
    
    % Was a valid problem entered?
    data.problemError_flag = 0;
    
    % Is the problem integrable?
    data.integrable = true;
    
    % Data collection
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
    
elseif problem == "pen"
    
    % Was a valid problem entered?
    data.problemError_flag = 0;
    
    % Is the problem integrable?
    data.integrable = false;
    
    data.initial = [1.5,1];
    data.h       = 0.01;
    data.dim     = 2;
    data.N       = 2000;
    data.H       = @(q,p) 0.5*p.^2 + 1 - cos(q);
    data.T       = @(p) 0.5*p.^2;
    data.gradT   = @(p) p;
    data.gradV   = @(q) sin(q);
    data.f       = @(x) [x(2); -sin(x(1))];
    
elseif problem == "kep"
    
    % Was a valid problem entered?
    data.problemError_flag = 0;
    
    % Is the problem integrable?
    data.integrable = true;
    
    % Data collection
    data.initial = [1-1/2; 0; 0; sqrt((1+1/2)/(1-1/2))];
    data.h       = 0.01;
    data.dim     = 4;
    data.N       = 10000;
    data.H       = @(q,p) 0.5*norm(p,2)^2 - 1/norm(q,2);
    data.T       = @(p) 0.5*norm(p,2);
    data.gradT   = @(p) [p(1); p(2)];
    data.V       = @(q) -1/norm(q,2);
    data.gradV   = @(q) [q(1); q(2)] / norm(q,2)^3;
    data.L       = @(q,p) q(1)*p(2) - q(2)*p(1);
    data.f       = @(x) [x(2); -x(1)/(x(1)^2 + x(3)^2)^1.5; x(4); ...
                        -x(3)/(x(1)^2 + x(3)^2)^1.5];
                    
elseif problem == "3bp"
    
    % Was a valid problem entered?
    data.problemError_flag = 0;
    
    % Is the problem integrable?
    data.integrable = false;
    
    % Constants
    data.G    = 1;
    data.m    = [1, 1, 1];
    data.N    = 20000;
    data.dim  = 18;
    data.h    = 0.01;
    
    % Initial data
    data.initial = [ 0.97000436; -0.24308753; 0; ...
                    -0.97000436;  0.24308753; 0; ...
                     0;           0;          0; ...
                     0.466203685; 0.43236573; 0; ...
                     0.466203685; 0.43236573; 0; ...
                    -0.93240737; -0.86473146; 0];
    
    % Function handles
    data.accel = @(r1, r2, r3, m1, m2, m3) ...
             data.G*(m2*(r2-r1)/norm(r2-r1)^3+m3*(r3-r1)/norm(r3-r1)^3);
    data.f     = @(state) f_rhs(state, data);
    data.H     = @(q,p) 1/(2*data.m(1)) *norm(p(1:3),2)^2+...
                        1/(2*data.m(2)) *norm(p(4:6),2)^2+...
                        1/(2*data.m(3)) *norm(p(7:9),2)^2-...
                        data.G*(data.m(1)*data.m(2)/norm(q(1:3)-q(4:6),2)-...
                        data.m(2)*data.m(3) /norm(q(7:9)-q(4:6),2)-...
                        data.m(1)*data.m(3) /norm(q(7:9)-q(1:3),2));
    data.T     = @(p)   1/(2*data.m(1)) *norm(p(1:3),2)^2+...              
                        1/(2*data.m(2)) *norm(p(4:6),2)^2+...
                        1/(2*data.m(3)) *norm(p(7:9),2)^2;
    data.V     = @(q)-data.G*(data.m(1)*data.m(2) /norm(q(1:3)-q(4:6),2)...
                     -data.m(2)*data.m(3) /norm(q(7:9)-q(4:6),2)...
                     -data.m(1)*data.m(3) /norm(q(7:9)-q(1:3),2));
    data.gradT = @(p) [p(1)/data.m(1); p(2)/data.m(1); p(3)/data.m(1);...
                      p(4)/data.m(2); p(5)/data.m(2); p(6)/data.m(2);...
                      p(7)/data.m(3); p(8)/data.m(3); p(9)/data.m(3);];
    helpGradV  = @(qSelf,qOtherx,qOthery,mSelf,mOtherx,mOthery,index) ... 
                    data.G*mSelf*(mOtherx*(qSelf(index)-qOtherx(index))/...
                    norm(qSelf-qOtherx,2)^3 + mOthery*(qSelf(index)-qOthery(index))/...
                  norm(qSelf-qOthery,2)^3);
    data.gradV = @(q) [helpGradV(q(1:3), q(4:6), q(7:9), data.m(1), data.m(2), data.m(3), 1);
                       helpGradV(q(1:3), q(4:6), q(7:9), data.m(1), data.m(2), data.m(3), 2);
                       helpGradV(q(1:3), q(4:6), q(7:9), data.m(1), data.m(2), data.m(3), 3);
                       helpGradV(q(4:6), q(1:3), q(7:9), data.m(2), data.m(1), data.m(3), 1);
                       helpGradV(q(4:6), q(1:3), q(7:9), data.m(2), data.m(1), data.m(3), 2);
                       helpGradV(q(4:6), q(1:3), q(7:9), data.m(2), data.m(1), data.m(3), 3);
                       helpGradV(q(7:9), q(1:3), q(4:6), data.m(3), data.m(1), data.m(2), 1);
                       helpGradV(q(7:9), q(1:3), q(4:6), data.m(3), data.m(1), data.m(2), 2);
                       helpGradV(q(7:9), q(1:3), q(4:6), data.m(3), data.m(1), data.m(2), 3);
                      ];
        
else
    data.problemError_flag = 1;
    fprintf("The variable 'problem' only takes on the values 'osc', 'kep' and 'rtb'.\n")
end

end