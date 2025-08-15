function amom = angMom3BP(approx)
% This function computes the angular momentum associated to each iteration
% of an integrator. It only takes a matrix containing the approximations as
% input and it produces a matrix of size 3 by N + 1, where N is the maximal
% amount of iterations to be computed. This computed matrix amom contains
% column-wise the angular momentum vector associated to an iteration.

N = size(approx,2);
amom = sparse(3, N);

for index = 1:N
    q1 = approx(1:3,index);
    q2 = approx(4:6,index);
    q3 = approx(7:9,index);
    p1 = approx(10:12,index);
    p2 = approx(13:15,index);
    p3 = approx(16:18,index);
    
    amom(:, index) = cross(q1,p1) + cross(q2,p2) + cross(q3,p3);
end

end

