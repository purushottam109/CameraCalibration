function [R, T, K] = calibration2Dto3D(X, x)

[r c] = size(X);

A = zeros(2*c, 12);

for i = 1 : c
	A(2*i-1, 1) = X(1, i);
	A(2*i-1, 2) = X(2, i);
	A(2*i-1, 3) = X(3, i);
	A(2*i-1, 4) = 1;
	A(2*i-1, 9) = -x(1, i)*X(1, i);
	A(2*i-1, 10) = -x(1, i)*X(2, i);
	A(2*i-1, 11) = -x(1, i)*X(3, i);
	A(2*i-1, 12) = -x(1, i);
	
	A(2*i, 5) = X(1, i);
	A(2*i, 6) = X(2, i);
	A(2*i, 7) = X(3, i);
	A(2*i, 8) = 1;
	A(2*i, 9) = -x(2, i)*X(1, i);
	A(2*i, 10) = -x(2, i)*X(2, i);
	A(2*i, 11) = -x(2, i)*X(3, i);
	A(2*i, 12) = -x(2, i);
end

[u, s, v] = svd(A);

p = v(:, 12);

P34 = reshape(p, 4, 3)';

P33 = P34(1:3, 1:3); b = P34(:, 4);
a1 = P33(1, :)'; a2 = P33(2, :)'; a3 = P33(3, :)';

ro = 1.0/norm(a3);
u0 = ro*ro*(a1'*a2);
v0 = ro*ro*(a2'*a3);
t23 = cross(a2, a3);
beta = ro*ro*norm(t23);
s = (ro*ro*(a1'*a2) - u0*v0)/beta;
alpha = sqrt(ro*ro*(a1'*a1) - s*s - u0*u0);

r1 =t23/norm(t23);
r3 = ro*a3;
r2 = cross(r3, r1);

K = [alpha s u0; 0 beta v0; 0 0 1];
R = [r1'; r2'; r3'];
T = ro*K^(-1)*b;
