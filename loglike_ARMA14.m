%% ARMA(1,4) Process
% negative of the log likeliood for ARMA(1,4)
% input: x = [phi1 mu psi(1) psi(2) psi(3) psi(4) sig]
% input: y = data; y0 = lags;
function ell = loglike_ARMA14(x,y0,y)
phi = zeros(2,1);
phi(1) = x(1);
phi(2) = x(2);  %mu
psi (1) = x(3);
psi (2) = x(4);
psi (3) = x(5);
psi (4) = x(6);
sig2 = x(7);
N = length(y); m = length(y0);
A = speye(N);
B = sparse(2:N, 1:N-1, ones(1,N-1), N,N);
C = sparse(3:N, 1:N-2, ones(1,N-2), N,N);
D = sparse(4:N, 1:N-3, ones(1,N-3), N,N);
E = sparse(5:N, 1:N-4, ones(1,N-4), N,N);
gam = A + B*psi(1)+C*psi(2) + D*psi(3) + E*psi(4);
gam2 = gam*gam';
X = [[y0(m);y(1:N-1)] ones(N,1)];
ell = -(y-X*phi)'*(gam2\(y-X*phi));
ell = -ell;

