% negative of the log likeliood for ARMA(1,3)
% input: x = [phi1 mu psi1 psi2 psi3]
% input: y = data; y0 = lags;
function ell = loglike_ARMA13(x,y0,y)
phi = zeros(2,1);
phi(1) = x(1);
phi(3) = x(3); % phi(3) = mu
psi1 = x(4); psi2 = x(5);
N = length(y); m = length(y0);
A = speye(N);
B = sparse(2:N, 1:N-1, ones(1,N-1), N,N);
C = sparse(3:N, 1:N-2, ones(1,N-2), N,N);
D = sparse(4:N, 1:N-3, ones(1,N-3), N,N);
gam = A + B*psi1 + C*psi2 + D*psi3;
gam2 = gam*gam';
X = [[y0(m);y(1:N-1)] ones(N,1)];
ell = -(y-X*phi)'*(gam2\(y-X*phi));
ell = -ell;