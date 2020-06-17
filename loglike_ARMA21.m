% negative of the log likeliood for ARMA(2,2)
% input: x = [phi1 phi2 mu psi1 psi2]
% input: y = data; y0 = lags;
function ell = loglike_ARMA22(x,y0,y)
phi = zeros(3,1);
phi(1) = x(1); phi(2) = x(2);
phi(3) = x(3); % phi(3) = mu
psi1 = x(4); psi2 = x(5);
N = length(y); m = length(y0);
A = speye(N);
B = sparse(2:N, 1:N-1, ones(1,N-1), N,N);
C = sparse(3:N, 1:N-2, ones(1,N-2), N,N);
gam = A + B*psi1 + C*psi2;
gam2 = gam*gam';
X = [[y0(m);y(1:N-1)] [y0(m-1:end);y(1:N-2)] ones(N,1)];
ell = -(y-X*phi)'*(gam2\(y-X*phi));
ell = -ell;
