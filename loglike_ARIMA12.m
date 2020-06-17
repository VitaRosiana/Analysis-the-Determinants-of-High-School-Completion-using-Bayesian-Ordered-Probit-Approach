%% negative of the log-lilkelihood for ARMA(1,2)
%% input: x = [mu, rho, psi(1), psi(2) sigma2];
%% input: y = data; y0 = lags;
function ell = loglike_ARIMA12(x,dely0,dely)
bet = zeros(2,1);
bet(1) = x(1); bet(2) = x(2); 
psi(1) = x(3);psi(2)=x(4); sigma2 = x(5);
T = length(dely); m = length(dely0);
A = speye(T); B = sparse(2:T,1:T-1, ones(1,T-1), T, T);
C = sparse(3:T,1:T-2, ones(1,T-2), T, T); 
Gam = A + B*psi(1) + C*psi(2); Gam2 = Gam*Gam';
X = [ones(T,1) [dely0(m);dely(1:T-1)] ];
%transformed ell removing constant term
ell = -log(det(sigma2*Gam2))-(1/sigma2)*(dely-X*bet)'*(Gam2\(dely-X*bet));
ell = -ell;