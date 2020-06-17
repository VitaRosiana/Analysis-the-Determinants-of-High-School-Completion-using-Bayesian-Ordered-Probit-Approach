%% ARIMA (1,1,1)
clear
GDP = readtable ('Consumption.csv');
y = GDP (:,3) ; y = table2array (y);

dely = y(2:end)-y(1:end-1);
m = 4; 
T0 = 30-m; h=1;

N = length(y);
dely0 = dely(1:m); dely1 = dely(m+1:end); T = length(dely1);
ydelhat_ARIMA11 = zeros(T-h-T0+1,1);
yhat_ARIMA11 = zeros(T-h-T0+1,1);
ytph_ARIMA11 = y(m+T0+h: end);

f = @(x) loglike_ARIMA11(x,dely0,dely1(1:T0));
xhat_ARIMA11 = fminsearch(f,[0 ; 0 ; 0 ; 0]);

for t = T0:T
    delyt = dely1(h:t);
    f = @(x) loglike_ARIMA11(x,dely0,delyt);
    xhat_ARIMA11 = fminsearch(f, xhat_ARIMA11);
    % make uhat
    A = speye(t); B = sparse(2:t,1:t-1,ones(1,t-1), t, t);
    Gam = A + B*xhat_ARIMA11(3);
    X = [ones(t,1) [dely0(m);dely1(1:t-h)] ];
    uhat = Gam\(delyt-X*xhat_ARIMA11(1:2));
    ydelhat_ARIMA11(t-T0+h,:) = xhat_ARIMA11(1) + xhat_ARIMA11(2)*dely1(t)+ xhat_ARIMA11(3)*uhat(end);
    yhat_ARIMA11(t-T0+h,:) = ydelhat_ARIMA11(t-T0+h,:) + y(t+m);
end

err=(ytph_ARIMA11 - yhat_ARIMA11);
MSFE_ARIMA11 = mean(err.^2);
n=length(ytph_ARIMA11);
plot(1:n,ytph_ARIMA11,1:n,yhat_ARIMA11)
% plot(1:n,err)