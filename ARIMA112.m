%% ARIMA (1,1,1)
clear
GDP = readtable ('Consumption.csv');
y = GDP (:,3) ; y = table2array (y);

dely = y(2:end)-y(1:end-1);
m = 4; 
T0 = 30-m; h=1;

N = length(y);
dely0 = dely(1:m); dely1 = dely(m+1:end); T = length(dely1);
ydelhat_ARIMA12 = zeros(T-h-T0+1,1);
yhat_ARIMA12 = zeros(T-h-T0+1,1);
ytph_ARIMA12 = y(m+T0+h: end);

f = @(x) loglike_ARIMA11(x,dely0,dely1(1:T0));
xhat_ARIMA12 = fminsearch(f,[0.5 ; 0.5 ; 0.5 ; 0.5; 0.5]);

for t = T0:T
    delyt = dely1(h:t);
    f = @(x) loglike_ARIMA12(x,dely0,delyt);
    xhat_ARIMA12 = fminsearch(f, xhat_ARIMA12);
    % make uhat
    A = speye(t); B = sparse(2:t,1:t-1,ones(1,t-1), t, t);
    C = sparse(3:t,1:t-2,ones(1,t-2), t, t);
    Gam = A + B*xhat_ARIMA12(3) + C*xhat_ARIMA12(4) ;
    X = [ones(t,1) [dely0(m);dely1(1:t-h)] ];
    uhat = Gam\(delyt-X*xhat_ARIMA11(1:2));
    ydelhat_ARIMA12(t-T0+h,:) = xhat_ARIMA12(1) + xhat_ARIMA12(2)*dely1(t)+ xhat_ARIMA12(3)*uhat(end)+ xhat_ARIMA12(4)*uhat(end-1);
    yhat_ARIMA12(t-T0+h,:) = ydelhat_ARIMA12(t-T0+h,:) + y(t+m);
end

MSFE_MA2 = mean((ytph_ARIMA12-yhatMA_ARIMA12).^2);
plot(y, 'k');
title('GDP Index Consumption');
xlabel('t');
ylabel('Consumption Index');
hold on
plot((T0:T), yhat_ARIMA12);
legend('GDP Consumption', '1-step ahead forecast ');
set(legend,'location','best')

err=ytph_ARIMA12-yhatMA_ARIMA12;