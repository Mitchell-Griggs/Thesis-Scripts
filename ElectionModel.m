%% ElectionModel simulates the election model from Section 8.3, for selected parameters.

%% Parameters:a
hR = 2^-5;  % Initial step size for the Wiener processes.
d1 = 0.1; d2 = 0.2; % Time delays.
X0 = [0.5;0];   % Initial popularity and initial expense.
strat = @(t,x2) x2; % Campaigning strategy.
opp = @(t,x1d2) 1-heaviside(x1d2-0.2);  % Opposition strategy. 
fin = @(t,x1d1) 1-heaviside(x1d1-0.5);  % Financial costs.
sigma1 = 0.2; sigma2 = 0.9; % Volatility terms.

%% Code:

T = 1; % Terminal time.

% Create the augmented mesh:
mesh0 = hR:hR:T;
NrValuesPotentially = 1:100; % The upper bound Nr on this needs to have Nr*r >= T.
N1 = find(NrValuesPotentially*d1>=T,1);
N2 = find(NrValuesPotentially*d2>=T,1);
TO = [mesh0,(1:N1)*d1,(1:N2)*d2];
TO = sort(TO);    TO = TO(0<TO&TO<=T);    TO = unique(TO);
TC = TO;
for i = 0:N1
    for j = 0:N2
        TC = [TC,TC - i*d1 - j*d2];
        TC = sort(TC);    TC = TC(0<TC&TC<=T);    TC = unique(round(TC,10));
    end
end
TC = sort(TC);    TC = TC(0<TC&TC<=T);    TC = unique(TC);
t = TC; % Time.
N = length(t);

rng('default')

% Generate Wiener paths:
dW1mesh0initial = sqrt(hR)*randn;
W1mesh0 = zeros(1,length(mesh0));
W1mesh0(1,1) = dW1mesh0initial;
W1mesh0(1,2:end) = dW1mesh0initial + cumsum(sqrt(hR)*randn(1,length(mesh0)-1));
dW1mesh0 = zeros(1,length(mesh0)-1);
for n = 1:length(mesh0)-1
    dW1mesh0(n) = W1mesh0(n+1)-W1mesh0(n);
end
% To test this:
% figure(7); plot(mesh0(1:end),cumsum([dW1mesh0initial,dW1mesh0]),':r',[0,mesh0],[0,W1mesh0],'--b'); legend('sum of $dW_1(t)$','$W_1(t)$','Interpreter','latex')
W1 = zeros(1,N);
for n = 1:N
    if (t(n)==mesh0)==zeros(1,length(mesh0))
        if n == 1
            indexnext = find(round(t,5)==round(mesh0(1,1),5),1);
            mu = 0 + (t(1)-0)/(t(indexnext)-0)*(W1mesh0(1,1)-0);
            sigma = sqrt((t(indexnext) - t(1)) * (t(1) - 0) / (t(indexnext) - 0));
            W1(1,n) = normrnd(mu,sigma);
        elseif n < find(round(t,5)==round(mesh0(1),5),1)
            mu = W1(1,n-1) + (t(n)-t(n-1))/(t(indexnext)-t(n-1))*(W1mesh0(1,1)-W1(1,n-1));
            sigma = sqrt((t(indexnext) - t(n)) * (t(n) - t(n-1)) / (t(indexnext) - t(n)));
            W1(1,n) = normrnd(mu,sigma);
        elseif n == find(round(t,5)==round(mesh0(1),5),1) 
            W1(1,n) = W1mesh0(1,1);
        elseif n > find(round(t,5)==round(mesh0(1),5),1) 
        if isempty(find(round(mesh0,5)==round(t(n),5),1))
            meshindexnext = find(mesh0>t(n),1);
            indexnext = find(round(t,5)==round(mesh0(1,meshindexnext),5),1);
            mu = W1(1,n-1) + (t(n)-t(n-1))/(t(indexnext)-t(n-1))*(W1mesh0(1,meshindexnext)-W1(1,n-1));
            sigma = sqrt((t(indexnext) - t(n)) * (t(n) - t(n-1)) / (t(indexnext) - t(n)));
            W1(1,n) = normrnd(mu,sigma);
        else
            W1(1,n) = W1mesh0(1,meshindexnext);
        end
        end
    else
        index = find(round(mesh0,5)==round(t(n),5),1);
        W1(1,n) = W1mesh0(1,index);
    end
end
dW1 = zeros(1,length(t));
for n = 1:N-1
    dW1(1,n) = W1(1,n+1) - W1(1,n);
end

% W2:
dW2mesh0initial = sqrt(hR)*randn;
W2mesh0 = zeros(1,length(mesh0));
W2mesh0(1,1) = dW2mesh0initial;
W2mesh0(1,2:end) = dW2mesh0initial + cumsum(sqrt(hR)*randn(1,length(mesh0)-1));
dW2mesh0 = zeros(1,length(mesh0)-1);
for n = 1:length(mesh0)-1
    dW2mesh0(n) = W2mesh0(n+1)-W2mesh0(n);
end
% To test this:
% figure(8); plot(mesh0(1:end),cumsum([dW2mesh0initial,dW2mesh0]),':r',[0,mesh0],[0,W2mesh0],'--b'); legend('sum of $dW_2(t)$','$W_2(t)$','Interpreter','latex')
W2 = zeros(1,N);
for n = 1:N
    if (t(n)==mesh0)==zeros(1,length(mesh0))
        if n == 1
            indexnext = find(round(t,5)==round(mesh0(1,1),5),1);
            mu = 0 + (t(1)-0)/(t(indexnext)-0)*(W2mesh0(1,1)-0);
            sigma = sqrt((t(indexnext) - t(1)) * (t(1) - 0) / (t(indexnext) - 0));
            W2(1,n) = normrnd(mu,sigma);
        elseif n < find(round(t,5)==round(mesh0(1),5),1)
            mu = W2(1,n-1) + (t(n)-t(n-1))/(t(indexnext)-t(n-1))*(W2mesh0(1,1)-W2(1,n-1));
            sigma = sqrt((t(indexnext) - t(n)) * (t(n) - t(n-1)) / (t(indexnext) - t(n)));
            W2(1,n) = normrnd(mu,sigma);
        elseif n == find(round(t,5)==round(mesh0(1),5),1) 
            W2(1,n) = W2mesh0(1,1);
        elseif n > find(round(t,5)==round(mesh0(1),5),1) 
        if isempty(find(round(mesh0,5)==round(t(n),5),1))
            meshindexnext = find(mesh0>t(n),1);
            indexnext = find(round(t,5)==round(mesh0(1,meshindexnext),5),1);
            mu = W2(1,n-1) + (t(n)-t(n-1))/(t(indexnext)-t(n-1))*(W2mesh0(1,meshindexnext)-W2(1,n-1));
            sigma = sqrt((t(indexnext) - t(n)) * (t(n) - t(n-1)) / (t(indexnext) - t(n)));
            W2(1,n) = normrnd(mu,sigma);
        else
            W2(1,n) = W2mesh0(1,meshindexnext);
        end
        end
    else
        index = find(round(mesh0,5)==round(t(n),5),1);
        W2(1,n) = W2mesh0(1,index);
    end
end
dW2 = zeros(1,length(t));
for n = 1:N-1
    dW2(1,n) = W2(1,n+1) - W2(1,n);
end

% Construct X(t) = (X1(t),X2(t)):

X = zeros(2,length(t));
for n = 0:N-1
    if n == 0
        tn = 0;
        h = t(1);
        Xn = X0;
        Xd1 = X0;
        Xd2 = X0;
        f = [strat(tn,Xn(2))*(1-Xn(1)) + opp(tn,Xd2(1))*Xn(1)*(Xn(1)-1) ; fin(tn,Xd1(1))];
        g1 = [sigma1*(1-Xn(1)) ; 0];
        g2 = [sigma2*Xn(1)*(1-Xn(1)) ; 0];
        X(:,n+1) = Xn + f*h + g1*dW1(1,n+1) + g2*dW2(1,n+1);
    else
        tn = t(n);
        h = t(n+1)-t(n);
        Xn = X(:,n);
        if tn <= d1
            Xd1 = X0;
        else
            td1index = find(round(t,5)==round(tn-d1,5),1); Xd1 = X(:,td1index);
        end
        if tn <= d2
            Xd2 = X0;
        else
            td2index = find(round(t,5)==round(tn-d2,5),1); Xd2 = X(:,td2index);
        end
        f = [strat(tn,Xn(2))*(1-Xn(1)) + opp(tn,Xd2(1))*Xn(1)*(Xn(1)-1) ; fin(tn,Xd1(1))];
        g1 = [sigma1*(1-Xn(1)) ; 0];
        g2 = [sigma2*Xn(1)*(1-Xn(1)) ; 0];
        X(:,n+1) = Xn + f*h + g1*dW1(1,n+1) + g2*dW2(1,n+1);
    end
end


figure(1)
hold on
plot([0,t],[X0(1),X(1,:)],'k','LineWidth',2)
plot([0,t],[X0(2),X(2,:)],'r','LineWidth',1)
plot([-0.2,1],[0.5,0.5],'k:')
legend('Polling ($X_1(t)$)','Expenses ($X_2(t)$)','Location','Northeast','Interpreter','latex')
set(gca,'xtick',[0,T],'FontSize',16)
xticklabels({'0','Election','Interpreter','latex'})
xlabel('$t$','Interpreter','latex','FontSize',16)
set(gca,'ytick',[0,0.25,0.5,0.75,1],'FontSize',16)
ylabel('Polling','Interpreter','latex','FontSize',16)
title('SDDE Model of Popularity','Interpreter','latex','FontSize',16)
axis([0,1.05,-0.05,1.05])
set(gcf,'position',[200,100,500,400])
