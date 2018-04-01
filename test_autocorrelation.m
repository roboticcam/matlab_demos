% ---------------------------
% Simulation of a Autocorrelation
% written by Richard Xu
% March 2013
% yida.xu@uts.edu.au
% --------------------------
% clear;
% m = 4; k=3; n = 5; A = randn(m,k); B = randn(k,n); % values
% a1 = A * B; % method 1
% display(a1);
% 
% a2 = zeros(m,n);
% for i=1:k;  % method 2
%     a2 = a2+ A(:,i) * B(i,:);  
% end;
% display(a2);
% 
% return;


clc;
clf;
clear;


% mu = [2 3];
% Sigma = [ 3  2; 2 5];

D = 5;

mu_full   = randn(1, D) * 2;
Sigma_full  = wishrnd(eye(D),5)/2;

display(mu_full);
display(Sigma_full);


% Sigma  = eye(D);
N = 2000;

X_full = mvnrnd(mu_full,Sigma_full,N);

figure(1);

for i = 2:D
    
    X       = X_full(:,1:i);
    mu      = mu_full(1:i);
    Sigma   = Sigma_full(1:i,1:i);
    
  
    % --------------------------------
    % sample it directly
    % --------------------------------
    
    subplot(D-1,2,(i-1)*2);
    autocorr(mvnpdf(X,mu,Sigma))
    title(['direct sampling of ' num2str(i) 'D Gaussian']);
    
    
    % --------------------------------
    % using Gibbs sampling
    % --------------------------------
    
    
    X = zeros(N,i);
    X(1,:) = mu;

    for j = 2:N
        
        X(j,:) = X(j-1,:);
        
        for k = 1:i
            
            
            Sigma_12 = Sigma(k,:);
            Sigma_12(k) = [];

            Sigma_21 = Sigma(:,k);
            Sigma_21(k) = [];
            
            Sigma_22 = Sigma;
            Sigma_22(k,:) = [];
            Sigma_22(:,k) = [];
            
            mu_2    = mu;
            mu_2(k) = [];
            
            
            X2      = X(j,:);
            X2(k)   = [];
            
            mu_k    = mu(k) + Sigma_12*inv(Sigma_22)*( X2 - mu_2)';
            Sigma_k = Sigma(k,k) - Sigma_12*inv(Sigma_22)*Sigma_21;

            X(j,k) = randn * sqrt(Sigma_k) + mu_k;
 
        end
        
    end 
    
    subplot(D-1,2,(i-1)*2-1);
    autocorr(mvnpdf(X,mu,Sigma))
    title(['gibbs sampling of ' num2str(i) 'D Gaussian']);

    
    display(mean(X));
    display(cov(X));
    
    
end


return;


figure(1);

N = 2000;

init_sample = [ 1 1];

directX = mvnrnd(mu,Sigma,N);

figure(1);
subplot(2,2,1);
plot(directX(:,1),directX(:,2),'.','color',[1 0 0]);

subplot(2,2,2);
autocorr(mvnpdf(directX,mu,Sigma));

SHOW_STEP = true;

if SHOW_STEP == true

    X = zeros(N,2);
    X(1,:) = init_sample;

    for i = 2:N

        mu_1    = mu(1) + Sigma(1,2)/Sigma(2,2)*( X(i-1,2) - mu(2));
        Sigma_1 = Sigma(1,1) - Sigma(1,2)/Sigma(2,2)*Sigma(2,1);

        X(i,1) = randn * sqrt(Sigma_1) + mu_1;

        mu_2 = mu(2) + Sigma(2,1)/Sigma(1,1)*(X(i,1) - mu(1));
        Sigma_2 = Sigma(2,2) - Sigma(2,1)/Sigma(1,1)*Sigma(1,2);

        X(i,2) = randn * sqrt(Sigma_2) + mu_2;
        

    end

    subplot(2,2,3);
    plot(X(:,1),X(:,2),'.');

    %---------------------------------------

else
    
    
    subplot(2,2,3);
    
    X = zeros(N,2);
    X(1,:) = [1 1];

    for i = 2:N

        mu_1    = mu(1) + Sigma(1,2)/Sigma(2,2)*( X(i-1,2) - mu(2));
        Sigma_1 = Sigma(1,1) - Sigma(1,2)/Sigma(2,2)*Sigma(2,1);

        X(i,1) = randn * sqrt(Sigma_1) + mu_1;
        
        plot( [X(i-1,1) X(i,1)], [X(i-1,2) X(i-1,2)],'LineWidth',1);
        axis([-6 8 -6 12]);
        hold on;

        waitforbuttonpress;
        %waitfor(100);
        mu_2 = mu(2) + Sigma(2,1)/Sigma(1,1)*(X(i,1) - mu(1));
        Sigma_2 = Sigma(2,2) - Sigma(2,1)/Sigma(1,1)*Sigma(1,2);

        X(i,2) = randn * sqrt(Sigma_2) + mu_2;      
        plot( [X(i,1) X(i,1)], [X(i-1,2) X(i,2)],'LineWidth',1);
        axis([-6 8 -6 12]);
        hold on;

        waitforbuttonpress;
        %waitfor(100);
    end


    plot(X(:,1),X(:,2),'.');

end

subplot(2,2,4);
autocorr(mvnpdf(X,mu,Sigma));

