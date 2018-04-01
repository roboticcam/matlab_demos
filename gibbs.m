% ---------------------------
% Simulation of a Gibbs sampling
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

mu = [2 3];
Sigma = [ 3  2; 2 5];


figure(1);

N = 2000;

init_sample = [ 1 1];

directX = mvnrnd(mu,Sigma,N);

minX = min(directX(:,1));
minY = min(directX(:,2));
maxX = max(directX(:,1));
maxY = max(directX(:,2));

figure(1);

subplot(2,2,1);
plot(directX(:,1),directX(:,2),'.','color',[1 0 0]);
axis([minX maxX minY maxY]);

subplot(2,2,2);
autocorr(mvnpdf(directX,mu,Sigma));



HIDE_GIBBS_STEP = false;

if HIDE_GIBBS_STEP == true

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
    axis([minX maxX minY maxY]);
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
        axis([minX maxX minY maxY]);
        hold on;

        waitforbuttonpress;
        %waitfor(100);
        mu_2 = mu(2) + Sigma(2,1)/Sigma(1,1)*(X(i,1) - mu(1));
        Sigma_2 = Sigma(2,2) - Sigma(2,1)/Sigma(1,1)*Sigma(1,2);

        X(i,2) = randn * sqrt(Sigma_2) + mu_2;      
        plot( [X(i,1) X(i,1)], [X(i-1,2) X(i,2)],'LineWidth',1);
        axis([minX maxX minY maxY]);

        hold on;

        waitforbuttonpress;
        %waitfor(100);
    end


    plot(X(:,1),X(:,2),'.');

end

subplot(2,2,4);
autocorr(mvnpdf(X,mu,Sigma));

