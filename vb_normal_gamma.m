%--------------------------------------------
% Demo for Variational Bayes: 
% use VB to approximate Normal-Gamma 
% written by Richard Xu
% July,2014
%--------------------------------------------

clc;
clf;
clear;

% ---------------------------------------------
% ground-truth parameters
% ---------------------------------------------

mu_0 = 0;
lambda_0 = 1;
a_0 = 1.5;
b_0 = 1;





T = 10;
n = 0;
k = a_0;
theta = 1/b_0;
tau = gamrnd(k, theta,[n 1]);


%Generate from a normal distribution with mu and standard deviation sigma:
%r = mu + sigma.*randn;

if n >0
   X = mu_0 + (tau * lambda_0).^(-0.5) .* randn(n,1);
else
   X = mu_0;
end   

mu_n = (lambda_0 * mu_0 + n * mean(X) )/ (lambda_0 + n);
lambda_n = lambda_0 + n;

a_n = a_0 + n/2;
b_n = b_0 + 1/2 * sum((X - mean(X)).^2) + (lambda_0*n * (mean(X) - mu_0)^2)/(2*(lambda_0 + n));


S_t = 100;
S_m = 100;

tau_pdf   = zeros(S_t,S_m);
mu_pdf    = zeros(S_t,S_m);
[mu_axis, tau_axis] = meshgrid(linspace(-5,5,S_m), linspace(0, 3 ,S_t));

%posterior for Normal-Gamma
%tau_n = gampdf(tau_axis, a_n,b_n);

tau_pdf = tau_axis.^(a_n -1).*exp(-b_n * tau_axis);

for i = 1:S_m
   mu_pdf (i,:) = normpdf(mu_axis(i,:),mu_n,(tau_pdf(i,1) * lambda_n).^(-0.5));
end

Z = mu_pdf .* tau_pdf ;


f1 = figure(1);
set(f1,'Name','variational Inferece to approximate Normal-Gamma distribution');

% -----------------------------------------------------
% plot ground truth distribution
% ------------------------------------------------------

subplot(1,2,1)
contour(tau_axis, mu_axis, Z);
title('ground truth distribution')
xlabel('\tau')
ylabel('\mu')

% -----------------------------
% Variational Approximation
% -----------------------------

init_mu_0      = -4;       
init_lambda_0  = 4;
init_a_0       = 1;       
init_b_0       = 3; 

% -------------------------

mu_prev     = init_mu_0;
lambda_prev = init_lambda_0;
a_prev      = init_a_0;
b_prev      = init_b_0;


for t =1:T
   
   if t ==1
       % the initial distribution for mu
       mu_current       = mu_prev;
       lambda_current   = lambda_prev;
       a_current = a_prev;
       b_current = b_prev;
       
   else
       
       %--------------------------------------
       % update q_tau (a_current, b_current)
       %--------------------------------------
       
       a_current  = a_0 + n/2;
       
       % E[mu^2] = var(mu) + (E[mu])^2
       E_mu_square = inv(lambda_0) + mu_prev^2;
       % E[mu]
       E_mu = mu_prev;
        
       %sum [(x_i - mu)^2]
       first = sum( X.^2 - 2 * X .* repmat(E_mu,size(X)) + repmat(E_mu_square, size(X)));
       
       %lambda_0 (mu - mu_0)^2
       second = lambda_0 *(E_mu_square - 2*mu_0*E_mu + mu_0^2);
       
       b_current = b_0 + (first + second)/2;
       
       %----------------------------------------------
       % update q_normal (mu_current, lambda_current)
       %----------------------------------------------
       
       %E[tau] = a/b;
       
       E_tau = a_current/b_current; 
       
       mu_current       = (lambda_0 * mu_0 + n * mean(X))/(lambda_0 + n);
       lambda_current   = (lambda_0 + n) * E_tau;
       
       % ------------------------------
       
       mu_prev      = mu_current;
       lambda_prev  = lambda_current;
       a_prev       = a_current;
       b_prev       = b_current;
       
   end
   

   tau_pdf = tau_axis.^(a_current -1).*exp(-b_current * tau_axis);

   mu_pdf = normpdf(mu_axis,mu_current, lambda_current^(-0.5)); 
   
   Z = mu_pdf .* tau_pdf ;
      
   subplot(1,2,2)
   contour(tau_axis, mu_axis, Z);
   title('Press to next iteration')
   xlabel('\tau')
   ylabel('\mu')
   
   waitforbuttonpress;
   
end

