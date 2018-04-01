% This is a little visualisation program, shows how one obtains efficient
% sampling using Hybrid-Monte Carlo:
%
% ( Duane, Simon; A.D. Kennedy, Brian J. Pendleton, and Duncan, Roweth (3).
% "Hybrid Monte Carlo". Physics Letters B 195 (2): 216–222.)
%
% The demo is sample Gaussian Mixture Model (GMM) using HMC. Note that you can sample
% GMM exactly, and there is absolutely no reason to use HMC. This is just
% for demo and visulisation only.
%
% Coded by Richard Xu @UTS (yida.xu@uts.edu.au)
%
% April 2012

function hybrid_gmm()

clear;
clc;


% GMM parameters
sigma(1) = 5.0;
sigma(2) = 3.0;
mu(1) = -6.0;
mu(2) = 6.0;
w(1) = 0.3;
w(2) = 0.7;


% start the sampling
N = 100;
q_samples = mvnrnd( mu(2), sigma(2), N);
Z = 1/N * sum( mvnpdf( q_samples, mu(2)-3, sigma(2)+2 ) ./ mvnpdf( q_samples, mu(2), sigma(2) ) );

PROPOSAL_SIGMA = 5.0;
i = 1;
accept(i) = 0.0;
rho = 3;

N = 500;

L = 5;

for t = 1:N
    
    x_0 = accept(i);
    
    u_star = mvnrnd( 0, 1);
   
    delta = each_step( x_0, mu, sigma, w);                              
    
    u_0 = u_star + rho/2 * delta; 
                                          
    for l=1:L
       
        if l == 1
           x(l) = x_0 + rho * u_0;
        else
           x(l) = x(l-1) + rho * u(l-1); 
        end
        
        if l == 1
            u(l) = u_0;
        else
            u(l) = u(l-1);
        end
        
        delta = each_step( x(l), mu, sigma, w);
        
        u(l) = u(l) + rho * delta;
                          
    end
    
    p_x_L   = w(1) * mvnpdf( x(L), mu(1), sigma(1)) + w(2) * mvnpdf( x(L), mu(2), sigma(2));
    p_x     = w(1) * mvnpdf( accept(i), mu(1), sigma(1)) + w(2) * mvnpdf( accept(i), mu(2), sigma(2));
    
    mh_ratio = min(1, ( p_x_L / p_x ) * exp ( -1/2 * ( u(L)^2 - u_star^2  )));

    v = rand;
    
    if v < mh_ratio
        
       accept(i+1) = x(L);
       i = i + 1;
    end

end

accept = accept(1,1:size(accept,2) * 0.9);

figure(1);

x_func = -15:0.1:12;
g1 = w(1) * 1/sigma(1) * exp( -1/2 * ( x_func - mu(1)) .^ 2 / sigma(1)^2 );
g2 = w(2) * 1/sigma(2) * exp( -1/2 * ( x_func - mu(2)) .^ 2 / sigma(2)^2 );
plot(x_func,g1+g2);

hold on;

first_term = (w(1)*(2*mu(1) - 2*x_func))./(2*sigma(1)^3*exp((mu(1) - x_func).^2/(2*sigma(1)^2))) + (w(2)*(2*mu(2) - 2*x_func))./( 2*sigma(2)^3*exp((mu(2) - x_func).^2/(2*sigma(2)^2)));
final_term = w(1)./(sigma(1)*exp((mu(1) - x_func).^2 /(2*sigma(1)^2))) + w(2)./(sigma(2)*exp((mu(2) - x_func).^2 /(2*sigma(2)^2)));

log_diff = first_term ./ final_term;

plot(x_func, log_diff, 'color',[1 0 0]);

hold on;
plot( accept(:), zeros([1 size(accept,2)])','.');
hold off;
fprintf('acceptance ratio = %f\n', size(accept,2)/((N+1)*0.9));


legend('f(x)', 'D log(f(x)) / D x', 'samples')

end

function step_value = each_step( xx, mu, sigma, w)
    step_value = ((w(1)*(2*mu(1) - 2*xx))/(2*sigma(1)^3*exp((mu(1) - xx)^2/(2*sigma(1)^2))) + (w(2)*(2*mu(2) - 2*xx))/(2*sigma(2)^3*exp((mu(2) - xx)^2/(2*sigma(2)^2))))/(w(1)/(sigma(1)*exp((mu(1) - xx)^2/(2*sigma(1)^2))) + w(2)/(sigma(2)*exp((mu(2) - xx)^2/(2*sigma(2)^2))));
end


