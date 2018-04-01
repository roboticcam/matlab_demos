% ----------------------------------
% A simple illustration of LDA
% (1) generate W according to the LDA's generative model
% (2) Perform Gibbs Sampling
% (3) part of teaching for USTC computation statistics
% 
% By Richard Xu
% yida.xu@uts.edu.au
%
% 2017/04/10
%
% ----------------------------------

function lda_gibbs_example

clc;
clear;

%check_dirrand

% K is number of topics
K = 5;

% alpha = [2 ... 2] with K elements
alphas = ones(1,K)*2;

% V is the total number of words in a corpus.
V = 20;

% eta = [ 2 ...  2] with V elements
etas = ones(1,V)*2;

% D is the number of documents
D = 6;

% N is number of words in a document
N = 10;

% ------------------------------------------------
% Generate W, the words data
% ------------------------------------------------

[betas, thetas, Z, W] = generate_LDA(K,V,D,N, alphas, etas);


% ------------------------------------------------
% Now we have W, start Gibbs Sampling
% ------------------------------------------------

[betas, thetas, Z, temp] = generate_LDA(K,V,D,N, alphas, etas);

    

    for g=1:100
        
        % sample theta
      
        for d = 1:D
            
            I = zeros(1,K);
            for i = 1:K
               I(i) = length(find(Z(d,:) == i));
            end

            thetas(d,:) = dirrnd(alphas + I) ;  
          
        end
        
        
        % sample beta       
        I = zeros(1,V);
        for i = 1:V
           I(i) = length(find(W == i));
        end
        
        for i = 1:K
            betas(i,:) = dirrnd(etas + I);
        end
        % sample Z
        
        for d =1:D
            
            for n =1:N
                
                probs = zeros(1,K);
                for k =1:K
                   
                    %Pr( Z_{d,n} = k)
                   
                    probs(k) = thetas(d,k) * betas(k,W(d,n));                   
                end
                
                probs = probs / sum(probs);
                
                Z(d,n) = find(mnrnd(1,probs));
                
            end
            
        end
        
        fprintf('at %d iterations\n',g);
        display(thetas);
        display(betas);
        display(Z);
        

    end

end

% ------------------------------------------------------------------
% Generate LDA according to the generative model
% ------------------------------------------------------------------
function [betas, thetas, Z, W] = generate_LDA(K,V,D,N,alphas, etas)
    
    % beta is the word distribution per topic
    betas = zeros(K,V);
    for k = 1:K
       betas(k,:) = dirrnd(etas);
    end

    Z = zeros(D,N);
    W = zeros(D,N);
    thetas = zeros(D,K);

    for d = 1:D

       % for each document d, theta is the topic distribution for that document
       thetas(d,:) = dirrnd(alphas);

       temp = mnrnd(1,thetas(d,:),N);

       for n =1:N
           Z(d,n) = find(temp(n,:));
           W(d,n) = find(mnrnd(1,betas(Z(d,n),:)),1);
       end

    end
end

% ------------------------------------------------------------------
% check to see if dirichlet distribution is generated properly
% ------------------------------------------------------------------
function check_dirrand
    
    alphas = [1 2 3];

    for i =1:10000
        g(i,:) = dirrnd(alphas);
    end

    % mean of Dirichlet distribution should be alpha/sum(alpha)
    display(mean(g));

    % This is the emperical variance.
    display(var(g));

    % This is the theoretical variance.
    emp_var = alphas .* ( repmat(sum(alphas),[1 3]) - alphas)/sum(alphas)^2/(sum(alphas) + 1); 

    display(emp_var);

end

% ------------------------------------------------------------------
% Since MATLAB does not provide any method for sampling dirichlet
% distribution, we need to sample it using Gamma distribution.
% ------------------------------------------------------------------


function g = dirrnd(alphas)

    g = gamrnd(alphas,1);
    g = g/sum(g);

end
