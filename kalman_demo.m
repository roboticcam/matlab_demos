% --------------------------------------
% Kalman Filter demos
% written by Richard Xu
% yida.xu@uts.edu.au, July, 2014
% --------------------------------------

clc;
clear;


% --------------------------------------------------------------------
% Synthetically generate data from nominated Kalman Filter parameters
% --------------------------------------------------------------------


% for p(x_1)
init_mu  = [ 0 0]';
init_var = [ 5 0; 0 5];      


% transition probability
A = [1 0; 0 1];

B = [0.2 0.1]';

Q = [0.1 0.05; 0.05 0.1];
%Q = [1E-5 0; 0 1E-5];

% fixed parameters
H = [1 0; 0 1];
R = [0.5 -0.05; -0.05 0.5];
%R = [1E-5 0 ; 0 1E-5];

T = 16;

xs(1,:) = [0 0];

for t=2:T
   xs(t,:) = (A * xs (t-1,:)')' + B' + mvnrnd([0 0], Q );
end

% data is a sequence of y
data = xs + mvnrnd([0 0],R,T);

data = data';


%update = Kalman_posterior (data, A, B, Q, H, R, T, init_mu, init_var);
   

update_mu       = zeros(2,T);
update_Sigma    = zeros(2,2,T);

predict_mu       = zeros(2,T);
predict_Sigma    = zeros(2,2,T);



for t=1:T
    
    
    % --------------------------------
    % Kalman Filter update
    % --------------------------------
    
    if t == 1
        predict_mu(:,1)         = init_mu;
        predict_Sigma(:,:,1)    = init_var;
    else
        predict_mu(:,t)         = A * update_mu(:,t-1) + B;
        predict_Sigma(:,:,t)    = A * update_Sigma(:,:,t-1) * A' + Q;
    end

    S = H * predict_Sigma(:,:,t) * H' + R;
    K = predict_Sigma(:,:,t) * H' * inv(S);
   
    update_mu(:,t)      = predict_mu(:,t) + K * (data(:,t) - H * predict_mu(:,t));
    update_Sigma(:,:,t) = (eye(size(R,2)) - K * H ) * predict_Sigma(:,:,t);
    
    %display(update_mu(:,t));
    display(update_Sigma(:,:,t));
    
   
    % --------------------------------
    % The plot
    % --------------------------------
   
    
    gmm_pdf = gmdistribution(update_mu(:,t)',update_Sigma(:,:,t),1); 
    ezcontour(@(u,v)pdf(gmm_pdf,[u v]));
    hold on;
    
    plot(data(1,1:t),data(2,1:t),'--','color',[0 0 1]);
    axis([min(data(1,:))-1 max(data(1,:))+1 min(data(2,:))-1 max(data(2,:))+1]); 
    hold on;
    
    plot(update_mu(1,1:t),update_mu(2,1:t),'--','color',[1 0 0]);
    text(update_mu(1,t),update_mu(2,t), num2str(t), 'color', [ 1 0 0]);
    
    axis([min(data(1,:))-1 max(data(1,:))+1 min(data(2,:))-1 max(data(2,:))+1]); 
    hold on;
    
    plot(data(1,t:T),data(2,t:T),'--','color',[0.7 0.7 0.7]);
    axis([min(data(1,:))-1 max(data(1,:))+1 min(data(2,:))-1 max(data(2,:))+1]); 
    hold on;
    
    plot(data(1,t),data(2,t),'*');
    
    for s=1:T
        text(data(1,s),data(2,s), num2str(s), 'color',[ 0 0 1]);
        hold on;
    end
    
    
    
    hold off;
    waitforbuttonpress;
end



