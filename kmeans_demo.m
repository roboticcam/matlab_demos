function kmeans_demo()

% ---------------------------------------------------
% K-means demos: This is a Octave friendly function
% written by Richard Xu
% yida.xu@uts.edu.au
% ---------------------------------------------------

clear all;
clc;

% --------------------------------------
% acts like a multinomial distribution
% --------------------------------------
N       = 1000;
weights = [0.33 0.33 0.34];
M       = length(weights);

ind = rand([N 1]);
Z = zeros([N 1]);
acc_weights = cumsum(weights);

for t = 1:M
   Z(ind > acc_weights(t)) = t;
end
Z = Z + 1;




% ------------------------------------------
% Generate data from Gaussian Mixture model
% ------------------------------------------

mu = [ -3 4; 0 0; 3 4];
Sigma(:,:,1) = [ 2 0; 0 1];
Sigma(:,:,2) = [ 2 0; 0 1];
Sigma(:,:,3) = [ 2 0; 0 3];
data = zeros([N 2]);

for t = 1:M
    comp_size = length(find(Z == t));
    data(Z == t,:) = mvnrnd(mu(t,:), Sigma(:,:,t), comp_size);
end



% ------------------------------------------
% K-means starts from here
% ------------------------------------------

Max_Iter = 12;

% ------------------------------------------
% draw before the first iteration
% ------------------------------------------

plot(data(:,1), data(:,2),'.', 'color',[1 0 1]);
hold on;
ctrs = [1 3; 1 1; 0 -2];
plot(ctrs(:,1),ctrs(:,2),'kx','MarkerSize',12,'LineWidth',2)
plot(ctrs(:,1),ctrs(:,2),'ko','MarkerSize',12,'LineWidth',2)
hold off;
waitforbuttonpress;


colors = {'r.', 'g.', 'b.'};

for i = 1:Max_Iter
    
    % ----------------------------------------------------
    %    MATLAB has a kmeans(),but we implement ourselves
    % ----------------------------------------------------

    %    [idx,ctrs] = kmeans(data,3,'Distance','city','MaxIter',1,'start',ctrs);
    
    for j =1:size(ctrs,1)
        mean_val = repmat(ctrs(j,:),[size(data, 1) 1]);
        distance(j,:) = diag((data - mean_val)*(data - mean_val)');
    end
    
    [temp idx] = min(distance);
    
    for j =1:size(ctrs,1)
        ctrs(j,:) = mean(data(idx==j,:));
        plot(data(idx==j,1),data(idx==j,2),colors{j},'MarkerSize',12);
        hold on;
    end
    
    plot(ctrs(:,1),ctrs(:,2),'kx','MarkerSize',12,'LineWidth',2)
    plot(ctrs(:,1),ctrs(:,2),'ko','MarkerSize',12,'LineWidth',2)
    
    legend('type 1','type 2','type 3', 'Centroids', 'Location','NW')
    
    waitforbuttonpress;
    hold off;

end
end