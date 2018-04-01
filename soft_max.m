% --------------------------------------------
% Demo for softmax (or multinomial) regression
% written by Richard Xu
% yida.xu@uts.edu.au, 
% April, 2016
% made it bare-minimum Octave-friendly
% --------------------------------------------

function soft_max()

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
data = zeros([N 3]);

for t = 1:M
   comp_size = length(find(Z == t));
   data(Z == t,1:2) = mvnrnd(mu(t,:), Sigma(:,:,t), comp_size);
%    data(Z == t,1:2) = mvnrnd2( Sigma(:,:,t), mu(t,:),comp_size,2);
   data(Z == t,3) = t;
end

Max_Iter = 12;
 
[x,y] = meshgrid(   min(data(:,1)):0.3:max(data(:,1)), ...
                    min(data(:,2)):0.3:max(data(:,2)));

                
data_mesh = [reshape(x,[],1),reshape(y,[],1)];
data_mesh = [data_mesh ones(size(data_mesh,1),1) ];


colors = [ 1 0 1; 0 0 1; 0 1 0];

  
fig1 = figure('KeyPressFcn',@keyPress, 'WindowButtonDownFcn', @mouseClicked_fig1);    
clf(fig1);

figure(fig1);

 
 for i =1:3
                   
    subset = data(data(:,3) == i,:);
    hh(i) = plot(subset(:,1), subset(:,2),'o', 'MarkerFaceColor', colors(i,:),'MarkerEdgeColor', colors(i,:));
 
    hold on;
 end 
 
legend(hh,{'class 1: [1 0 0]' 'class 2: [0 1 0]' 'class 3: [0 0 1]'});



 h_sel = [];

 function mouseClicked_fig1(objectHandle , eventData )
     
     if ~isempty(h_sel)
        delete(h_sel); 
     end
     
     coord = get(gca,'CurrentPoint');
     coord = coord(1,1:2);
       
     distances = sqrt(sum(bsxfun(@minus, data(:,1:2), coord).^2,2));
            
    [m sel] = min(distances);
     
    h_sel = plot(data(sel,1), data(sel,2),'o', 'MarkerEdgeColor', [1 0 0], 'LineWidth',5);
    
    str = '[1 0 0]';
    
    if data(sel,3) == 1
        str = '[1  0  0]';
    elseif data(sel,3) == 2
        str = '[0  1  0]';
    else
        str = '[0  0  1]';
    end
    
    str = ['cost(' str ' , '];
     
    a1 = exp([data(sel,1:2) 1]*thetas');
    probs1 = a1 ./ repmat(sum(a1,2),[1 3]);     
    
    str = [str [ '[' num2str(probs1(1),'%.3f') '  ' num2str(probs1(2),'%.3f') '  ' num2str(probs1(3),'%.3f')]];
    str = [str '] )'];
    h_sel = [h_sel; text(min(data(:,1)), max(data(:,2)),str)];
    
 end


thetas = [7 1 -1; 8 2 2; 6 3 3];


alpha = 0.001;


fig2 = figure();    
clf(fig2);


for n = 1:30
    
%         figure(2);
%         data_mesh = [reshape(x,[],1),reshape(y,[],1)];
%         data_mesh = [data_mesh ones(size(data_mesh,1),1) ];
            
        a = exp(data_mesh*thetas');
        probs = a ./ repmat(sum(a,2),[1 3]);     
        [val ind] = sort(probs','descend');
        
        for i =1:3
           
            subset = data_mesh(ind(1,:)== i,:);
            plot(subset(:,1), subset(:,2),'o', 'color', colors(i,:)*0.9);                        
            hold on;                 
            subset = data(data(:,3) == i,:);                
            hh(i) = plot(subset(:,1), subset(:,2),'o', 'MarkerFaceColor', colors(i,:),'MarkerEdgeColor', colors(i,:));
            hold on;           
        end
        
        legend(hh,{'class 1: [1 0 0]' 'class 2: [0 1 0]' 'class 3: [0 0 1]'});

        
%         for i = 1:3
%         thetas_mesh = reshape(thetas)';
%         a = exp([x,y]*thetas');
%         
%         end
% 
%         h =1;
        
%     data_t = [data(:,1:2) ones(size(data,1),1)];
%     a = exp(data_t*thetas');    
%     for i =1:3
%        prob = a(:,i)./sum(a,2);     
%        [val ind] = sort(prob,'descend');    
%        ind = ind(1:length(find(data(:,3)==i)));
%        subset = data(ind,:);
%        gg(i) = plot(subset(:,1), subset(:,2),'o', 'color', colors(i,:));
%        hold on; 
%     end
    
%     legend(gg,{'class 1: [1 0 0]' 'class 2: [0 1 0]' 'class 3: [0 0 1]'});
    
%     Cp = get(gca,'CurrentPoint');
    
    waitforbuttonpress;
    
    
    for i =1:3

        r_data = data(data(:,3) == i,:);
        r_data(:,3) = 1;
        a = exp(r_data*thetas');
        prob = a(:,i)./sum(a,2);   
        thetas(i,:) = thetas(i,:) - alpha * ( sum(r_data .* repmat((prob -1),[1,3])) + thetas(i,:));    
        
    end
    
    
end
end

