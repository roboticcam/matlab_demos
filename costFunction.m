% ------------------------------------------------------
% Demonstraton of loss function using linear regression
% written by Richard Xu
% Jan 2017
% yida.xu@uts.edu.au
% ------------------------------------------------------

function costFunction

    
    % ------------------------------------------------------------
    % Generate linear regression dataset, adding noise to both x and y
    % ------------------------------------------------------------

    m = 0.5;
    b = 2;
    
    
    X = 0:0.5:10;
    X = X + (randn(length(X),1))';    
    Y = m* X + b + (randn(length(X),1))'* 0.5;
    
    
    
    figure(1);
    
    % ------------------------------------------------------------
    % Plot loss function across the whole range of values of (m,b)
    % ------------------------------------------------------------
    
    [val_m,val_b] = meshgrid(0:0.02:1, -3:0.02:3);
    val_m_data = repmat(val_m,[1 1 length(X) ]);
    val_b_data = repmat(val_b,[1 1 length(X) ]);
    
    val_X = repmat(X', [1 size(val_m,1) size(val_m,2)]);
    
    % ------------------------------------------------------------
    % put the dimensionality in the correct order
    % ------------------------------------------------------------
    
    val_X=permute(val_X,[2 3 1]); 
    
    
    val_Y = repmat(Y', [1 size(val_m,1) size(val_m,2)]);
    val_Y = permute(val_Y,[2 3 1]); 
    
    
    cost_func = val_Y - (val_m_data .* val_X + val_b_data);
    cost_func = sum(cost_func.^2,3) / length(X);
    mesh(val_m,val_b,cost_func);
    
    xlabel('m');
    ylabel('b');
    zlabel('average loss');

    
    figure(2);
    plot(X,Y,'+','color',[1 0 0]);
    
    xlabel('project duration');
    ylabel('cost');
    axis equal;
    axis( [ -5 15 0 10]);
    
    
    while(true)
        
        min_Y = m*min(X)+b;
        max_Y = m*max(X)+b;
        
        h1 = line([min(X) max(X)],[min_Y max_Y]);
        h1 = [h1; line([X' X']',[Y' (m*X+b)']','color',[0 1 0])];
        
        average_loss = sum((Y- (m*X+b)).^2)/length(X);
        xlabel(['project duration: Parameter: m = ' num2str(m,'%0.1f') ', b = ' num2str(b,'%0.1f') ...
                 ', average square loss = ' num2str(average_loss,'%0.1f') ]);
        
        m = rand;
        b = rand*5 - 2.5;
        
        waitforbuttonpress; 
        delete(h1);
    end
    
    
end