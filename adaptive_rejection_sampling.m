%--------------------------------------------
% Demo to illustrate Adaptive Rejection Sampling
%
% written by Richard Xu
% Yida.Xu@uts.edu.au
% July,2012
%--------------------------------------------


function s = adaptive_rejection_sampling()

clc;
clear;

% --------------------------------------------------
% this is to sample Gaussian distribution
% note that although gaussian is log-concave, ARS is an over-kill for 
% sampling it, this is for pure demo purpose
% --------------------------------------------------

mu = 1;
sigma = 3;
log_func    = @(x) (-1/2 * (x - mu).^2/sigma);
log_dev     = @(x) (- (x-mu)/sigma);
%func        = @(x) exp(-1/2 * (x - mu).^2/sigma);
func        = @(x) exp( (-1/2 * (x - mu).^2/sigma) );
init_points = [4 ; -1.5; 2; -3];

bounds(1) = -inf;
bounds(2) = inf; 

plot_bounds(1) = -10;
plot_bounds(2) = 10;


% ------------------------------------------------------------------
% please adjust this parameter to see how the envelop was constructed
% ------------------------------------------------------------------

MAX_ACCEPTED = 10;
IS_PLOT = 1;


% -------------------------------------------
% this is to sample posteior of DP's concentration parameter
% -------------------------------------------

% w = 4;
% s = ceil(rand([10 1]) * 10);
% w_s_j_term_log =  sum(log(s) + log(w) - s * w);
% k = length(s);
% 
% w_s_j_term_log = - 22;
% %k = 20;
% 
% offset = 0;
% 
% log_func     = @(b)(    - k * gammaln(b/2) ...
%                          - 1 ./(2*b)...
%                          + (k * b - 3)/2 .* log(b/2)... 
%                          + b/2 * w_s_j_term_log + offset );
% 
% log_dev      = @(b)(    - k/2 * psi(b/2) ...
%                          + 1/2 * b.^(-2) ... % 1/2 * b^-1 => -1/2 * b^-2
%                          + k/2 * log(b/2)...
%                          + (k * b - 3)./(2 * b)... log(b/2) => 2/b * 1/2 = 1/b
%                          + w_s_j_term_log /2 );
%                      
%     func     = @(b)(exp( - k * gammaln(b/2) ...
%                          - 1 ./(2*b)...
%                          + (k * b - 3)/2 .* log(b/2)... 
%                          + b/2 * w_s_j_term_log + offset));
%  
% log_func_neg  = @(b)( -(   - k * gammaln(b/2) ...
%                          - 1 ./(2*b)...
%                          + (k * b - 3)/2 .* log(b/2)... 
%                          + b/2 * w_s_j_term_log + offset));
%                      
% max_pt = fminbnd(log_func_neg,0,100);
%                      
% bounds(1) = 0;
% bounds(2) = inf;
% 
% plot_bounds(1) = 0;
% plot_bounds(2) = 4* max_pt;
% 
% init_points = [max_pt/2;  max_pt + max_pt/2];
% MAX_ACCEPTED = 100;
% IS_PLOT = 1;                    

% -------------------------------------------                    

[samples num_sampled] = ARS(func, log_func, log_dev, init_points, MAX_ACCEPTED, IS_PLOT, bounds, plot_bounds);
fprintf('acceptance ratio = %f',MAX_ACCEPTED / num_sampled);


function [samples num_sampled] = ARS(func, log_func, log_dev, init_points, MAX_ACCEPTED, IS_PLOT, bounds, plot_bounds)


    % (x,f) are the points on the log_func curve
    % z is the intersection points
    % m, b are gradient and y-intecept of tangenet
    x = init_points;

    [m b x z f ] = get_piecewise_exp (x, log_func, log_dev);

    left_bounds     = [bounds(1);  z'];
    right_bounds    = [z';  bounds(2)];
    cdf = exp(b) ./ m .* ( exp( m .* right_bounds) - exp( m .* left_bounds)); 

    samples =[];

    num_accepted = 0;
    num_sampled  = 0;

    while num_accepted < MAX_ACCEPTED

        u_star = rand * sum(cdf,1);

        % compute inverse CDF
        c = cumsum(cdf);

        ind = find(c>=u_star, 1);

        if ind == 1
            L_int = 0;
        else
            L_int = c(ind-1);
        end

        L = left_bounds(ind);

        % int^x_{L} { exp (m x + b) }   = exp(b)/m [ exp( m x) ]^x_{L}
        %                               = exp(b)/m [ exp( m x) - exp( mL) ]
        % u - L_int                     = exp(b)/m [ exp( m x) - exp( mL) ]

        u_k     = u_star - L_int;

        %                  u_k      =  exp(b)/m [ exp( m x) - exp( mL) ]
        %   exp(m x) - exp(m L )    = u_k * m / exp(b)
        %               exp(m x)    = u_k * m / exp(b) + exp(m L )

        if m(ind) <0

            % there is log(-ve) problem, so we re-write m => - m, and m >0
            %           exp(- mx)   = - u_k * m /exp(b) + exp(- m L)
            %               -mx     = log [- u_k * m /exp(b) + exp(-mL) ]
            %                mx     = log {[- u_k * m /exp(b) + exp(-mL) ]^(-1)}
            %                 x     = log {[- u_k * m /exp(b) + exp(-mL) ]^(-1)}/m

            m_pos   = - m(ind);   
            x_star  = log ((- u_k * m_pos /exp(b(ind)) + exp(-m_pos * L) )^(-1))/m_pos;

        else

            %                   x    = log [u_k * m / exp(b) + exp(m L )]/m

            %x_star = ( log(u+K) + log(m(ind)) - b(ind))/ m(ind);
            try
            x_star = log (u_k * m(ind) / exp(b(ind)) + exp(m(ind) * L ))/m(ind);
            catch
               h = 1; 
            end
        end

        u_accept = rand;

        ratio = func(x_star) / exp(m(ind)*x_star + b(ind));

        if u_accept < ratio
           samples(end+1) = x_star; 
           num_accepted = num_accepted + 1;
        else
            [m b x z f ] = get_piecewise_exp ([x;x_star], log_func, log_dev);
            left_bounds     = [-inf; z'];
            right_bounds    = [z';  inf];
            cdf = exp(b) ./ m .* ( exp( m .* right_bounds) - exp( m .* left_bounds)); 
        end

        num_sampled = num_sampled + 1;
        
    end
    
    % ----------------------------------------------
    % Debug - plot the final envelop
    % ----------------------------------------------
    
    
    f1 = figure(1);
    set(f1,'Name','adapative rejection sampling');
    
    if IS_PLOT
        
        step = (plot_bounds(2) - plot_bounds(1))/100;
        t = plot_bounds(1):step:plot_bounds(2);

        d = log_dev(x);
        for i=1:size(x,1)-1
            f_z(i)  = (z(i) - x(i)) * d(i) + f(i);
        end

        subplot(1,2,1)

        y = log_func(t);
        plot(t,y);

        hold on;

        plot(x,f,'o');

        hold on;

        plot(z,f_z,'o','color',[1 0 0]);

        hold on;
        
        border_x(1) =  plot_bounds(1);
        border_y(1) = f(1) - (x(1) - border_x(1))*d(1);
        
        border_x(2) =  plot_bounds(2);
        border_y(2) = f(end) - ( x(end) - border_x(2)) * d(end);
        
%         border_y(1) = -18;
%         border_x(1) = (border_y(1) - f(1))/d(1) + x(1);
%         border_y(2) = -18;
%         border_x(2) = (border_y(end) - f(end))/d(end) + x(end);

        plot(border_x,border_y,'o','color',[1 0 1]);

        hold on;

        l_line(1,:) = [border_x(1),border_y(1)];
        l_line(2:2+size(x,1)-2,1) = z';
        l_line(2:2+size(x,1)-2,2) = f_z;
        l_line(end+1,:) = [border_x(2),border_y(2)];

        plot ( l_line(:,1), l_line(:,2), 'Color', [1 0 1], 'LineWidth',1);

        % for i = 1: size(l_line,1)-1
        %     plot ( l_line(i:i+1,1), l_line(i:i+1,2), 'Color', [0 0 1], 'LineWidth',1);
        %     hold on;
        % end

        hold off;
        title('log(P(x)) space')
        
        

        % ---------------------------------------------
        % plot q(.) in the original scale
        % ---------------------------------------------

        subplot(1,2,2)

        y = func(t);
        plot(t,y);

        hold on;
        plot(x,func(x),'o');
        hold on;

        max_f = max(func(t));
        min_f = min(func(t));

        for i=1:size(x,1)-1
          plot([z(i) z(i)], [max_f min_f],'LineWidth',1, 'Color', [ 0 1 1]);
        end
        hold on;

        for i=1:size(x,1)
            if i == 1
                t = border_x(1):0.01:z(i);
            elseif i == size(x,1)
                t = z(end):0.01:border_x(2);
            else
                t = z(i-1):0.01:z(i);
            end

            plot(t, exp(m(i)*t + b(i)), 'color',[1 0 1]);
            hold on;
        end

        hold on;

        plot(samples, zeros(length(samples),1),'*');

        hold off;
        title('P(x) space')
    
    end
    
    
    

function [m b x z f_log] = get_piecewise_exp (x, log_func, log_dev)

    try
    f_log = log_func(x); 
    catch
        h = 1;
    end
    
    d = log_dev(x);

    [u v] = sort(x);

    temp = f_log(v);
    f_log = temp;

    temp = d(v);
    d = temp;
    
    temp = x(v);
    x = temp;

    % get the meeting points
    for i=1:size(x,1)-1
        z(i)    = ( f_log(i+1) - f_log(i) - x(i+1) * d(i+1) + x(i) * d(i) )/(  d(i) - d(i+1) );
        %f_log_z(i)  = (z(i) - x(i)) * d(i) + f_log(i);
    end

    % u(x) = f_log(i) + (x - x(i)) * d(i)
    %      = f_log(i) + d(i) x - x(i) * d(i) 
    %      = [f_log(i) - x(i) * d(i)] + d(i) * x
    %      = b + m x


    for i=1:size(x,1)
        b(i) = f_log(i) - x(i) * d(i);
        m(i) = d(i);
    end

    b = b'; m = m';

