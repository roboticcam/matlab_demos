%--------------------------------------------------
% visualization of stick-breaking process
%
%
% Written by Richard Xu:  YiDa.Xu@uts.edu.au 
% July 2014
%--------------------------------------------------

function dirichlet_process()


    % alpah         :   DP's concentration factor
    % K             :   number of sticks (Large K}
    % N_Iter        :   number of times G ~ DP(alpha, H) is drawn
    % init_parts    :   some initial partiion (just so partitions not too small
    % x_range       :   x range (for plotting)
    % max_p         :   max_p (maximum probability, for plotting)


    clear all;
    clc;

    alpha = 3;
    K = 50;
    N_Iter = 1000;
    
    CASE = 1;
     
    if CASE == 1
      
        %-----------------------------------------------
        % H is gaussian
        %-----------------------------------------------
        
        init_parts = -3:2:3;
        x_range = [-4 4];
        max_p = 0.7;
        H_pdf   = @normpdf;
        H_cdf   = @normcdf;
        H_rand  = @nrand_func;
    
    else

        %-----------------------------------------------
        % H is gamma
        %-----------------------------------------------

        init_parts = 0:2:8;
        x_range = [0 10];
        max_p = 0.7;
        H_pdf   = @gampdf_func;
        H_cdf   = @gamcdf_func;
        H_rand  = @gamrnd_func;

    end
    
    % ------------------------------------------------------------------
    % function starts here
    % ------------------------------------------------------------------

    % create random set of partitions  
    r_parts = rand([1 length(init_parts)-1 ]) .* ( init_parts(2:end) - init_parts(1:end-1)) + init_parts(1:end-1);

    P = length(r_parts);

    mean_weights = zeros(N_Iter, P + 1);

    % theretical mean from CDF
    p = H_cdf(r_parts);
    cdf_regions = [p 1] - [0 p];
    
    theory_mean = cdf_regions;

    % theretical variances
    theory_variance = theory_mean .* (ones(1, P+1) - theory_mean) / (alpha + 1);    


    for i = 1:N_Iter

        % draw random partitions
        plot ( repmat(r_parts, [ 2 1]), [zeros(1,P); ones(1,P)*max_p], 'LineWidth',2,'color',[1 0 0]);
        hold on;

        % plot f(theta)
        x_data = x_range(1):0.01:x_range(2);
        plot(x_data,H_pdf(x_data),'.');  
        hold on;
        
        % sample G ~ DP (alpha, H)
        
            
        [sticks_weights, thetas] = Stick_breaking_process(alpha, K, H_rand);

        num_samples = length(thetas);


        % draw G
        lineX = repmat(thetas, [1 2]);
        lineY = [zeros([1 num_samples]); sticks_weights' ];
        plot ( lineX', lineY, 'LineWidth',1,'color',[0 0 1]);
        hold on;

        for j=1:P+1
            if j == 1
                w_region = find(thetas < r_parts(j));
            elseif j == P+1
                w_region = find(thetas > r_parts(j-1)); 
            else
                w_region = find(thetas > r_parts(j-1) & thetas < r_parts(j) );
            end

            if ~isempty(w_region)
                mean_weights(i,j) = sum(sticks_weights(w_region));
            else
                mean_weights(i,j) = 0;
            end    
        end

        if ~isempty(w_region)
            mean_weights(i,j) = sum(sticks_weights(w_region));
        else
            mean_weights(i,j) = 0;
        end



            emperical_mean      = mean(mean_weights(1:i,:),1);
            if i > 1
                emperical_variance  = var(mean_weights(1:i,:),1);  
            end

            for j=1:P+1
                if j == 1
                    display_coord = (x_range(1)+ r_parts(j)       )/2;          
                elseif j == P+1
                    display_coord = (r_parts(j-1) + x_range(2)    )/2;
                else
                    display_coord = (r_parts(j-1) + r_parts(j)    )/2;
                end

                text( display_coord, max_p-0.1, num2str(theory_mean(j)      ,'%.3f'),   'color', [1 0 0]);
                text( display_coord, max_p-0.2, num2str(emperical_mean(j)   ,'%.3f'),   'color', [1 0 0]);
                text( display_coord, max_p-0.3, num2str(theory_variance(j)  ,'%.3f'),   'color', [1 0 1]);
                if i > 1
                    text( display_coord, max_p-0.4, num2str(emperical_variance(j)  ,'%.3f'),'color', [1 0 1]);
                else
                    text( display_coord, max_p-0.4, num2str(0  ,'%.3f'),'color', [1 0 1]);
                end
            end


        set(gca,'YTick',[0.5:0.1:max_p]);
        set(gca,'YTick',[max_p-0.4:0.1:max_p-0.1]);
        set(gca,'YTickLabel',[ 'emperical var '; 'theoretic var '; 'emperical mean'; 'theoretic mean']);


        hold on;

        axis([x_range(1) x_range(2) 0 max_p]);

        waitforbuttonpress;

        hold off;

        %sum(mean_weights(i,:))
    end

    % emperical mean
    emperical_mean = mean(mean_weights);
    display(emperical_mean);


    % emperical variances
    emperical_variance = var(mean_weights);    
    display(emperical_variance);

end


function probs = gampdf_func (x_data)
    a = 2; b = 2;
    probs = gampdf(x_data,a,b);
end

function probs = gamcdf_func (x_data)
    a = 2; b = 2;
    probs = gamcdf(x_data,a,b);
end

function x = gamrnd_func(N)
    a = 2; b = 2;
    x = gamrnd(a,b,[N 1]);
end

function x = nrand_func(N)
    x = randn([N 1]);
end


% stick-breaking process
% when H_rand is a (continous) density function handel, it draws a sample,
% otherwise, it performs a multinomial distribution draw

function [sticks_weights, thetas] = Stick_breaking_process(alpha, N, H_rand)

    raw_sticks = betarnd(1,alpha,[N 1]);

    sticks_weights = zeros(N,1);

    sticks_weights(1) = raw_sticks(1);

    for i = 2:N
        sticks_weights(i) = raw_sticks(i) * (1 - sum(sticks_weights(1:i-1))); 
    end

    % sampling thetas
    
    if isa(H_rand, 'function_handle')
        thetas = H_rand(N);
    else
        for t=1:N
            thetas(t) = find(rand < cumsum(H_rand),1);
        end
        h =1;    
    end
end
