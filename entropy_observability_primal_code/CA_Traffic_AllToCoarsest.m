clc;
clear all;

% Here we assume the group renormalization is scaled by 2 each time, and
% scaling_reps is the number of RG conducted on the systems
max_Scale_Reps = 4;   % max number in the scaling series
pixel_diff_accum = zeros(max_Scale_Reps,1);
rep = 100;
corr_val = zeros(5, 4, rep);



for m = 1:5
    for k = 1:rep
        numSite = 512;
        timeLen = 2048;      %time of simulations
        J = 0.2;    % interaction coefficient
        B = 1.5;    % external field
        
        scaling_Reps = 0;  % number of scaling times
        Config_init = randi(10,numSite,1)<= 8*ones(numSite,1); % density is 0.8
        Config_init = Config_init*2-ones(numSite,1);  % changes states to be +/-1
        Config_init(7*numSite/8:numSite,1) = -1;  % set a range to be empty sites
        
        
        
        
        
        % original system
        [Config_Orig, time_series_orig] = Stat_CA_Traffic_withField_Dec_IC(J, B, 0, numSite, timeLen, Config_init);
        scale = 2;
        
        
        i = 1;
        
        Config_pre = Config_init;
        Scale_Reps = max_Scale_Reps;
        j = 1;
        
        Config_true_scale = imresize(Config_Orig,0.5^4);
        
        
        subplot(max_Scale_Reps+1,1,1);
        imagesc(Config_true_scale);
        set(gca,'YDir','normal');
        colormap(gray)
        xlabel('Time(1 s/site)');ylabel('Position(5 m/site)');
        
        Scale_Reps = Scale_Reps+m-1;
        while(Scale_Reps > 0)
            
            Config_init = IC_Decimation(Config_pre);
            Config_pre = Config_init;
            
            [Config_RU, time_series] = Stat_CA_Traffic_withField_Dec_IC(J, B, i, numSite, timeLen, Config_init);
            Config_Reduced = imresize(Config_RU,0.5^(3+m-j));
            
            corr_val(m, i, k) = corr2(Config_true_scale, Config_Reduced);
            Scale_Reps = Scale_Reps - 1;
            j = j+1;
            i = i+1;
            subplot(max_Scale_Reps+1,1,i);
            imagesc(Config_Reduced);
            set(gca,'YDir','normal');
            colormap(gray)
            xlabel(['Time(' num2str(2^i/4) ' s/site)']);ylabel(['Position(' num2str(2^i*5/4) ' m/site)']);
        end
        
        
    end
    
end




% plot(1:4,pixel_diff_accum(:,1)');
% xlabel('Scale of RG');
% ylabel('Relative error in Norm 2');
% grid on;



%% ===================================================



%% ========================================================================

% Try to use downsampling in x and y direction respectively to reduce the
% matrix size of original system, making they are comparable

% scale = 2;
%
% Config_pst_down = zeros(numSite/scale,timeLen,max_Scale_Reps+1);
% Config_pst_down2 = zeros(numSite/scale,timeLen/scale,max_Scale_Reps+1);
%
% figure;
% subplot(max_Scale_Reps+1,1,1);
% Config_pst_down(:,:,1) = -downsample(Config_RU,scale);
% Config_pst_down2(:,:,1) = downsample(Config_pst_down(:,:,1)',scale)';
%
% i = 2;
% imagesc(Config_pst_down2(:,:,1));
% set(gca,'YDir','normal');
% colormap(gray)
% xlabel('Time(1 s/site)');ylabel('Position(5 m/site)');
% i = 2;
% Scale_Reps = max_Scale_Reps;
%
%
% while(Scale_Reps > 0)
%     Scale_Reps = Scale_Reps - 1;
%     subplot(max_Scale_Reps+1,1,i);
%     i = i+1;
%     Config_pst_down(:,:,i-1) = -downsample(Config_pst(:,:,i-2),scale);
%     Config_pst_down2(:,:,i-1) = downsample(Config_pst_down(:,:,i-1)',scale)';
%     imagesc(Config_pst_down2(:,:,i-1));
%     set(gca,'YDir','normal');
%     colormap(gray)
%     xlabel(['Time(' num2str(2^i/4) ' s/site)']);ylabel(['Position(' num2str(2^i*5/4) ' m/site)']);
%
% end




%% ========================================================================








%%

% corr_series = sum(corr_val, 2)/rep;
%
% stdard = std(corr_val');
%
%
% figure;
%
% errorbar(1:4, corr_series, stdard);
% plot(1:4, corr_series);
% xlabel('Number of renormalization times');
% ylabel('Correlation Coefficients');
% grid on;




corr_surface = zeros(5,4);

for i = 1:5
    corr_surface(i,:) = sum(corr_val(i,:,:), 2)/rep;
end



