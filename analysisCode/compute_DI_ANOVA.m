%% Code to compute_DI_ANOVA

% compute the following metrics for 5 stim types
% - ANOVA
% - baseline subtracted DI

load ../dataFiles/cellData_sua.mat; 

for n=1:length(cellData_sua)

    % pre-computed response matrix
    % - rows: trials
    % - cols: conditions
    respMtx = cellData_sua(n).respMtx;

    % no stim response 
    baseline = nanmean(respMtx(:,end)); 

    % stim condition: five stim types
    condNames = [{'LRM-noise'},{'LRM-sinusoid'},{'Local'},...
                 {'LRM-sinusoid-Local-same'},{'LRM-sinusoid-Local-opp'}]; 


    dir_tuning = [];    % 8 x 5 matrix
    DI_base = [];    % baseline subtracted DI

    for m=1:5  % five stim types: LRM-noise, LRM-sinusoid, Local, 
               %                  LRM-sinusoid-Local-same, 
               %                  LRM-sinusoid-Local-opp
        
        % condition numbers
        condNums = (m-1)*8+1:m*8; 

        % ANOVA
        [ap1,tbl1,stats1] = anova1(respMtx(:,condNums), [1:8],'off'); 
        if isnan(ap1)
            cellData_sua(n).anova(:,m) = [0, 1];
        else
            cellData_sua(n).anova(:,m) = [tbl1{2,5}, ap1];   % fstat, p-val
        end

        % direction tuning curve
        dir_tuning(:,m) = nanmean(respMtx(:,condNums),1); 


        % baseline subtracted DI
        abs_modulation = abs(dir_tuning(:,m) - baseline); 
        pref_id = find(abs_modulation(:)==max(abs_modulation));         
        mod_depth = []; 
        for p=1:length(pref_id) % when multiple points were found
            if pref_id(p) <= 4
                non_pref_id(p) = pref_id(p) + 4; 
            else
                non_pref_id(p) = pref_id(p) - 4; 
            end

            mod_depth(p) = abs(dir_tuning(pref_id(p),m) - dir_tuning(non_pref_id(p),m)); 
        end
        pref_id = pref_id(find(mod_depth(:)==max(mod_depth))); 
        non_pref_id = non_pref_id(find(mod_depth(:)==max(mod_depth)));         

        numerator = abs(dir_tuning(pref_id,m) - dir_tuning(non_pref_id,m)); 
        denominator = abs(dir_tuning(pref_id,m)-baseline);         
        if(size(numerator/denominator)> 1)
            DI_base(m) = numerator(1)/denominator(1);
        else
            DI_base(m) = numerator / denominator; 
        end        
    end

    % replace NaNs with appropriate values
    DI_base(isnan(DI_base))=0;
    
    cellData_sua(n).dir_tuning = dir_tuning; 
    cellData_sua(n).DI_base = DI_base; 

    clearvars -except cellData_sua; 
end

