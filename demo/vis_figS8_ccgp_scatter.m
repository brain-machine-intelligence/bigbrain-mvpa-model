
clear all

%% Set working directory

base_path = fullfile('..'); % If your current working directory is 'demo', use this. Otherwise, set your own path.
addpath(fullfile(base_path, 'functions', 'utils'))

% If your current working directory is 'demo', use this. Otherwise, set your own path that includes 'data'.
cd('..'); 


%% load setting 

% ===== experiment =====
Exp = 'Lee2014';  N_sbj = 20; Exp_sfx = [];


% ===== ROI mask for multi-voxel pattern (MVP) =====
ROI = {'F3TL', 'F3TR', 'DLPFCL', 'DLPFCR', 'OFC', 'ACCL', 'ACCR', 'preSMAL', 'preSMAR', 'V1L', 'V1R', 'HIPPOL', 'HIPPOR', 'VentStrL', 'VentStrR'};
% ROI = {'F3TL', 'F3TR', 'DLPFCL', 'DLPFCR', 'OFC', 'V1L', 'V1R', 'HIPPOL', 'HIPPOR', 'VentStrL', 'VentStrR'};


% ====================== labelling by task variable =======================
% The labels for goal CCGP across uncertainty, goal SD in low uncertainty
% (uL), and goal SD in high uncertianty (uH), respectively
LABEL = {'Goal(6,7,8)'}; N_allclass = 3;


%% load setting 2

% ========================= MVP setting =========================
% ROI_types = {};
% main set
% ROI_types = {'AAL3', 'AAL3', 'AAL', 'AAL', 'AAL', 'AAL3', 'AAL3', 'AAL3', 'AAL3', 'AAL3', 'AAL3'};
% ROI_types = {'AAL3', 'AAL3', 'AAL3', 'AAL3'};
ROI_types = {'AAL3', 'AAL3', 'AAL', 'AAL', 'AAL', 'AAL3', 'AAL3', 'JuBrain', 'JuBrain', 'AAL3', 'AAL3', 'AAL3', 'AAL3', 'AAL3', 'AAL3'};

EVENT = [1 2 4 5 7 8 9 10];      
EVENT_NAME = {  'fix ', 'S1 ', 'A1 ', 'S2 ', 'A2 ', 'S3 ', ... 
                'fix` ', 'S1` '}; % f2 = A1, f3 = A2

% sanity check
% EVENT_map = containers.Map(1:8, EVENT_NAME);
if length(EVENT) ~= length(EVENT_NAME); error('EVENT_NAME error'); end

% strict = 0;
% interpolation = 0;


% ========================= measure =========================

% sanity_verbose = 1;
sanity_verbose = 0;

% --- shattering --- 
% measure = 'accs';
% measure = 'accs_r';
% measure = 'separable';
% measure = 'separable2';
% measure = 'separable2_r';
% measure = 't-stat';
measure = 'acc_mean';
% measure = 'PR_mean';
% measure = 'nPC_mean';
% measure = 'sample unbal';
% measure = 'sample_sizes';
% measure = 'diff acc_mean acc_mean';
% measure = 'cosine';
% measure = 'corr_sig';   % significance of weight vector correlation

% multiple measures (mainly for 'mean_msr_box' only, not for 'indiv_msr_vec')
% measure_cell = {'CCGP', 'acc_mean', 'acc_mean', 'acc_mean'}; 
% measure_cell = {'CCGP', 'acc_mean', 'acc_mean'}; 
measure_cell = {};

% For shattering CCGP
% CClabel_name = [];
% CClabel_name = 'CmplxCond';
CClabel_name = 'UncCond';
% CClabel_name = 'Goal(6,7,8)';



% ========================= save setting =========================
% --- shattering ---
save_path = 'data/ccgp_result';

if isempty(measure_cell)
    if contains(measure, 'diff')
        save_name = {['fitclinear_accbox' Exp_sfx], ['fitclinear_accbox_' CClabel_name 'CCGP' Exp_sfx]};
        measure_cell = {};
    elseif contains(measure, 'cosine') || contains(measure, 'corr_sig')
%         save_name = ['fitclinear_weights' Exp_sfx];
        save_name = ['fitclinear_accbox' Exp_sfx];
    else
        if isempty(CClabel_name) % SD
%             save_name = ['fitclinear_weights' Exp_sfx];
            save_name = ['fitclinear_accbox' Exp_sfx];
        else % CCGP
            save_name = ['fitclinear_accbox_' CClabel_name 'CCGP' Exp_sfx];
        end
        measure_cell = {};
    end
else
    save_name = cell(1, length(measure_cell));
    for mi = 1:length(measure_cell)
        temp_msr = measure_cell{mi};
        switch temp_msr
            case 'acc_mean'
                save_name{mi} = ['fitclinear_accbox' Exp_sfx];
            case 'CCGP'
                save_name{mi} = ['fitclinear_accbox_' CClabel_name 'CCGP' Exp_sfx];
        end
    end
end
% save_name = 'fitclinear_accbox_percent';
% save_name = 'fitclinear_accbox_strict';
% save_name = 'fitclinear_accbox_Heo';
% save_name = 'fitclinear_accbox_Kim';
% save_name = 'fitclinear_accbox_Kim_wa';




if contains(save_path, 'ccgp')
    chance_lvl = 0.5 * ones(size(N_allclass));
elseif contains(save_path, 'decoding')
    if contains(save_name, 'fitrlinear')
        chance_lvl = 0./N_allclass;
    else
        chance_lvl = 1./N_allclass;
    end
end


if isempty(measure_cell) && (contains(measure, 'diff') || contains(measure, 'cosine'))
    chance_lvl = zeros(size(N_allclass));
end



disp('step1 done')

%% data load

% newrun = 0;     % including symmetric duplicated dichotomy
newrun = 1;   % ignoring symmetric duplicated dichotomy

% multiple measures operation (mainly for 'mean_msr_box' only, not for 'indiv_msr_vec')
if isempty(measure_cell) && contains(measure, 'diff')
    measure_cell = split(measure);
    msr_operation = measure_cell{1};
    measure_cell = measure_cell(2:end);
    measure_ori = measure;
elseif isempty(measure_cell)
    measure_cell = {measure};
    msr_operation = [];
else
    msr_operation = [];
end
mean_msr_cell = cell(1, length(measure_cell)); 

% multiple save_name handling
if iscell(save_name); name_cell = save_name; else; name_cell = {save_name}; end

% measure loop
for mi = 1:length(measure_cell)
    measure = measure_cell{mi}; disp(' '); fprintf([measure ' '])
    save_name = name_cell{mi};
    
    % initializations
    mean_msr_box = nan * ones(N_sbj, length(LABEL), length(ROI), length(EVENT));
    if ~newrun
        vec_len = 2^max(N_allclass)-2;
    else
        vec_len = (2^max(N_allclass)-2)/2;
    end
    if contains(measure, 'cosine'); vec_len = 3000; end
    indiv_msr_vec = nan * ones(N_sbj, length(LABEL), length(ROI), length(EVENT), vec_len);
    posi_labs_box = cell(N_sbj, length(LABEL), length(ROI), length(EVENT));
    
    % main loading
    for id = 1:N_sbj
        
        fprintf('%d ', id)
        
        for roii = 1:length(ROI)

            if ~isempty(ROI_types)
                roi_type = ROI_types{roii};
            end
            
            for labi = 1:length(LABEL)
                
                % 2021-11-20 temp
                if length(LABEL) == length(measure_cell) && (labi ~= mi)
                    continue
                end
                
                for evi = 1:length(EVENT)
                    
                    % event name for saving
                    if iscell(EVENT)
                        ev_name = EVENT{evi};
                    else
                        ev_name = EVENT(evi);
                    end
                    

                    try
                        data_load = load([save_path '/' roi_type '/' ROI{roii} '/' , ...
                            LABEL{labi} '/' num2str(ev_name) '/', ...
                            save_name num2str(id) '.mat']);

                        % sanity check
                        invalid = sum(data_load.save_info.acc_box. ...
                            ShatteredClassesNumber{N_allclass(labi)}.posi_or_neg_lab_only) ...
                            + sum(data_load.save_info.acc_box. ...
                            ShatteredClassesNumber{N_allclass(labi)}.insufficient_for_CV);

                        if invalid
                            disp('invalid')

                        else
                            if strcmp(measure, 't-stat')
                                acc_means = data_load.save_info.acc_box. ...
                                    ShatteredClassesNumber{N_allclass(labi)}.accs;
                                acc_stes = data_load.save_info.acc_box. ...
                                    ShatteredClassesNumber{N_allclass(labi)}.stes;
                                acc_ts = (acc_means-0.5)./acc_stes;
                                if any(isnan(acc_ts)); error('tstat error'); end
                                mean_msr_box(id, labi, roii, evi) = mean(acc_ts);
                                indiv_msr_vec(id, labi, roii, evi, 1:2^N_allclass(labi)-2) = acc_ts;

                            elseif strcmp(measure, 'PR_mean')
                                msr_vec = data_load.save_info.acc_box. ...
                                    ShatteredClassesNumber{N_allclass(labi)}.PRs;
                                msr_vec = squeeze(mean(msr_vec, 3));
                                mean_msr_box(id, labi, roii, evi) = mean(msr_vec);
                                if ~newrun
                                    indiv_msr_vec(id, labi, roii, evi, 1:2^N_allclass(labi)-2) = msr_vec;
                                else
                                    indiv_msr_vec(id, labi, roii, evi, 1:(2^N_allclass(labi)-2)/2) = msr_vec;
                                end

                            elseif strcmp(measure, 'nPC_mean')
                                msr_vec = data_load.save_info.acc_box. ...
                                    ShatteredClassesNumber{N_allclass(labi)}.nPCs;
                                msr_vec = squeeze(mean(msr_vec, 3));
                                mean_msr_box(id, labi, roii, evi) = mean(msr_vec);
                                if ~newrun
                                    indiv_msr_vec(id, labi, roii, evi, 1:2^N_allclass(labi)-2) = msr_vec;
                                else
                                    indiv_msr_vec(id, labi, roii, evi, 1:(2^N_allclass(labi)-2)/2) = msr_vec;
                                end

                            elseif strcmp(measure, 'sample unbal')
                                msr_vec = data_load.save_info.acc_box. ...
                                    ShatteredClassesNumber{N_allclass(labi)}.sample_sizes;
                                for i = 1:length(msr_vec); msr_vec{i} = abs(diff(msr_vec{i})); end
                                msr_vec = cell2mat(msr_vec);
                                mean_msr_box(id, labi, roii, evi) = mean(msr_vec);
                                if ~newrun
                                    indiv_msr_vec(id, labi, roii, evi, 1:2^N_allclass(labi)-2) = msr_vec;
                                else
                                    indiv_msr_vec(id, labi, roii, evi, 1:(2^N_allclass(labi)-2)/2) = msr_vec;
                                end

                            elseif strcmp(measure, 'accs_r')
                                msr_vec = data_load.save_info.acc_box. ...
                                    ShatteredClassesNumber{N_allclass(labi)}.accs_r;
                                msr_vec = squeeze(mean(msr_vec, 3));
                                mean_msr_box(id, labi, roii, evi) = mean(msr_vec);
                                indiv_msr_vec(id, labi, roii, evi, 1:2^N_allclass(labi)-2) = msr_vec;

                            elseif strcmp(measure, 'separable2_r')
                                accs_r = data_load.save_info.acc_box. ...
                                    ShatteredClassesNumber{N_allclass(labi)}.accs_r;
                                msr_vec = zeros(1, 2^N_allclass(labi)-2);
                                for i = 1:size(accs_r, 2)
                                    accs_r_ = squeeze(accs_r(:,i,:));
                                    msr_vec(i) = ttest(accs_r_(~isnan(accs_r_)), 0.5 * ones(sum(~isnan(accs_r_)), 1), ...
                                        'Tail', 'right', 'Alpha', 0.0001);
                                end
                                mean_msr_box(id, labi, roii, evi) = mean(msr_vec);
                                indiv_msr_vec(id, labi, roii, evi, 1:2^N_allclass(labi)-2) = msr_vec;

                            elseif strcmp(measure, 'CCGP')
                                msr_vec = data_load.save_info.acc_box. ...
                                    ShatteredClassesNumber{N_allclass(labi)}.acc_mean;
                                mean_msr_box(id, labi, roii, evi) = mean(msr_vec);
                                if ~newrun
                                    indiv_msr_vec(id, labi, roii, evi, 1:2^N_allclass(labi)-2) = msr_vec;
                                else
                                    indiv_msr_vec(id, labi, roii, evi, 1:(2^N_allclass(labi)-2)/2) = msr_vec;
                                end

                            else
                                %% ===== load shattering measure =====
                                % msr_mat: acc. or sep. of all possible classification
                                % size: (N_allclass choose N_shatteredclass) ...
                                % x (2^N_shatteredclass - 2)
                                msr_vec = data_load.save_info.acc_box. ...
                                    ShatteredClassesNumber{N_allclass(labi)}.(measure);
                                % ***** original *****
                                mean_msr_box(id, labi, roii, evi) = mean(msr_vec);

                                if ~newrun
                                    indiv_msr_vec(id, labi, roii, evi, 1:2^N_allclass(labi)-2) = msr_vec;
                                else
                                    indiv_msr_vec(id, labi, roii, evi, 1:(2^N_allclass(labi)-2)/2) = msr_vec;
                                end
                            end
                            posi_labs_box{id, labi, roii, evi} = data_load. ...
                                save_info.acc_box. ...
                                ShatteredClassesNumber{N_allclass(labi)}. ...
                                positive_labels;
                        end

                    catch Err
                        disp([num2str(Err.stack.line) ' ' Err.message])
                        %                     disp('load problem')

                    end


                end % evi
            end % labi
            %         disp(' ' )  % 2021-03-25 temporary
        end % roii
    end % id
    
    mean_msr_cell{mi} = mean_msr_box;
    
end % measure

% ===== measure manipulation =====
if length(mean_msr_cell) > 1
    if isempty(msr_operation)
        % 2021-11-20 temp
        for tlabi = 1:length(LABEL)
            temp_mean_msr_box = mean_msr_cell{tlabi};
            mean_msr_box(:, tlabi, :, :) = temp_mean_msr_box(:, tlabi, :, :);
        end
        measure = 'acc_mean';
    elseif contains(msr_operation, 'diff')
        mean_msr_box = mean_msr_cell{1} - mean_msr_cell{2};
        measure = measure_ori;
    elseif contains(msr_operation, 'ratio')
        mean_msr_box = mean_msr_cell{1} / mean_msr_cell{2};
        measure = measure_ori;
    end
end

% ===== running information display =====
disp(' ')
disp(data_load.run_info)
disp(data_load.run_info.MVP_setting)
disp(data_load.run_info.stat_setting)

disp(' ')

disp('step2 done')

% for below step (save original mean_msr_box)
temp_msr_means = mean_msr_box; 
temp_msr_vecs = indiv_msr_vec; 
original_ROI = ROI;
ROI_modif = 0;

%% measure box reorganization

% +++++++++++++++++++++++ for specific clsf result ++++++++++++++++++++++++

% bilateral_integ = 0;
bilateral_integ = 1;

indiv_clsf_result = 0;
% indiv_clsf_result = 1;
% indiv_clsf_result = 2;
% indiv_clsf_result = 3;
% indiv_clsf_result = 4;

% --- blkCond clsf types ---
switch Exp
    case {'Lee2014', 'Heo2018'}
        if newrun
            interaction_clsf = 5;
            goalcond_clsf = 3;
            unccond_clsf = 6;
            other_clsf = setdiff(1:7, [3 5 6]);
            other1 = 1;
            other2 = 2;
            other3 = 4;
            other4 = 7;
        else
            interaction_clsf = [5 10];
            goalcond_clsf = [3 12];
            unccond_clsf = [6 9];
            other_clsf = setdiff(1:14, [3 5 6 9 10 12]);
            other1 = [1 14];
            other2 = [2 13];
            other3 = [4 11];
            other4 = [7 8];
        end
    case 'Kim2019' % assuming always newrun
        interaction_clsf = 6;
        goalcond_clsf = 5;      % in fact, complexity condition
        unccond_clsf = 3;
        other_clsf = setdiff(1:7, [3 5 6]);
        other1 = 1;
        other2 = 2;
        other3 = 4;
        other4 = 7;
end

temp_msr_box = temp_msr_means;     % original mean_msr_box


% temp_msr_box = squeeze(mean(temp_msr_vecs(:,:,:,:,interaction_clsf),5));   % [2 4] vs [1 3] blkCond
% temp_msr_box = squeeze(mean(temp_msr_vecs(:,:,:,:,goalcond_clsf),5));   % [3 4] vs [1 2] blkCond
% temp_msr_box = squeeze(mean(temp_msr_vecs(:,:,:,:,unccond_clsf),5));   % [1 4] vs [2 3] blkCond

% temp_msr_box = squeeze(mean(temp_msr_vecs(:,:,:,:,other1),5));   
% temp_msr_box = squeeze(mean(temp_msr_vecs(:,:,:,:,other2),5));   
% temp_msr_box = squeeze(mean(temp_msr_vecs(:,:,:,:,other3),5));   
% temp_msr_box = squeeze(mean(temp_msr_vecs(:,:,:,:,other4),5));   
% temp_msr_box = squeeze(mean(temp_msr_vecs(:,:,:,:,other_clsf),5));   

% temp_msr_box = squeeze(mean(temp_msr_vecs(:,:,:,:,[goalcond_clsf, unccond_clsf]),5));
% temp_msr_box = squeeze(mean(temp_msr_vecs(:,:,:,:,[goalcond_clsf, unccond_clsf, other_clsf]),5));
% temp_msr_box = squeeze(mean(temp_msr_vecs(:,:,:,:,[interaction_clsf, other_clsf]),5));

if indiv_clsf_result
    if length(EVENT)>1
        temp_msr_vec_lin = squeeze(mean(temp_msr_vecs(:,:,:,:,[goalcond_clsf, unccond_clsf, other_clsf]),5));
        temp_msr_vec_non = squeeze(mean(temp_msr_vecs(:,:,:,:,interaction_clsf),5));
        temp_msr_vec_main = squeeze(mean(temp_msr_vecs(:,:,:,:,[goalcond_clsf, unccond_clsf]),5));
        temp_msr_vec_other = squeeze(mean(temp_msr_vecs(:,:,:,:,other_clsf),5));
    else
        temp_msr_vec_lin = reshape(mean(temp_msr_vecs(:,:,:,:,[goalcond_clsf, unccond_clsf, other_clsf]),5), N_sbj, length(LABEL), length(ROI), length(EVENT));
        temp_msr_vec_non = reshape(mean(temp_msr_vecs(:,:,:,:,interaction_clsf),5), N_sbj, length(LABEL), length(ROI), length(EVENT));
        temp_msr_vec_main = reshape(mean(temp_msr_vecs(:,:,:,:,[goalcond_clsf, unccond_clsf]),5), N_sbj, length(LABEL), length(ROI), length(EVENT));
        temp_msr_vec_other = reshape(mean(temp_msr_vecs(:,:,:,:,other_clsf),5), N_sbj, length(LABEL), length(ROI), length(EVENT));
    end
end

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% for combining bilateral result

if length(original_ROI) == 8
    % F3T, V1, HPC, vStr
    roii1 = [1 2];
    % roii1 = 1:17;
    roii2 = [3 4];
    roii3 = [5 6];
    roii4 = [7 8];
    bilateral_mode = 1;
elseif length(original_ROI) == 11
    % % +DLPFC & R/L integrated OFC
    roii1 = [1 2];
    roii2 = [3 4];
    roii3 = 5;
    roii4 = [6 7];
    roii5 = [8 9];
    roii6 = [10 11];
    bilateral_mode = 2;
elseif length(original_ROI) == 5
    roii1 = [1 2]; % vlPFC
    roii2 = [3 4]; % dlPFC
    roii3 = 5; % OFC
    bilateral_mode = 3;
elseif length(original_ROI) == 7
    roii1 = [1 2]; % vlPFC
    roii2 = [3 4]; % dlPFC
    roii3 = 5; % OFC
    roii4 = [6 7]; % V1
    bilateral_mode = 4;
elseif length(original_ROI) == 15
    roii1 = [1 2]; % vlPFC
    roii2 = [3 4]; % dlPFC
    roii3 = 5; % OFC
    roii4 = [6 7]; % ACC
    roii5 = [8 9]; % SMA
    roii6 = [10 11]; % V1
    roii7 = [12 13]; % HPC
    roii8 = [14 15]; % vStr
    bilateral_mode = 5;
end

% clear mean_msr_box indiv_msr_vec

if bilateral_integ
    
    % ===== for bilateral averaging =====
    mean_msr_box = NaN(size(mean_msr_box));
    indiv_msr_vec = NaN(size(indiv_msr_vec));
    
    switch bilateral_mode
        case 1 % F3T, V1, HPC, vStr
            mean_msr_box(:,:,1,:) = mean(temp_msr_box(:, :, roii1, :), 3);
            mean_msr_box(:,:,2,:) = mean(temp_msr_box(:, :, roii2, :), 3);
            mean_msr_box(:,:,3,:) = mean(temp_msr_box(:, :, roii3, :), 3);
            mean_msr_box(:,:,4,:) = mean(temp_msr_box(:, :, roii4, :), 3);
            ROI = {'vlPFC', 'V1', 'HPC', 'vStr'};
            %     ROI = {'ilPFC', 'V1'};
            %     ROI = {'ilPFC', 'V1', 'HIPPO', 'VentStr'};
            %     ROI = {'vlPFC', 'V1', 'dlPFC', 'OFC'};
            %     ROI = {'PALL', 'Put', 'Cd', 'VentStr'};
            
            indiv_msr_vec(:,:,1,:,:) = mean(temp_msr_vecs(:, :, roii1, :, :), 3);
            indiv_msr_vec(:,:,2,:,:) = mean(temp_msr_vecs(:, :, roii2, :, :), 3);
            indiv_msr_vec(:,:,3,:,:) = mean(temp_msr_vecs(:, :, roii3, :, :), 3);
            indiv_msr_vec(:,:,4,:,:) = mean(temp_msr_vecs(:, :, roii4, :, :), 3);
            
        case 2 % F3T, V1, HPC, vStr +DLPFC & R/L integrated OFC
            mean_msr_box(:,:,1,:) = mean(temp_msr_box(:, :, roii1, :), 3);
            mean_msr_box(:,:,2,:) = mean(temp_msr_box(:, :, roii2, :), 3);
            mean_msr_box(:,:,3,:) = mean(temp_msr_box(:, :, roii3, :), 3);
            mean_msr_box(:,:,4,:) = mean(temp_msr_box(:, :, roii4, :), 3);
            mean_msr_box(:,:,5,:) = mean(temp_msr_box(:, :, roii5, :), 3); % +DLPFC & R/L integrated OFC
            mean_msr_box(:,:,6,:) = mean(temp_msr_box(:, :, roii6, :), 3); % +DLPFC & R/L integrated OFC
%             ROI = {'vlPFC', 'dlPFC', 'OFC', 'V1', 'HIPPO', 'VentStr'}; % +DLPFC & R/L integrated OFC
            ROI = {'vlPFC', 'dlPFC', 'OFC', 'V1', 'HPC', 'vStr'}; % +DLPFC & R/L integrated OFC
            
            indiv_msr_vec(:,:,1,:,:) = mean(temp_msr_vecs(:, :, roii1, :, :), 3);
            indiv_msr_vec(:,:,2,:,:) = mean(temp_msr_vecs(:, :, roii2, :, :), 3);
            indiv_msr_vec(:,:,3,:,:) = mean(temp_msr_vecs(:, :, roii3, :, :), 3);
            indiv_msr_vec(:,:,4,:,:) = mean(temp_msr_vecs(:, :, roii4, :, :), 3);
            indiv_msr_vec(:,:,5,:,:) = mean(temp_msr_vecs(:, :, roii5, :, :), 3); % +DLPFC & R/L integrated OFC
            indiv_msr_vec(:,:,6,:,:) = mean(temp_msr_vecs(:, :, roii6, :, :), 3); % +DLPFC & R/L integrated OFC
    
        case 3 % F3T (vlPFC), DLPFC, R/L integrated OFC
            mean_msr_box(:,:,1,:) = mean(temp_msr_box(:, :, roii1, :), 3);
            mean_msr_box(:,:,2,:) = mean(temp_msr_box(:, :, roii2, :), 3);
            mean_msr_box(:,:,3,:) = mean(temp_msr_box(:, :, roii3, :), 3);
            ROI = {'vlPFC', 'dlPFC', 'OFC'};
            
            indiv_msr_vec(:,:,1,:,:) = mean(temp_msr_vecs(:, :, roii1, :, :), 3);
            indiv_msr_vec(:,:,2,:,:) = mean(temp_msr_vecs(:, :, roii2, :, :), 3);
            indiv_msr_vec(:,:,3,:,:) = mean(temp_msr_vecs(:, :, roii3, :, :), 3);
            
        case 4 % F3T (vlPFC), DLPFC, R/L integrated OFC, V1
            mean_msr_box(:,:,1,:) = mean(temp_msr_box(:, :, roii1, :), 3);
            mean_msr_box(:,:,2,:) = mean(temp_msr_box(:, :, roii2, :), 3);
            mean_msr_box(:,:,3,:) = mean(temp_msr_box(:, :, roii3, :), 3);
            mean_msr_box(:,:,4,:) = mean(temp_msr_box(:, :, roii4, :), 3);
            ROI = {'vlPFC', 'dlPFC', 'OFC', 'V1'};
            
            indiv_msr_vec(:,:,1,:,:) = mean(temp_msr_vecs(:, :, roii1, :, :), 3);
            indiv_msr_vec(:,:,2,:,:) = mean(temp_msr_vecs(:, :, roii2, :, :), 3);
            indiv_msr_vec(:,:,3,:,:) = mean(temp_msr_vecs(:, :, roii3, :, :), 3);
            indiv_msr_vec(:,:,4,:,:) = mean(temp_msr_vecs(:, :, roii4, :, :), 3);

        case 5
            ROI = {'vlPFC', 'dlPFC', 'OFC', 'ACC', 'preSMA', 'V1', 'HPC', 'vStr'};

            mean_msr_box(:,:,1,:) = mean(temp_msr_box(:, :, roii1, :), 3);
            mean_msr_box(:,:,2,:) = mean(temp_msr_box(:, :, roii2, :), 3);
            mean_msr_box(:,:,3,:) = mean(temp_msr_box(:, :, roii3, :), 3);
            mean_msr_box(:,:,4,:) = mean(temp_msr_box(:, :, roii4, :), 3);
            mean_msr_box(:,:,5,:) = mean(temp_msr_box(:, :, roii5, :), 3);
            mean_msr_box(:,:,6,:) = mean(temp_msr_box(:, :, roii6, :), 3);
            mean_msr_box(:,:,7,:) = mean(temp_msr_box(:, :, roii7, :), 3);
            mean_msr_box(:,:,8,:) = mean(temp_msr_box(:, :, roii8, :), 3);
            
            indiv_msr_vec(:,:,1,:,:) = mean(temp_msr_vecs(:, :, roii1, :, :), 3);
            indiv_msr_vec(:,:,2,:,:) = mean(temp_msr_vecs(:, :, roii2, :, :), 3);
            indiv_msr_vec(:,:,3,:,:) = mean(temp_msr_vecs(:, :, roii3, :, :), 3);
            indiv_msr_vec(:,:,4,:,:) = mean(temp_msr_vecs(:, :, roii4, :, :), 3);
            indiv_msr_vec(:,:,5,:,:) = mean(temp_msr_vecs(:, :, roii5, :, :), 3);
            indiv_msr_vec(:,:,6,:,:) = mean(temp_msr_vecs(:, :, roii6, :, :), 3);
            indiv_msr_vec(:,:,7,:,:) = mean(temp_msr_vecs(:, :, roii7, :, :), 3);
            indiv_msr_vec(:,:,8,:,:) = mean(temp_msr_vecs(:, :, roii8, :, :), 3);
            
    end
    
else
    ROI = original_ROI;
end

switch indiv_clsf_result
    % ========== for indiv clsf ==========
    
    case 1
        if bilateral_integ
            mean_msr_box(:,:,1,:) = mean(temp_msr_vec_main(:, :, roii1, :), 3);
            mean_msr_box(:,:,2,:) = mean(temp_msr_vec_other(:, :, roii1, :), 3);
            mean_msr_box(:,:,3,:) = mean(temp_msr_vec_non(:, :, roii1, :), 3);
            ROI = {'ilPFC main', 'ilPFC one vs others', 'ilPFC nonlin'}; ROI_modif = 1; color_id = 1;
        else
%             mean_msr_box(:,:,1,:) = mean(temp_msr_vec_main(:, :, 1, :), 3);
%             mean_msr_box(:,:,2,:) = mean(temp_msr_vec_other(:, :, 1, :), 3);
%             mean_msr_box(:,:,3,:) = mean(temp_msr_vec_non(:, :, 1, :), 3);
%             mean_msr_box(:,:,4,:) = mean(temp_msr_vec_main(:, :, 2, :), 3);
%             mean_msr_box(:,:,5,:) = mean(temp_msr_vec_other(:, :, 2, :), 3);
%             mean_msr_box(:,:,6,:) = mean(temp_msr_vec_non(:, :, 2, :), 3);
%             ROI = {'F3TL main', 'F3TL one vs others', 'F3TL nonlin', ...
%                    'F3TR main', 'F3TR one vs others', 'F3TR nonlin'}; 
            mean_msr_box(:,:,1,:) = mean(temp_msr_vec_main(:, :, 1, :), 3);
            mean_msr_box(:,:,2,:) = mean(temp_msr_vec_other(:, :, 1, :), 3);
            mean_msr_box(:,:,3,:) = mean(temp_msr_vec_non(:, :, 1, :), 3);
            ROI = {'F3TLus1 main', 'F3TLus1 one vs others', 'F3TLus1 nonlin'}; 
            ROI_modif = 1; color_id = 1;
        end
       
    case 2
        if bilateral_integ
            mean_msr_box(:,:,1,:) = mean(temp_msr_vec_main(:, :, roii2, :), 3);
            mean_msr_box(:,:,2,:) = mean(temp_msr_vec_other(:, :, roii2, :), 3);
            mean_msr_box(:,:,3,:) = mean(temp_msr_vec_non(:, :, roii2, :), 3);
            ROI = {'V1 main', 'V1 one vs others', 'V1 nonlin'}; ROI_modif = 1; color_id = 2;
        else
%             mean_msr_box(:,:,1,:) = mean(temp_msr_vec_main(:, :, 3, :), 3);
%             mean_msr_box(:,:,2,:) = mean(temp_msr_vec_other(:, :, 3, :), 3);
%             mean_msr_box(:,:,3,:) = mean(temp_msr_vec_non(:, :, 3, :), 3);
%             mean_msr_box(:,:,4,:) = mean(temp_msr_vec_main(:, :, 4, :), 3);
%             mean_msr_box(:,:,5,:) = mean(temp_msr_vec_other(:, :, 4, :), 3);
%             mean_msr_box(:,:,6,:) = mean(temp_msr_vec_non(:, :, 4, :), 3);
%             ROI = {'V1L main', 'V1L one vs others', 'V1L nonlin', ...
%                    'V1R main', 'V1R one vs others', 'V1R nonlin'}; 
            mean_msr_box(:,:,1,:) = mean(temp_msr_vec_main(:, :, 2, :), 3);
            mean_msr_box(:,:,2,:) = mean(temp_msr_vec_other(:, :, 2, :), 3);
            mean_msr_box(:,:,3,:) = mean(temp_msr_vec_non(:, :, 2, :), 3);
            ROI = {'V1Lus1 main', 'V1Lus1 one vs others', 'V1Lus1 nonlin'}; 
            ROI_modif = 1; color_id = 2;
        end
        
    case 3
        if bilateral_integ
            mean_msr_box(:,:,1,:) = mean(temp_msr_vec_main(:, :, roii3, :), 3);
            mean_msr_box(:,:,2,:) = mean(temp_msr_vec_other(:, :, roii3, :), 3);
            mean_msr_box(:,:,3,:) = mean(temp_msr_vec_non(:, :, roii3, :), 3);
            ROI = {'HIPPO main', 'HIPPO one vs others', 'HIPPO nonlin'}; ROI_modif = 1; color_id = 3;
        else
%             mean_msr_box(:,:,1,:) = mean(temp_msr_vec_main(:, :, 5, :), 3);
%             mean_msr_box(:,:,2,:) = mean(temp_msr_vec_other(:, :, 5, :), 3);
%             mean_msr_box(:,:,3,:) = mean(temp_msr_vec_non(:, :, 5, :), 3);
%             mean_msr_box(:,:,4,:) = mean(temp_msr_vec_main(:, :, 6, :), 3);
%             mean_msr_box(:,:,5,:) = mean(temp_msr_vec_other(:, :, 6, :), 3);
%             mean_msr_box(:,:,6,:) = mean(temp_msr_vec_non(:, :, 6, :), 3);
%             ROI = {'HIPPOL main', 'HIPPOL one vs others', 'HIPPOL nonlin', ...
%                    'HIPPOR main', 'HIPPOR one vs others', 'HIPPOR nonlin'}; 
            mean_msr_box(:,:,1,:) = mean(temp_msr_vec_main(:, :, 3, :), 3);
            mean_msr_box(:,:,2,:) = mean(temp_msr_vec_other(:, :, 3, :), 3);
            mean_msr_box(:,:,3,:) = mean(temp_msr_vec_non(:, :, 3, :), 3);
            ROI = {'HIPPOLus1 main', 'HIPPOLus1 one vs others', 'HIPPOLus1 nonlin'}; 
            ROI_modif = 1; color_id = 3;
        end
        
    case 4
        if bilateral_integ
            mean_msr_box(:,:,1,:) = mean(temp_msr_vec_main(:, :, roii4, :), 3);
            mean_msr_box(:,:,2,:) = mean(temp_msr_vec_other(:, :, roii4, :), 3);
            mean_msr_box(:,:,3,:) = mean(temp_msr_vec_non(:, :, roii4, :), 3);
            ROI = {'VentStr main', 'VentStr one vs others', 'VentStr nonlin'};
            ROI_modif = 1; color_id = 4;
        else
%             mean_msr_box(:,:,1,:) = mean(temp_msr_vec_main(:, :, 7, :), 3);
%             mean_msr_box(:,:,2,:) = mean(temp_msr_vec_other(:, :, 7, :), 3);
%             mean_msr_box(:,:,3,:) = mean(temp_msr_vec_non(:, :, 7, :), 3);
%             mean_msr_box(:,:,4,:) = mean(temp_msr_vec_main(:, :, 8, :), 3);
%             mean_msr_box(:,:,5,:) = mean(temp_msr_vec_other(:, :, 8, :), 3);
%             mean_msr_box(:,:,6,:) = mean(temp_msr_vec_non(:, :, 8, :), 3);
%             ROI = {'VentStrL main', 'VentStrL one vs others', 'VentStrL nonlin', ...
%                 'VentStrR main', 'VentStrR one vs others', 'VentStrR nonlin'};
            mean_msr_box(:,:,1,:) = mean(temp_msr_vec_main(:, :, 4, :), 3);
            mean_msr_box(:,:,2,:) = mean(temp_msr_vec_other(:, :, 4, :), 3);
            mean_msr_box(:,:,3,:) = mean(temp_msr_vec_non(:, :, 4, :), 3);
            ROI = {'VentStrL main', 'VentStrL one vs others', 'VentStrL nonlin'};
            ROI_modif = 1; color_id = 4;
        end
        
end

% ROI = {'LPFC'};
% ROI = {'LPFC', 'V1'};
% ROI = {'LPFC', 'V1', 'VentStr'};

% ROI = {'PALL', 'Put', 'Cd', 'VentStr'};

disp('=================== mean_msr_box manipulation done! ===================')


%%

% for vEVENTi_global = 1:8 % for vs BHV per trial event


%% visualization setting

% fdr = 0; % no multiple testing correction
fdr = 'Benjamini-Hochberg';


LoadMyColor;
LoadDefaultSettings;
if bilateral_integ && length(ROI) == 8
    region_color_key = [4, ...  % ilPFC
        9, ...  % V1
        16,...  % HIPPO
        19];    % VentStr
else
    region_color_key = [];
end
close all
% axis square; box off;
% set(0,'DefaultAxesFontSize',16,'DefaultTextFontSize',16);
set(0,'DefaultAxesFontSize',14,'DefaultTextFontSize',14);
% set(0,'DefaultAxesFontSize',10,'DefaultTextFontSize',10);
% set(0,'DefaultAxesFontSize',10,'DefaultTextFontSize',8);

% ===== visualization mode =====
% vis_mode = 'regional_info';
% vis_mode = 'temporal_info'; temp_mode = 'group';     set(0,'DefaultAxesFontSize',12,'DefaultTextFontSize',12);
% vis_mode = 'temporal_info'; temp_mode = 'individual';     set(0,'DefaultAxesFontSize',12,'DefaultTextFontSize',12);
vis_mode = 'corr_behav';  corr_mode = 'original';  set(0,'DefaultAxesFontSize',15,'DefaultTextFontSize',15);
% vis_mode = 'corr_behav';  corr_mode = 'ROI effects';  set(0,'DefaultAxesFontSize',10,'DefaultTextFontSize',10);
% vis_mode = 'corr_behav';  corr_mode = 'dichotomy effects';  set(0,'DefaultAxesFontSize',10,'DefaultTextFontSize',10);
% vis_mode = 'exhaust_classif';       set(0,'DefaultAxesFontSize',14,'DefaultTextFontSize',14);
% vis_mode = 'corr_pca'; pca_seed = 2021; set(0,'DefaultAxesFontSize',10,'DefaultTextFontSize',10);

% --- PCA setting ---
% LABEL_pca = 'Goal(6,7,8)'; restrict_pca = 'G'; vis_class_pca = 'all'; method_pca = {'PR', 'ER', 'nPC'}; strict_pca = 0; norm_pca = 'basic'; shuffle_pca = []; 
% LABEL_pca = 'Goal(6,7,8)'; restrict_pca = 'G'; vis_class_pca = 'all'; method_pca = {'MV'}; strict_pca = 0; norm_pca = 'none'; shuffle_pca = []; 
% LABEL_pca = 'Goal(6,7,8)'; restrict_pca = 'G'; vis_class_pca = 'all'; method_pca = {'pat_size', 'N_eigs'}; strict_pca = 0; norm_pca = 'basic'; shuffle_pca = []; 

% LABEL_pca = 'GoalCond'; restrict_pca = []; vis_class_pca = 1; method_pca = {'PR', 'ER', 'nPC'}; strict_pca = 0; norm_pca = 'basic'; shuffle_pca = []; % GoalCond:1 -- specific goal block
% LABEL_pca = 'GoalCond'; restrict_pca = []; vis_class_pca = 1; method_pca = {'PR', 'ER', 'nPC'}; strict_pca = 0; norm_pca = 'undersample100'; shuffle_pca = []; % GoalCond:1 -- specific goal block
% LABEL_pca = 'GoalCond'; restrict_pca = []; vis_class_pca = 1; method_pca = {'MV'}; strict_pca = 0; norm_pca = 'none'; shuffle_pca = []; % GoalCond:1 -- specific goal block
% LABEL_pca = 'GoalCond'; restrict_pca = []; vis_class_pca = 1; method_pca = {'pat_size', 'N_eigs'}; strict_pca = 0; norm_pca = 'basic'; shuffle_pca = []; % GoalCond:1 -- specific goal block
% LABEL_pca = 'GoalCond'; restrict_pca = []; vis_class_pca = 1; method_pca = {'pat_size', 'N_eigs'}; strict_pca = 0; norm_pca = 'none'; shuffle_pca = []; % GoalCond:1 -- specific goal block

% LABEL_pca = 'GoalCond'; restrict_pca = []; vis_class_pca = 2; method_pca = {'PR', 'ER', 'nPC'}; strict_pca = 0; norm_pca = 'basic'; shuffle_pca = []; % GoalCond:2 -- flexible goal block
% LABEL_pca = 'GoalCond'; restrict_pca = []; vis_class_pca = 2; method_pca = {'pat_size', 'N_eigs'}; strict_pca = 0; norm_pca = 'basic'; shuffle_pca = []; % GoalCond:2 -- flexible goal block

region_sort = 1;
% region_sort = 0;
BHV_load = 1;
% BHV_load = 0;
% scatter_color = 1;
scatter_color = 0;
disp_clsf_N = inf;
% disp_clsf_N = 30;
fprintf('----------- %s ----------- \n', vis_mode)

% ===== measure name =====
if contains(measure, 'diff')
    % 2021-11-17 temp
%     msr_name = ['SD - CCGP across ' CClabel_name];
    msr_name = ['CCGP across ' CClabel_name ' - SD'];
%     msr_name = [measure_cell{1} ' - ' measure_cell{2}];
%     msr_name = [measure_cell{1} ' - ' measure_cell{2} ' (trained:', trained_mdl_info, ')'];
elseif contains(measure, 'accs')
    if contains(save_name, 'fitrlinear')
        msr_name = 'CV correlation';
    else
        msr_name = 'Mean accuracy';
    end
elseif contains(measure, 'acc_mean')
    if isempty(CClabel_name)
        if ~isempty(SepClass) && ~iscell(SepClass)
            msr_name = [SepClass ' separability'];
        else
            msr_name = 'Separability dimension';
        end
    else
        msr_name = ['CCGP across ' CClabel_name];
    end
elseif contains(measure, 'separable')
    msr_name = 'Sensitivity dimension';
elseif contains(measure, 'separable2')
    msr_name = 'Sensitivity dimension';
elseif contains(measure, 'CCGP')
    msr_name = ['CCGP (trained:', trained_mdl_info, ')'];
elseif contains(measure, 'cosine')
    msr_name = 'Cosine similarity';
    if length(LABEL) == 2
        LABEL = {[LABEL{1} ' vs ' LABEL{2}]};
    end
elseif contains(measure, 'corr_sig')
    msr_name = 'Correlation significance';
    if length(LABEL) == 2
        LABEL = {[LABEL{1} ' vs ' LABEL{2}]};    
    end
elseif contains(measure, 'GI')
    msr_name = {'Generalization index', [' across ' CClabel_name]};
elseif contains(measure, 'ratio')
%     msr_name = [measure_cell{1}, '/', measure_cell{2}];
    msr_name = ['Ratio(CCGP across ' CClabel_name ', SD)'];
else
    msr_name = measure;
end

% ===== visualized id index =====
vIDi = 1:N_sbj;

% ===== visualized label index =====
vLABELi = 1:length(LABEL);
% vLABELi = 1;
% vLABELi = 1:5;
% vLABELi = [1 2 3 5];
% vLABELi = 6:10;
% vLABELi = [6 7 8 10];
% vLABELi = 11:15;
% vLABELi = [11 12 13 15];
% vLABELi = 16:20;
% vLABELi = [16:18 20];
% vLABELi = 1:5:16;
% vLABELi = 2:5:17;
% vLABELi = 3:5:18;
% vLABELi = 4:5:19;
% vLABELi = 5:5:20;
% vLABELi = 3;
% vLABELi = 4:7;
% vLABELi = [1 2];
% vLABELi = [1 3];
% vLABELi = [1 4];
% vLABELi = [2 3];
% vLABELi = [3 4];
% vLABELi = [3 5];
% vLABELi = [2 5];
% vLABELi = [4 5];
% vLABELi = [5 6];
% vLABELi = [6 7];
% vLABELi = [7 8];

% ===== visualized region index =====
vROIi = 1:length(ROI);
% vROIi = 1:4;
% vROIi = 1;
% vROIi = [1 2];
% vROIi = [1 3];
% vROIi = [1 4];
% vROIi = [1 2 4];
% vROIi = [2 3];
% vROIi = [2 4];
% vROIi = 1:5;

% vROIi = [1 2 5 6];
% vROIi = [3 4 7 8];
% vROIi = [1 5];
% vROIi = [3 7];
% vROIi = [2 6];
% vROIi = [4 8];


% ===== visualized event index =====
if strcmp(vis_mode, 'temporal_info')
    vEVENTi = 1:length(EVENT);
else
    
%     vEVENTi = 1:length(EVENT);
    
%     --- S3 (1 2 4 5 7 8 9) ---
    % vEVENTi = vEVENTi_global;
%     vEVENTi = 1;
%     vEVENTi = 2;
%     vEVENTi = 3;
%     vEVENTi = 4;
%     vEVENTi = 5;
%     vEVENTi = 6;
%     vEVENTi = 7;
%     vEVENTi = 8;
%     vEVENTi = [2 3];
%     vEVENTi = [4 5];
%     vEVENTi = [6 7];
%     vEVENTi = 6;
%     vEVENTi = 2:6;
%     vEVENTi = 2:7;
%     vEVENTi = 2:5;
%     vEVENTi = 1:6;
%     vEVENTi = 1:7;

%     --- S3 (1 2 4 5 7 8) ---
%     vEVENTi = 2:6;
    vEVENTi = 2:5;
%     vEVENTi = [2 3];
%     vEVENTi = [4 5];
%     vEVENTi = 6;
    
end

if iscell(EVENT)
    ev_title = num2str([EVENT{vEVENTi}]);                           % title
    ev_names = cellfun(@num2str, EVENT, 'UniformOutput', false);    % ticks label
else
%     ev_title = [EVENT_NAME{EVENT(vEVENTi)}];
%     ev_names = EVENT_NAME(EVENT(vEVENTi));
    ev_title = [EVENT_NAME{vEVENTi}];
    ev_names = EVENT_NAME(vEVENTi);     
end

% ===== target task variables =====
% BEHAV = {'R', 'ChoOpt2'};    % cosyne21 fig. C
% ylims =  {[12 22], [0.5 0.8]};

% bhv_cond = []; bhv_sfx = [];   % all trials
% bhv_cond = 1; bhv_sfx = ' in G';  % only G trials
% bhv_cond = 2; bhv_sfx = ' in H';  % only H trials
% bhv_cond = 3; bhv_sfx = 'in uL';  % only uL trials
% bhv_cond = 4; bhv_sfx = 'in uH';  % only uH trials
bhv_cond = 'G'; bhv_sfx = ' (G)';   % G trials
% bhv_cond = 'H'; bhv_sfx = ' (H)';   % G trials
% bhv_cond = 'G uL'; bhv_sfx = ' (G uL)';   % G trials, low uncertainty
% bhv_cond = 'G uH'; bhv_sfx = ' (G uH)';   % G trials, high uncertainty
% bhv_cond = {'G uL', 'G uH'}; bhv_sfx = ' (G uL-uH)';
% bhv_cond = 'uL'; bhv_sfx = ' (uL)';   % low uncertainty
% bhv_cond = 'uH'; bhv_sfx = ' (uH)';   % high uncertainty
% bhv_cond = {'uL', 'uH'}; bhv_sfx = ' (uL-uH)';   
% bhv_cond = 'G uEarly'; bhv_sfx = ' (G uEarly)';   % G trials, UC early phase
% bhv_cond = 'G uLate'; bhv_sfx = ' (G uLate)';   % G trials, UC late phase
% bhv_cond = {'G uLate', 'G uEarly'}; bhv_sfx = ' (G uL-uE)';   % G trials, uLate - uEarly
% bhv_cond = 'binPMB0'; bhv_sfx = ' (binPMB0)';   % binPMB0 trials
% bhv_cond = 'binPMB1'; bhv_sfx = ' (binPMB1)';   % binPMB1 trials
% bhv_cond = 'binPMB1 uEarly'; bhv_sfx = ' (binPMB1 uEarly)';   % binPMB1 trials, UC early phase
% bhv_cond = 'binPMB1 uLate'; bhv_sfx = ' (binPMB1 uLate)';   % binPMB1 trials, UC late phase
% bhv_cond = {'binPMB1 uLate', 'binPMB1 uEarly'}; bhv_sfx = ' (binPMB1 uL-uE)';   % binPMB1, uLate - uEarly
% bhv_cond = 'binPMB1 uL'; bhv_sfx = ' (binPMB1 uL)';   
% bhv_cond = 'binPMB1 uH'; bhv_sfx = ' (binPMB1 uH)';   
% bhv_cond = {'binPMB1 uL', 'binPMB1 uH'}; bhv_sfx = ' (binPMB1 uL-uH)';   
% bhv_cond = 'Goal(6,7,8)xUCSwitch'; bhv_sfx = ' (Switch)';   % goal x UC context switch
% bhv_cond = 'Goal(6,7,8)xUCStay'; bhv_sfx = ' (Stay)';   % goal x UC context stay
% bhv_cond = {'Goal(6,7,8)xUCSwitch', 'Goal(6,7,8)xUCStay'}; bhv_sfx = ' (Switch - Stay)';

% bhv_horz = 6;
% bhv_horz = 5;
% bhv_horz = 4;
% bhv_horz = 3;
% bhv_horz = 2;
bhv_horz = 1;

% % --- exhaustive BHV measures ---
switch Exp
    case {'Lee2014', 'Lee2014pp2'}
        
%         BEHAV = {'R'};
        
        % optimality, flexibility, stability (trial-event averaged)
        % BEHAV = {'ChoSwitch1n2(GoalSwitch)', 'ChoConsist1n2', 'ChoOptBinNaN1n2', ...
        %     'mb tau uL', 'mb tau uH', 'mb tau uH-uL', ...
        %     'mb lr uL', 'mb lr uH', 'mb lr uH-uL'};
        BEHAV = {'ChoSwitch1n2(GoalSwitch)', 'ChoConsist1n2', 'ChoOptBinNaN1n2'};
        % BEHAV = {'ChoSwitch1n2(GoalSwitch)', 'ChoConsist1n2', 'ChoOptBin1n2'};
        
% %         For BHV vs goal SD
%         BEHAV = {'ChoOptBin1', 'ChoOptBin2', 'ChoOptBin1n2', 'R' ...
%             'ChoiceSwitch1(GoalSwitch)', 'ChoiceSwitch2(GoalSwitch)', 'ChoiceSwitch1n2(GoalSwitch)', 'PMB', ...
%             'ChoiceSwitch1(GoalStay)', 'ChoiceSwitch2(GoalStay)', 'ChoiceSwitch1n2(GoalStay)', 'lr', ...
%             'ChoConsist1', 'ChoConsist2', 'ChoConsist1n2', 'invTemp'};
        
        % For parameters vs UC SD
%         BEHAV = {'ThrSPE', 'lrRPE', 'MB2MF', 'MF2MB', 'invTemp', 'lr', ...
%             'SPE2', 'SPE3', 'SPE', 'RPE2', 'RPE3', 'RPE', ... 
%             'relMB2', 'relMB3', 'relMB', 'relMF2', 'relMF3', 'relMF', ...
%             'relMAX2', 'relMAX3', 'relMAX', 'relDiff2', 'relDiff3', 'relDiff'};

        % UC choice efx vs UC SD
%         BEHAV = {'UncCond--A1(anovan)', 'UncCond--A2 state2(anovan)', ... 
%             'UncCond--A2 state3(anovan)', 'UncCond--A2 state4(anovan)', 'UncCond--A2 state5(anovan)'};
        
%         BEHAV = {'ChoOpt1', 'ChoOpt2', 'ChoOpt1n2'};
        
%         BEHAV = {'ChoOpt1', 'ChoOpt2', 'ChoOpt1n2', 'PMB'};
%         BEHAV = {'ChoOpt1', 'ChoOpt2', 'ChoOpt1n2', 'PMB' ...
%             'ChoOptVar1', 'ChoOptVar2', 'ChoOptVar', 'lr', ...
%             'ChoConsist1', 'ChoConsist2', 'ChoConsist1n2', 'StrictChoCon', ...
%             'ChoiceSwitch1', 'ChoiceSwitch2', 'ChoiceSwitch1n2', 'invTemp'};

%         BEHAV = {'ChoOpt1', 'ChoOpt2', 'ChoOpt1n2', 'PMB' ...
%             'ChoOptVar1', 'ChoOptVar2', 'ChoOptVar', 'lr', ...
%             'ChoConsist1', 'ChoConsist2', 'ChoConsist1n2', 'StrictChoCon', ...
%             'ChoCon1 ignoring context', 'ChoCon2 ignoring context', 'ChoCon ignoring context', 'invTemp'}; % for Lee 2014

%         BEHAV = {'ChoOpt1', 'ChoOpt2', 'ChoOpt1n2', 'R' ...
%             'ChoOptVar1', 'ChoOptVar2', 'ChoOptVar', 'PMB', ...
%             'ChoConsist1', 'ChoConsist2', 'ChoConsist1n2', 'StrictChoCon', ...
%             'ChoCon1 ignoring context', 'ChoCon2 ignoring context', 'ChoCon ignoring context'}; % for Lee 2014
        % BEHAV = {'ChoOpt1', 'ChoOpt2', 'ChoOpt1n2', 'R' ...
%     'ChoOptVar1', 'ChoOptVar2', 'ChoOptVar', 'PMB', ...
%     'ChoConsist1', 'ChoConsist2', 'ChoConsist1n2', 'StrictChoCon'}; % for Lee 2014
% BEHAV = {'ChoOpt1', 'ChoOpt2', 'ChoOpt1n2', 'R' ...
%     'ChoConsist1', 'ChoConsist2', 'ChoConsist1n2', 'PMB', ...
%     'ChoCon1 ignoring context', 'ChoCon2 ignoring context', 'ChoCon ignoring context', ...
%     'StrictChoCon'}; % for Lee 2014
% BEHAV = {'ChoOpt1', 'ChoOpt2', 'ChoOpt1n2', 'R' ...
%     'ChoConsist1', 'ChoConsist2', 'ChoConsist1n2', 'PMB', ...
%     'ChoCon1 ignoring context', 'ChoCon2 ignoring context', 'ChoCon ignoring context', ...
%     'StrictChoCon'}; % for Lee 2014
    case 'Kim2019'
%         NLLname = 'OptimalA2SoftmaxNLLUncCond1';
        NLLname = 'OptimalA2SoftmaxNLL';
        BEHAV = {'ChoOpt NLL', 'MaxGoal NLL', 'AchGoal NLL', ... 
            'ChoOpt - MaxGoal NLL', 'ChoOpt - AchGoal NLL', 'MaxGoal - AchGoal NLL'};
%         BEHAV = {'ChoOpt1', 'ChoOpt2', 'ChoOpt1n2', 'PMB' ...
%             'ChoOptVar1', 'ChoOptVar2', 'ChoOptVar', 'PMBvar', ...
%             'ChoCon1 ignoring context', 'ChoCon2 ignoring context', 'ChoCon ignoring context'}; % for Kim 2019
        % BEHAV = {'ChoOpt1', 'ChoOpt2', 'ChoOpt1n2', 'R' ...
%     'ChoCon1 ignoring context', 'ChoCon2 ignoring context', 'ChoCon ignoring context', ...
%     'PMB'}; % for Kim 2019
end

% BEHAV = {'ChoOpt1', 'ChoOpt2', 'ChoOpt1n2', 'Hit' ...
%     'ChoConsist1', 'ChoConsist2', 'ChoConsist1n2', 'R', ...
%     'ChoCon1 ignoring context', 'ChoCon2 ignoring context', 'ChoCon ignoring context', ...
%     'StrictChoCon'};
% BEHAV = {'Hit', 'ChoOpt1', 'ChoCon1', 'ChoCon1 per goal', 'RTA1'...
%     'R', 'ChoOpt2', 'ChoCon2', 'ChoCon2 per goal', 'RTA2'};
% BEHAV = {'relMB', 'relMF', 'relMAX', 'relDiff'};
% BEHAV = {'Hit', 'ChoOpt1', 'ChoConsist1', 'RTA1'...
%     'R', 'ChoOpt2', 'ChoConsist2', 'RTA2'};
% BEHAV = {'Hit', 'ChoOpt1', 'ChoCon1', 'RTA1'...
%     'R', 'ChoOpt2', 'ChoCon2', 'RTA2'};
% BEHAV = {'Hit', 'ChoOpt1', 'RTA1'...
%     'R', 'ChoOpt2', 'RTA2'};
% BEHAV = {'CESD'};
% BEHAV = {'lr'};
ylims = {};

% BEHAV = {'Hit', 'R', 'RTA1', 'RTA2', 'ChoOpt1n2', ... 
%     'SPE3', 'RPE3', 'Qarb2', 'relMB3', 'relMF3', 'relMAX3'};
% BEHAV = {'Hit', 'R', 'RTA1', 'RTA2', 'ChoOpt1n2', 'ChoOptVar'};
% BEHAV = {'Hit', 'R', 'RTA1', 'RTA2', 'ChoOpt1n2', 'ChoOptVar', ...
%     'relMAX'};
% BEHAV = {'Hit', 'R', 'RTA1', 'RTA2', 'ChoOpt1', 'ChoOptVar'};
% BEHAV = {'Hit', 'R', 'RTA1', 'RTA2', 'ChoOpt2', 'ChoOptVar2'};

sig_lvl = 0.05;

% ===== individual result subplot parameter =====
vert_slot = 4;  % regional_info, exhaust_classif
horz_slot = 6;  % temporal_info

% ===== ylim =====
if strcmp(measure, 'separable') || strcmp(measure, 'separable2')
    
    y_lower = 0;
%     y_lower = 0.2;
%     y_lower = 0.4;
%     y_lower = [];
%     y_upper = [];
    % y_upper = 0.5;
%     y_upper = 0.6;
    y_upper = 1;
    
    x_lower = 0.2;
%     x_lower = 0.4;
%     x_lower = 0.5;
%     x_lower = 0;
%     x_upper = 0.6;
    x_upper = 1;
    
elseif strcmp(measure, 'accs') || strcmp(measure, 'acc_mean')
    
%     y_lower = 0;
%     y_lower = 0.25;
%     y_lower = 0.4;
%     y_lower = 0.45;
%     y_lower = 0.5;
%     y_lower = chance_lvl;
%     y_lower = chance_lvl - 0.1;
%     y_lower = chance_lvl - 0.05;
    y_lower = chance_lvl - 0.025;
%     y_lower = 0.55;
%     y_lower = 0.6;
%     y_lower = [];
%     y_upper = [];
%     y_upper = 0.45;
%     y_upper = 0.5;
%     y_upper = 0.55;
%     y_upper = 0.56;
%     y_upper = 0.58;
%     y_upper = 0.6;
%     y_upper = 0.62;
%     y_upper = 0.65;
%     y_upper = chance_lvl + 0.05;
%     y_upper = chance_lvl + 0.06;
%     y_upper = chance_lvl + 0.075;
    y_upper = chance_lvl + 0.08;
%     y_upper = chance_lvl + 0.1;
%     y_upper = chance_lvl + 0.15;
%     y_upper = chance_lvl + 0.2;
%     y_upper = chance_lvl + 0.25;
%     y_upper = chance_lvl + 0.4;
%     y_upper = 0.7;
%     y_upper = 0.75;
%     y_upper = 0.8;

%     x_lower = chance_lvl - 0.1;
%     x_upper = chance_lvl + 0.2;
    
    x_lower = 0.45;
%     x_lower = 0.5;
%     x_upper = 0.55;
    x_upper = 0.6;
%     x_upper = 0.65;
%     x_upper = 0.67;
%     x_upper = 0.7;
%     x_upper = 0.75;

elseif strcmp(measure, 't-stat')
    y_lower = 1;
%     y_lower = [];
%     y_upper = [];
    y_upper = 25;
    
    x_lower = 0;
%     x_lower = [];
%     x_upper = [];
    x_upper = 50;
    
elseif strcmp(measure, 'PR_mean') 
    
%     y_lower = 5;
    y_lower = 10;
%     y_upper = 25;
    y_upper = 30;
    
%     x_lower = 5;
    x_lower = 10;
    x_upper = 30;
    
elseif strcmp(measure, 'nPC_mean')
    
    y_lower = 0;
%     y_upper = 100;
    y_upper = 120;
    x_lower = 0;
%     x_upper = 100;
    x_upper = 120;
    
elseif strcmp(measure, 'CCGP')
    
%     y_lower = chance_lvl * 0;
%     y_lower = 0.25;
%     y_lower = chance_lvl;
%     y_lower = chance_lvl - 0.2;
    y_lower = chance_lvl - 0.1;
%     y_lower = 0.45;
%     y_upper = 0.45;
%     y_upper = 0.5;
%     y_upper = chance_lvl + 0.2;
    y_upper = chance_lvl + 0.15;
%     y_upper = chance_lvl + 0.5;
%     y_upper = 0.65;
%     y_upper = 1;
    
    x_lower = 0.2;
    x_upper = 1;
    
elseif contains(measure, 'diff')
%     y_lower = chance_lvl - 0.05;
    y_lower = chance_lvl - 0.025;
    y_upper = chance_lvl + 0.025;
%     y_upper = chance_lvl + 0.2;
    x_lower = [];
    x_upper = [];
    
elseif contains(measure, 'cosine')
    y_lower = -.1;
    y_upper = .1;
    x_lower = [];
    x_upper = [];
    
elseif contains(measure, 'ratio')
    y_lower = 0.95;
    y_upper = 1.03;
    x_lower = [];
    x_upper = [];
       
else 
    
    y_lower = [];
    y_upper = [];
    x_lower = [];
    x_upper = [];
    
end

if strcmp(vis_mode, 'temporal_info')
    if strcmp(temp_mode, 'individual')
        y_lower = [];
        y_upper = [];
        x_lower = [];
        x_upper = [];
    end
end

% ===== significance level =====
sgnf_lvl = 0.05;


disp('step3 done')

%% visualization

close all

switch vis_mode
    
    % ---------------------------------------------------------------------
    case 'corr_behav'
        
        % load mean behavioral variable for each sbj
        if BHV_load
            fprintf('Behavior variable loading: ')
            BHVs = cell(1, length(BEHAV));
            
            if iscell(bhv_cond)
                bhv_conds_cell = bhv_cond;
            else
                bhv_conds_cell = {bhv_cond};
            end
            
            for bhvi = 1:length(BEHAV)
                
                bhvs_mat = [];
                for bhvcondi = 1:length(bhv_conds_cell)
                    bhv_cond = bhv_conds_cell{bhvcondi};

                    % bhv loop: load
                    fprintf('%d/%d ', bhvi, length(BEHAV))

                    bhvs = arbMBMF_load_var(Exp, BEHAV{bhvi}, 0, bhv_cond);

                    bhvs_mat = [bhvs_mat; bhvs];
                    
                end % bhvcondi
                
                % store
                if size(bhvs_mat, 1) == 2
                    bhvs_vec = bhvs_mat(1,:) - bhvs_mat(2,:);
                else
                    bhvs_vec = bhvs_mat;
                end
                BHVs{bhvi} = bhvs_vec;
                
            end % bhvi
            
            disp(' ')
        end
        
        % 201110
        rois = ROI(vROIi);
        labels = LABEL(vLABELi);
        msrvecs = cell(length(vROIi), length(vLABELi));
        msrmat = mean(mean_msr_box(vIDi, vLABELi, vROIi, vEVENTi), 4);
        
        switch corr_mode
            case 'original' % corr. btw two variables: corr(BHV, neural)

                r_stats = nan(length(BEHAV), length(vLABELi), length(vROIi));
                t_stats = nan(length(BEHAV), length(vLABELi), length(vROIi));
                t_dfs = nan(length(BEHAV), length(vLABELi), length(vROIi));
                p_vals = nan(length(BEHAV), length(vLABELi), length(vROIi));
                sig_texts = cell(length(BEHAV), length(vLABELi), length(vROIi));
                
                % scatter for each ROI & LABEL
                for vroii = 1:length(vROIi)
                    for vlabi = 1:length(vLABELi)
                        figure('Name', [ROI{vROIi(vroii)} ', Event: ' ev_title], ... 
                            'Position', [(vroii-1)/length(vROIi)*1920 1080/2-80 1920/12 1080/2]) % [left bottom width height] 
                        
                        % averaged decoding performance
                        msrvec_id = squeeze(nanmean( ...
                            mean_msr_box(vIDi, vLABELi(vlabi), vROIi(vroii), vEVENTi), ...
                            4) );
                        msrvecs{vroii, vlabi} = msrvec_id;  % 201110
                        
                        for bhvi = 1:length(BEHAV)
                            
                            bhvs = BHVs{bhvi};
                            
                            % scatter plot
                            subplot(ceil(length(BEHAV)/bhv_horz), bhv_horz, bhvi)
                            if scatter_color
                                scatter(msrvec_id, bhvs(vIDi), ...
                                    20, ...                             % size
                                    linspace(1,10,length(vIDi)), ...    % colored points
                                    'filled')
                                p_val_form = '%2.2f';
                            else
                                if bilateral_integ
                                    if ~isempty(region_color_key)
                                        % COSYNE 2023 poster fig 4
                                        scatter(msrvec_id, bhvs(vIDi), ...
                                            11, ...
                                            [.5 .5 .5], ...
                                            'filled')
%                                         scatter(msrvec_id, bhvs(vIDi), ...
%                                             30, ...
%                                             colors{region_color_key(vroii)}, ...
%                                             'filled')
                                    else
                                        scatter(msrvec_id, bhvs(vIDi), ...
                                            11, ...
                                            'filled')
                                    end
                                    %                         axis square
                                    p_val_form = '%2.4f';
                                    %                         xlim([y_lower y_upper]);
                                    if ~isempty(ylims); ylim(ylims{bhvi}); end
                                else
                                    if ~isempty(region_color_key)
                                        scatter(msrvec_id, bhvs(vIDi), ...
                                            11, ...
                                            colors{region_color_key(ceil(vroii/2))}, ...
                                            'filled')
                                    else
                                        scatter(msrvec_id, bhvs(vIDi), ...
                                            11, ...
                                            'filled')
                                    end
                                    %                         axis square
                                    p_val_form = '%2.4f';
                                    %                         xlim([y_lower y_upper]);
                                    if ~isempty(ylims); ylim(ylims{bhvi}); end
                                end
                            end
                            hold on;
                            
                            % significant correlation check
                            [R, p] = corrcoef(bhvs(vIDi), msrvec_id);

                            % Extract the correlation coefficient
                            r = R(1,2);  % off-diagonal element
                            n = length(msrvec_id);  % or y, assuming x and y are the same length
                            % Compute t-statistic
                            tstat = r * sqrt((n - 2) / (1 - r^2));
                            % Degrees of freedom
                            df = n - 2;
                            stats.tstat = tstat;
                            stats.df = df;
                            stats.p = 2 * (1 - tcdf(abs(tstat), df));

                            X = msrvec_id; Y = bhvs(vIDi)';
                            ttl_sfx = [];
                            if p(2,1) < sig_lvl
                                ttl_sfx = [' r=' num2str(R(2,1),'%2.2f') ', p=' num2str(p(2,1), p_val_form)];
                                %                         X = bhvs(vIDi)'; Y = msrvec_id;
                                
                                % COSYNE 2023 poster fig 4
                                tbl = table(X, Y); mdl = fitlm(tbl,'Y~X');
                                crr_fit = plotAdded(mdl,'X', 'Marker', 'none');
%                                 crr_fit(2).Color = 'k';
%                                 crr_fit(3).Color = 'k';
                                legend('off')
                                %                         disp([ROI{vROIi(vroii)} ' ', ...
                                %                             LABEL{vLABELi(vlabi)} ' ', ...
                                %                             BEHAV{bhvi}, ...
                                %                             ' r=' sprintf('%.2f', R(2,1)) ', p=' sprintf('%.2e', p(2,1))])
                                % disp([ev_title ' ', ...
                                %     ROI{vROIi(vroii)} ' ', ...
                                %     LABEL{vLABELi(vlabi)} ' ', ...
                                %     BEHAV{bhvi}, bhv_sfx, ...
                                %     ' r=' sprintf('%.3f', R(2,1)) ', p=' sprintf('%.2e', p(2,1))])
                            end
                            hold off

                            % To print statistic and p value
                            p = p(1,2);
                            if p < 0.0001
                                sig_text = '****';
                            elseif p < 0.001
                                sig_text = '***';
                            elseif p < 0.01
                                sig_text = '**';
                            elseif p < 0.05
                                sig_text = '*';
                            else
                                sig_text = 'n.s.';
                            end


                            t_stats(bhvi, vlabi, vroii) = tstat;
                            t_dfs(bhvi, vlabi, vroii) = df;
                            r_stats(bhvi, vlabi, vroii) = r;
                            p_vals(bhvi, vlabi, vroii) = p;
                            sig_texts{bhvi, vlabi, vroii} = sig_text;
                            
                            % COSYNE 2023 poster fig 3
                            title(ttl_sfx)
%                             title({[ROI{vROIi(vroii)} ' ' LABEL{vLABELi(vlabi)}], ...
%                                 ['Event: ' ev_title], ttl_sfx})
                            %                     title({[ROI{vROIi(vroii)} ' ' LABEL{vLABELi(vlabi)}, ...
                            %                         ' ' ttl_sfx], ['Event: ' ev_title]})
                            %                     title([ROI{vROIi(vroii)} ' ' LABEL{vLABELi(vlabi)} ' ' ttl_sfx])
                            %                     title([ROI{vROIi(vroii)} ' ' LABEL{vLABELi(vlabi)}])
                            
                            %                     xlabel(['Average ' BEHAV{bhvi} ' ' bhv_sfx])
                            %                     ylabel(['Mean ' measure], 'Interpreter', 'none')
                            
                            % COSYNE 2023 poster fig 3
                            ylabel([])
                            xlabel([])
%                             ylabel(['Average ' BEHAV{bhvi} ' ' bhv_sfx])
%                             xlabel(msr_name)

                            %                     xlabel(['Mean ' measure], 'Interpreter', 'none')
                            
                            LABEL_temp = LABEL{vLABELi(vlabi)};
                            ROI_temp = ROI{vROIi(vroii)};
                            BEHAV_temp = BEHAV{bhvi};
                            EVENT_temp = EVENT(vEVENTi);
                            
                            % save results
                            if iscell(msr_name)
                                msr_name_save = [msr_name{:}];
                            else
                                msr_name_save = msr_name;
                            end
                            
                            try
                                save(['data/ccgp_result/corr_behav_original/' erase(num2str(EVENT_temp), ' ') '/' Exp '_' msr_name_save '_' LABEL_temp '_' ROI_temp '_' BEHAV_temp bhv_sfx], ...
                                    'R', 'p', 'X', 'Y', 'BEHAV_temp', 'bhv_conds_cell', 'Exp', 'msr_name', 'vIDi', 'LABEL_temp', 'ROI_temp', 'EVENT', 'vEVENTi')
                            catch
                                mkdir(['data/ccgp_result/corr_behav_original/' erase(num2str(EVENT_temp), ' ')])
                                save(['data/ccgp_result/corr_behav_original/' erase(num2str(EVENT_temp), ' ') '/' Exp '_' msr_name_save '_' LABEL_temp '_' ROI_temp '_' BEHAV_temp bhv_sfx], ...
                                    'R', 'p', 'X', 'Y', 'BEHAV_temp', 'bhv_conds_cell', 'Exp', 'msr_name', 'vIDi', 'LABEL_temp', 'ROI_temp', 'EVENT', 'vEVENTi')
                            end
                                        
                        end
                    end
                end
                



                % #########################################################
                % FDR
                % #########################################################

                if fdr
                    p_vals_ori = p_vals;
                    for bhvi = 1:length(BEHAV)
                        for vlabi = 1:length(vLABELi)

                            pVals = squeeze(p_vals_ori(bhvi, vlabi, :)); % (1, nROI)
                            switch fdr
                                case 'Benjamini-Hochberg'
                                    [cpVals, cAlpha] = fdr_BH(pVals, .5);
                                    p_vals(bhvi, vlabi, :) = cpVals;
                                case 'Bonferroni'
                                    [cpVals, cAlpha] = fwer_bonf(pVals, .5);
                                    p_vals(bhvi, vlabi, :) = cpVals;

                            end
                        end
                    end

                    for vroii = 1:length(vROIi)
                        for vlabi = 1:length(vLABELi)
                            figure('Name', [ROI{vROIi(vroii)} ', Event: ' ev_title], ...
                                'Position', [(vroii-1)/length(vROIi)*1920 1080/2-80 1920/12 1080/2]) % [left bottom width height]

                            % averaged decoding performance
                            msrvec_id = squeeze(nanmean( ...
                                mean_msr_box(vIDi, vLABELi(vlabi), vROIi(vroii), vEVENTi), ...
                                4) );
                            msrvecs{vroii, vlabi} = msrvec_id;  % 201110

                            for bhvi = 1:length(BEHAV)

                                bhvs = BHVs{bhvi};

                                % scatter plot
                                subplot(ceil(length(BEHAV)/bhv_horz), bhv_horz, bhvi)
                                if scatter_color
                                    scatter(msrvec_id, bhvs(vIDi), ...
                                        20, ...                             % size
                                        linspace(1,10,length(vIDi)), ...    % colored points
                                        'filled')
                                    p_val_form = '%2.2f';
                                else
                                    if bilateral_integ
                                        if ~isempty(region_color_key)
                                            % COSYNE 2023 poster fig 4
                                            scatter(msrvec_id, bhvs(vIDi), ...
                                                11, ...
                                                [.5 .5 .5], ...
                                                'filled')
                                            %                                         scatter(msrvec_id, bhvs(vIDi), ...
                                            %                                             30, ...
                                            %                                             colors{region_color_key(vroii)}, ...
                                            %                                             'filled')
                                        else
                                            scatter(msrvec_id, bhvs(vIDi), ...
                                                11, ...
                                                'filled')
                                        end
                                        %                         axis square
                                        p_val_form = '%2.4f';
                                        %                         xlim([y_lower y_upper]);
                                        if ~isempty(ylims); ylim(ylims{bhvi}); end
                                    else
                                        if ~isempty(region_color_key)
                                            scatter(msrvec_id, bhvs(vIDi), ...
                                                11, ...
                                                colors{region_color_key(ceil(vroii/2))}, ...
                                                'filled')
                                        else
                                            scatter(msrvec_id, bhvs(vIDi), ...
                                                11, ...
                                                'filled')
                                        end
                                        %                         axis square
                                        p_val_form = '%2.4f';
                                        %                         xlim([y_lower y_upper]);
                                        if ~isempty(ylims); ylim(ylims{bhvi}); end
                                    end
                                end
                                hold on;

                                % significant correlation check
                                [R, ~] = corrcoef(bhvs(vIDi), msrvec_id);

                                % Extract the correlation coefficient
                                r = R(1,2);  % off-diagonal element
                                n = length(msrvec_id);  % or y, assuming x and y are the same length
                                % Compute t-statistic
                                tstat = r * sqrt((n - 2) / (1 - r^2));
                                % Degrees of freedom
                                df = n - 2;
                                stats.tstat = tstat;
                                stats.df = df;
                                stats.p = 2 * (1 - tcdf(abs(tstat), df));

                                X = msrvec_id; Y = bhvs(vIDi)';
                                p = p_vals(bhvi, vlabi, vroii);
                                ttl_sfx = [];
                                if p < sig_lvl
                                    ttl_sfx = [' r=' num2str(R(2,1),'%2.2f') ', p=' num2str(p, p_val_form)];
                                    %                         X = bhvs(vIDi)'; Y = msrvec_id;

                                    % COSYNE 2023 poster fig 4
                                    tbl = table(X, Y); mdl = fitlm(tbl,'Y~X');
                                    crr_fit = plotAdded(mdl,'X', 'Marker', 'none');
                                    crr_fit(2).Color = 'r';
                                    crr_fit(3).Color = 'r';
                                    %                                 crr_fit(2).Color = 'k';
                                    %                                 crr_fit(3).Color = 'k';
                                    legend('off')
                                    %                         disp([ROI{vROIi(vroii)} ' ', ...
                                    %                             LABEL{vLABELi(vlabi)} ' ', ...
                                    %                             BEHAV{bhvi}, ...
                                    %                             ' r=' sprintf('%.2f', R(2,1)) ', p=' sprintf('%.2e', p(2,1))])
                                    % disp([ev_title ' ', ...
                                    %     ROI{vROIi(vroii)} ' ', ...
                                    %     LABEL{vLABELi(vlabi)} ' ', ...
                                    %     BEHAV{bhvi}, bhv_sfx, ...
                                    %     ' r=' sprintf('%.3f', R(2,1)) ', p=' sprintf('%.2e', p(2,1))])
                                end
                                hold off

                                % To print statistic and p value
                                if p < 0.0001
                                    sig_text = '****';
                                elseif p < 0.001
                                    sig_text = '***';
                                elseif p < 0.01
                                    sig_text = '**';
                                elseif p < 0.05
                                    sig_text = '*';
                                else
                                    sig_text = 'n.s.';
                                end

                                sig_texts{bhvi, vlabi, vroii} = sig_text;

                                % COSYNE 2023 poster fig 3
                                title(ttl_sfx)
                                
                                ylabel([])
                                xlabel([])
                                
                            end
                        end
                    end                    
                end


                % if length(measure_cell) > 1 
                %     save(['msr_vec/' [measure_cell{:}] '_' Exp '_' [LABEL{:}] '_' num2str(vEVENTi) '_' [bhv_conds_cell{:}]], ...
                %         'msrvecs', 'Exp', 'vIDi', 'LABEL', 'vLABELi', 'ROI', 'vROIi', 'EVENT', 'vEVENTi', 'BEHAV', 'bhv_conds_cell', 'BHVs')
                % else
                %     save(['msr_vec/' msr_name_save '_' Exp '_' [LABEL{:}] '_' num2str(vEVENTi) '_' [bhv_conds_cell{:}]], ...
                %         'msrvecs', 'Exp', 'vIDi', 'LABEL', 'vLABELi', 'ROI', 'vROIi', 'EVENT', 'vEVENTi', 'BEHAV', 'bhv_conds_cell', 'BHVs')
                % end

                % To print statistic and p value
                for bhvi = 1:length(BHVs)

                    disp(BHVs{bhvi})

                    for vlabi = 1:length(vLABELi)

                        disp(LABEL{vLABELi(vlabi)})

                        for roii = 1:length(vROIi)
                            
                            fprintf([ROI{vROIi(roii)} ': r=%.3f, t(%d)=%.3f\t p=%.3e (' ...
                                sig_texts{bhvi, vlabi, roii} ')\n'], r_stats(bhvi, vlabi, roii), ...
                                t_dfs(bhvi, vlabi, roii), t_stats(bhvi, vlabi, roii), p_vals(bhvi, vlabi, roii))
                        end


                        disp(' ')

                    end

                    disp(' ')

                end
                
            case 'ROI effects' % BHV = f(ROI1, ROI2, ...)
                
                EVENT_temp = EVENT(vEVENTi);
                
                % for each LABEL
                for vlabi = 1:length(vLABELi)
                    LABEL_temp = LABEL{vLABELi(vlabi)};
                    figure('Name', ['Event: ' ev_title])
                    
                    % design matrix
                    X = squeeze(nanmean( ...
                        mean_msr_box(vIDi, vLABELi(vlabi), vROIi, vEVENTi), 4)); % [N_sbj x N_roi]
                    
                    for bhvi = 1:length(BEHAV)
                        bhvs = BHVs{bhvi};
                        BEHAV_temp = BEHAV{bhvi};
                        y = bhvs(vIDi)';
                        tbl = array2table([X, y], 'VariableNames', cat(2, ROI(vROIi), {'y'})); 
                        temp_format = join(ROI(vROIi), '+');
                        mdl_format = ['y~' temp_format{:}];
                        mdl = fitlm(tbl, mdl_format, 'RobustOpts', 'ols');
%                         mdl = fitrlinear(X, y, 'Learner', 'leastsquares');
%                         mdl = fitlm(tbl);
                        disp(mdl)
                        subplot(ceil(length(BEHAV)/bhv_horz), bhv_horz, bhvi)
                        plotAdded(mdl)
                        
                        % save results
                        try
                            save(['data/ccgp_result/corr_behav_ROIeffects/' erase(num2str(EVENT_temp), ' ') '/' Exp '_' msr_name '_' LABEL_temp '_' BEHAV_temp bhv_sfx], ...
                                'mdl', 'tbl', 'X', 'y', 'BEHAV_temp', 'bhv_conds_cell', 'Exp', 'msr_name', 'vIDi', 'LABEL_temp', 'ROI', 'vROIi', 'EVENT', 'vEVENTi')
                        catch
                            mkdir(['data/ccgp_result/corr_behav_ROIeffects/' erase(num2str(EVENT_temp), ' ')])
                            save(['data/ccgp_result/corr_behav_ROIeffects/' erase(num2str(EVENT_temp), ' ') '/' Exp '_' msr_name '_' LABEL_temp '_' BEHAV_temp bhv_sfx], ...
                                'mdl', 'tbl', 'X', 'y', 'BEHAV_temp', 'bhv_conds_cell', 'Exp', 'msr_name', 'vIDi', 'LABEL_temp', 'ROI', 'vROIi', 'EVENT', 'vEVENTi')
                        end
                        % writematrix: not available before MATLAB 2019
%                         writematrix(tbl, ['ccgp_result/corr_behav_ROIeffects/' Exp '_' msr_name '_' LABEL_temp '_' BEHAV_temp '.csv'])
                    end
                    
                end
                
            case 'dichotomy effects'
                
                EVENT_temp = EVENT(vEVENTi);
                
                % for each LABEL x ROI
                for vlabi = 1:length(vLABELi)
                    LABEL_temp = LABEL{vLABELi(vlabi)};
                    
                    vec_len = (2^max(N_allclass(vLABELi(vlabi)))-2)/2;
                    if contains(LABEL_temp, 'Goal(6,7,8) x Unc')
                        goal_ef = [9 18 27]; % including Goal(6,7,8) dichotomies only
                        UC_ef = 7; % including the UC_G dichotomy only
                        nonlinear_it = [10 12 17 20 29 30 14 21 28]; % including nonlinear interactions only
                        linear_it = 1:vec_len;
                        linear_it([9 18 27 7 10 12 17 20 29 30 14 21 28]) = []; % including linear interactions only
                        
%                         targ_msr_vec = indiv_msr_vec;
                        targ_msr_vec = nan(size(indiv_msr_vec, 1), ...
                            size(indiv_msr_vec, 2), ...
                            size(indiv_msr_vec, 3), ...
                            size(indiv_msr_vec, 4), ...
                            4);
                        targ_msr_vec(:,:,:,:,1) = nanmean(indiv_msr_vec(:,:,:,:,goal_ef), 5);
                        targ_msr_vec(:,:,:,:,2) = nanmean(indiv_msr_vec(:,:,:,:,UC_ef), 5);
                        targ_msr_vec(:,:,:,:,3) = nanmean(indiv_msr_vec(:,:,:,:,nonlinear_it), 5);
                        targ_msr_vec(:,:,:,:,4) = nanmean(indiv_msr_vec(:,:,:,:,linear_it), 5);
                    else
                        targ_msr_vec = indiv_msr_vec;
                    end
                    
                    for vroii = 1:length(vROIi)
                        ROI_temp = ROI{vROIi(vroii)};
                        
                        figure('Name', ['Event: ' ev_title])
                    
                        % design matrix
                        X = squeeze(nanmean( ...
                            targ_msr_vec(vIDi, vLABELi(vlabi), vROIi(vroii), vEVENTi, :), 4)); % [N_sbj x N_dichotomies]
                        
                        for bhvi = 1:length(BEHAV)
                            bhvs = BHVs{bhvi};
                            BEHAV_temp = BEHAV{bhvi};
                            y = bhvs(vIDi)';
                            
                            tbl = array2table([X, y], 'VariableNames', cat(2, {'Goal', 'UC', 'nonlinearIT', 'linearIT'}, {'y'})); 
                            temp_format = join({'Goal', 'UC', 'nonlinearIT', 'linearIT'}, '+');
                            mdl_format = ['y~' temp_format{:}];
                            mdl = fitlm(tbl, mdl_format, 'RobustOpts', 'ols');
                            
%                             tbl = array2table([X, y]);
%                             mdl = fitlm(tbl, 'RobustOpts', 'ols');
                          
                            disp(mdl)
                            subplot(ceil(length(BEHAV)/bhv_horz), bhv_horz, bhvi)
                            plotAdded(mdl)
                            
                            % save results
                            try
                                save(['data/ccgp_result/corr_behav_dichotomy/' erase(num2str(EVENT_temp), ' ') '/' Exp '_' msr_name '_' LABEL_temp '_' ROI_temp '_' BEHAV_temp bhv_sfx], ...
                                    'mdl', 'tbl', 'X', 'y', 'BEHAV_temp', 'bhv_conds_cell', 'Exp', 'msr_name', 'vIDi', 'LABEL_temp', 'ROI_temp', 'EVENT', 'vEVENTi')
                            catch
                                mkdir(['data/ccgp_result/corr_behav_dichotomy/' erase(num2str(EVENT_temp), ' ')])
                                save(['data/ccgp_result/corr_behav_dichotomy/' erase(num2str(EVENT_temp), ' ') '/' Exp '_' msr_name '_' LABEL_temp '_' ROI_temp '_' BEHAV_temp bhv_sfx], ...
                                    'mdl', 'tbl', 'X', 'y', 'BEHAV_temp', 'bhv_conds_cell', 'Exp', 'msr_name', 'vIDi', 'LABEL_temp', 'ROI_temp', 'EVENT', 'vEVENTi')
                            end
                            
                        end
                        
                    end
                end
                
        end

        
end

disp('step4 done')

%% 

% end % for vs BHV per trial event


disp('END:ccgp_result')