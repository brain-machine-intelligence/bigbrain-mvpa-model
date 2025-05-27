
clear all

%% Set working directory

base_path = fullfile('..'); % If your current working directory is 'demo', use this. Otherwise, set your own path.
addpath(fullfile(base_path, 'functions', 'utils'))

% If your current working directory is 'demo', use this. Otherwise, set your own path that includes 'data'.
cd('..'); 


%% load setting 

% ===== Experiment =====
Exp = 'Lee2014';  N_sbj = 20; Exp_sfx = [];


% ===== ROI mask for multi-voxel pattern (MVP) =====
ROI = {'F3TL', 'F3TR', 'DLPFCL', 'DLPFCR', 'OFC', 'ACCL', 'ACCR', 'preSMAL', 'preSMAR', 'V1L', 'V1R', 'HIPPOL', 'HIPPOR', 'VentStrL', 'VentStrR'};
% ROI = {'F3TL', 'F3TR', 'DLPFCL', 'DLPFCR', 'OFC', 'V1L', 'V1R', 'HIPPOL', 'HIPPOR', 'VentStrL', 'VentStrR'};


% ====================== labelling by task variable =======================
% The labels for goal CCGP across uncertainty, goal SD in low uncertainty
% (uL), and goal SD in high uncertianty (uH), respectively
LABEL = {'Goal(6,7,8)', 'Goal(6,7,8) uL', 'Goal(6,7,8) uH'}; N_allclass = 3 * ones(1, length(LABEL)); LABEL_temp = {'CCGP', 'low uncert. SD', 'high uncert. SD'};


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
% measure = 'accs_rr';
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
measure_cell = {'CCGP', 'acc_mean', 'acc_mean'}; 
% measure_cell = {'CCGP', 'CCGP'}; 
% measure_cell = {'CCGP'}; 
% measure_cell = {};

% For shattering CCGP
% CClabel_name = [];
% CClabel_name = 'CmplxCond';
CClabel_name = 'UncCond'; ccSelect = [];
% CClabel_name = 'Goal(6,7,8)';



% ========================= save setting =========================
% --- shattering ---
save_path = 'data/ccgp_result';

if isempty(measure_cell)
    if contains(measure, 'diff')
        save_name = {['fitclinear_accbox' Exp_sfx], ['fitclinear_accbox_' CClabel_name ccSelect 'CCGP' Exp_sfx]};
        measure_cell = {};
    elseif contains(measure, 'cosine') || contains(measure, 'corr_sig')
%         save_name = ['fitclinear_weights' Exp_sfx];
        save_name = ['fitclinear_accbox' Exp_sfx];
    else
        if isempty(CClabel_name) % SD
%             save_name = ['fitclinear_weights' Exp_sfx];
            save_name = ['fitclinear_accbox' Exp_sfx];
        else % CCGP
            save_name = ['fitclinear_accbox_' CClabel_name ccSelect 'CCGP' Exp_sfx];
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
                save_name{mi} = ['fitclinear_accbox_' CClabel_name ccSelect 'CCGP' Exp_sfx];
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

                            elseif strcmp(measure, 'accs_rr')
                                msr_vec = data_load.save_info.acc_box. ...
                                    ShatteredClassesNumber{N_allclass(labi)}.accs_rr; % (accs_rr: nchoosek(c,m), (2^m - 2)/2, Repeat, rRepeat)
                                msr_vec = squeeze(mean(msr_vec, [3 4])); % (msr_vec: nchoosek(c,m), (2^m - 2)/2)
                                mean_msr_box(id, labi, roii, evi) = mean(msr_vec(:));
                                indiv_msr_vec(id, labi, roii, evi, 1:vec_len) = msr_vec;

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

%% visualization setting

LoadMyColor;
LoadDefaultSettings;
region_color_key = [];
% region_color_key = [4, ...  % ilPFC
%                     9, ...  % V1
%                     16,...  % HIPPO
%                     19];    % VentStr
% axis square; box off;
% set(0,'DefaultAxesFontSize',16,'DefaultTextFontSize',16);
set(0,'DefaultAxesFontSize',14,'DefaultTextFontSize',14);
% set(0,'DefaultAxesFontSize',13,'DefaultTextFontSize',13);
% set(0,'DefaultAxesFontSize',10,'DefaultTextFontSize',10);
% set(0,'DefaultAxesFontSize',10,'DefaultTextFontSize',8);

% ===== visualization mode =====
vis_mode = 'regional_info';
% vis_mode = 'temporal_info';     set(0,'DefaultAxesFontSize',12,'DefaultTextFontSize',12);
% vis_mode = 'corr_behav';    set(0,'DefaultAxesFontSize',10,'DefaultTextFontSize',10);
% vis_mode = 'exhaust_classif';       set(0,'DefaultAxesFontSize',14,'DefaultTextFontSize',14);
region_sort = 1;
% region_sort = 0;
BHV_load = 1;
% BHV_load = 0;
% scatter_color = 1;
scatter_color = 0;
disp_clsf_N = inf;
% disp_clsf_N = 30;
fprintf('----------- %s ----------- \n', vis_mode)

% ===== visualized id index =====
vIDi = 1:N_sbj;

% ===== visualized label index =====
vLABELi = 1:length(LABEL);
% vLABELi = 1;
% vLABELi = 2;
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
% vLABELi = [1 3];
% vLABELi = [2 3];
% vLABELi = [1 2];
% vLABELi = [3 4];
% vLABELi = [2 3];
% vLABELi = [5 6];
% vLABELi = [6 7];
% vLABELi = [7 8];

% ===== visualized region index =====
vROIi = 1:length(ROI);
% vROIi = 1:2;
% vROIi = 1;
% vROIi = [1 2];
% vROIi = [1 3];
% vROIi = [1 4];
% vROIi = [1 2 4];
% vROIi = [2 3];
% vROIi = [2 4];
% vROIi = 1:4;

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
    
    % --- S3 (1 2 4 5 7 8 9) ---
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
    % vEVENTi = 2:6;
%     vEVENTi = 2:7;
%     vEVENTi = 2:5;
%     vEVENTi = 1:6;
%     vEVENTi = 1:7;

    % --- S3 (1 2 4 5 7 8) ---
    vEVENTi = 2:7;
%     vEVENTi = 2:6;
%     vEVENTi = 2:5;
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

% BEHAV = {'PMB'};
% BEHAV = {'SPE'};
% BEHAV = {'RPE'};
% BEHAV = {'Qarb'};
% BEHAV = {'Hit', 'R', 'ChoOpt1', 'ChoOpt2'};

% BEHAV = {'Hit', 'R', 'ChoOpt2'};
% BEHAV = {'PMB', 'R', 'ChoOpt', 'ChoOpt2'};
% BEHAV = {'PMB', 'R', 'PMBvar', 'relMB'};
% BEHAV = {'PMB', 'R', 'ChoOpt', 'ChoOpt2'};
% BEHAV = {'R', 'ChoOpt', 'PMB', 'PMBvar'};

% BEHAV = {'PMB', 'R', 'ChoOpt1', 'ChoOpt2', 'ChoConsist1', 'ChoConsist2'};
% BEHAV = {'PMB', 'R', 'ChoOpt1', 'ChoOpt2', 'ChoCon1', 'ChoCon2'};
% BEHAV = {'Hit', 'R', 'ChoOpt1n2', 'ChoOpt2'};
% BEHAV = {'Hit', 'R', 'ChoOpt1n2', 'ChoConsist1n2'};



% % --- exhaustive BHV measures ---
switch Exp
    case 'Lee2014'
        % 2022-06-01 new
        BEHAV = {'ChoOpt1n2'}; % (C)
%         BEHAV = {'ChoOpt1', 'ChoOpt2', 'ChoOpt1n2', 'PMB' ...
%             'ChoiceSwitch1(GoalSwitch)', 'ChoiceSwitch2(GoalSwitch)', 'ChoiceSwitch1n2(GoalSwitch)', 'lr', ...
%             'ChoiceSwitch1(GoalStay)', 'ChoiceSwitch2(GoalStay)', 'ChoiceSwitch1n2(GoalStay)', 'StrictChoCon', ...
%             'ChoiceSwitch1(GoalSw-St)', 'ChoiceSwitch2(GoalSw-St)', 'ChoiceSwitch1n2(GoalSw-St)', 'invTemp'}; % for Lee 2014
%         BEHAV = {'ChoOpt1', 'ChoOpt2', 'ChoOpt1n2', 'PMB' ...
%             'ChoOptVar1', 'ChoOptVar2', 'ChoOptVar', 'lr', ...
%             'ChoConsist1', 'ChoConsist2', 'ChoConsist1n2', 'StrictChoCon', ...
%             'ChoiceSwitch1', 'ChoiceSwitch2', 'ChoiceSwitch1n2', 'invTemp'}; % for Lee 2014
        % original
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
        BEHAV = {'ChoOpt1', 'ChoOpt2', 'ChoOpt1n2', 'PMB' ...
            'ChoOptVar1', 'ChoOptVar2', 'ChoOptVar', 'PMBvar', ...
            'ChoCon1 ignoring context', 'ChoCon2 ignoring context', 'ChoCon ignoring context'}; % for Kim 2019
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

% bhv_cond = []; bhv_sfx = [];   % all trials
% bhv_cond = 1; bhv_sfx = ' in G';  % only G trials
% bhv_cond = 2; bhv_sfx = ' in H';  % only H trials
% bhv_cond = 3; bhv_sfx = 'in uL';  % only uL trials
% bhv_cond = 4; bhv_sfx = 'in uH';  % only uH trials
bhv_cond = 'G'; bhv_sfx = ' (G)';   % G trials
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

% bhv_horz = 5;
bhv_horz = 4;
% bhv_horz = 3;
% bhv_horz = 2;
% bhv_horz = 1;

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
    
elseif strcmp(measure, 'accs_rr')

    y_lower = chance_lvl - 0.01;
    y_upper = chance_lvl + 0.01;
    % y_upper = chance_lvl + 0.02;
    % y_upper = chance_lvl + 0.055;

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
%     y_upper = chance_lvl + 0.04;
%     y_upper = chance_lvl + 0.05;
    y_upper = chance_lvl + 0.06;
%     y_upper = chance_lvl + 0.075;
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
       
else 
    
    y_lower = [];
    y_upper = [];
    x_lower = [];
    x_upper = [];
    
end

% ===== significance level =====
sgnf_lvl = 0.05;

% ===== measure name =====
if contains(measure, 'diff')
    % 2021-11-17 temp
    msr_name = ['SD - CCGP across ' CClabel_name];
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
%         msr_name = 'Separability dimension';
        msr_name = [];
    else
        msr_name = ['CCGP across ' CClabel_name];
    end
elseif contains(measure, 'separable')
    msr_name = 'Sensitivity dimension';
elseif contains(measure, 'separable2')
    msr_name = 'Sensitivity dimension';
elseif contains(measure, 'CCGP')
    if isempty(CClabel_name)
%         msr_name = 'Separability dimension';
        msr_name = [];
    else
        msr_name = ['CCGP across ' CClabel_name];
    end
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
else
    msr_name = measure;
end


disp('step3 done')

%% visualization

close all

switch vis_mode
    
    % ---------------------------------------------------------------------
    case 'regional_info'
        
        % averageing event
        msrmat = squeeze(mean(mean_msr_box(vIDi, vLABELi, vROIi, vEVENTi), 4));
        msrmat_sbjmean = squeeze(mean(msrmat, 1)); % LEN_LABEL x LEN_ROI
        msrmat_sbjste = squeeze(std(msrmat, 0, 1) ... 
            ./sqrt(sum(~isnan(msrmat), 1))); % LEN_LABEL x LEN_ROI

        % Single label vs chance level
        if isscalar(vLABELi)
            
            % Significance check & marking
            [h, p] = ttest(msrmat - chance_lvl(vLABELi)); % msrmat: [N_sbj x N_roi]
            colors_temp = 0.4 * ones(length(vROIi), 3); 
            if any(h)
                colors_temp(logical(h), :) = repmat([0.7 0.2 0.3], sum(logical(h)), 1); 
            end
            
            % record significant event
            sig_roi = nan * ones(1, length(vROIi));
            sig_text = cell(1, length(vROIi));
            for vroii = 1:length(vROIi)
%                 [h, p] = ttest(msrmat1(:,vroii)', msrmat2(:,vroii)', ...
%                     'Alpha', sgnf_lvl);
                [h, p, ~, stats] = ttest(msrmat(:,vroii), chance_lvl(vLABELi), 'Alpha', sgnf_lvl);
                if h
                    sig_roi(vroii) = h; 
                    fprintf(ROI{vROIi(vroii)})
                    fprintf(' p=%.2e \n', p)
                    disp(stats)
                end
                if p < 0.0001
                    sig_text{vroii} = '****';
                elseif p < 0.001
                    sig_text{vroii} = '***';
                elseif p < 0.01
                    sig_text{vroii} = '**';
                elseif p < 0.05
                    sig_text{vroii} = '*';
                end
            end
                        
            % Visualization
            figure
            plot(msrmat', 'LineWidth', 0.01, 'Marker', 'o', 'Color', [.9 .9 .9]); hold on
            text(find(~isnan(sig_roi)), ...
                max(msrmat(:, ~isnan(sig_roi)), [], 1), ...
                sig_text(~isnan(sig_roi)), 'FontSize', 20, ...
                'HorizontalAlignment', 'center')
            boxplot(msrmat, 'Colors', colors_temp); hold off
            title({[msr_name ' (red:sig, gray:n.s.)'], ...
                [LABEL{vLABELi} ' (' ev_title ')'], ... 
                ['N(>chance): ' num2str(sum(msrmat > chance_lvl(vLABELi), 1))]}, ...
                'Interpreter', 'none')
            xticklabels(ROI(vROIi))
            % if length(vROIi) > 6; xtickangle(45); end
            
        end
        
        % visualization
        figure
        % subplot(2, length(vROIi), 1:length(vROIi))
        B = bar(msrmat_sbjmean', 'EdgeColor', 'none'); axis square; box off;
        hold on
        
        if length(vLABELi) == 3
            % diff_erbr = .175;
            % diff_erbr = .2; % rldm 2022
            diff_erbr = .225;
            % diff_erbr = .25;
        elseif isscalar(vLABELi)
            diff_erbr = 0;
        else 
            diff_erbr = .2;
        end
        xdiff_erbr = linspace(-diff_erbr, diff_erbr, length(vLABELi));
        xtick_erbr = reshape((xdiff_erbr)' + (1:length(vROIi)), 1, []);     % LAB x REG -> 1 * (LAB * REG)
        E = errorbar(xtick_erbr, reshape(msrmat_sbjmean,1,[]), ... 
            reshape(msrmat_sbjste,1,[]), ...
            'LineStyle', 'none', ...
            'LineWidth', 1.5, ...
            'Color', 'red'); % rldm 2022
        
%         E.Color = colors{21}; % rldm 2022

        
        % Draw the chance level line
        line([0 length(vROIi)+1], [chance_lvl(1) chance_lvl(1)], ...
            'Color', [.1 .1 .1], ...
            'LineStyle', "--", ...
            'LineWidth', 1.5)
        
        hold off
        
        title(ev_title)
%         if iscell(EVENT)
%             title(['Event ' num2str([EVENT{vEVENTi}])])
%         else
%             title(['Event ' num2str(EVENT(vEVENTi))])
%         end
        xticks(1:length(vROIi))
        % xtickangle(0)
        xticklabels(ROI(vROIi))
        if ~isempty(y_lower); ylim([min(y_lower) max(y_upper)]); end
%         ylabel(['Mean ' measure], 'Interpreter', 'none')
        ylabel(msr_name, 'Interpreter', 'none')
        % rldm 2022 flag
        B(1).FaceColor = [.5 .5 .5];
        try
            B(2).FaceColor = [1 1 1];
            B(2).EdgeColor = [0 0 0];
            B(3).FaceColor = [0 0 0];
        catch
        end
        try
            legend(B(:), LABEL{vLABELi})
        catch
            legend(LABEL{vLABELi})
        end

        % ===== significance check B =====
        if length(vLABELi) > 1

            for vroii = 1:length(vROIi)

                msrmat_roi = squeeze(msrmat(:, :, vroii));      % N_sbj x LEN_LABEL

                lab_pair_sig = ones(length(vLABELi), length(vLABELi));
                lab_pairs = nchoosek(1:length(vLABELi), 2);
                %                 pvals = nan * ones(size(lab_pairs, 1), 1);
                pvals = cell(size(lab_pairs, 1), 1);
                disp(ROI{vROIi(vroii)})
                for pairi = 1:size(lab_pairs, 1)
                    lab_pair = lab_pairs(pairi, :);
                    % pairwise t-test for 2 labels
                    [h, p, ~, stats] = ttest(msrmat_roi(:, lab_pair(1)), ...
                        msrmat_roi(:, lab_pair(2)), ...
                        'Alpha', sgnf_lvl);
                    % pvals(pairi) = floor(p*10000)/10000;
                    %                     pvals{pairi} = sprintf('%.2e', p);
                    if p < 0.0001
                        pvals{pairi} = '****';
                    elseif p < 0.001
                        pvals{pairi} = '***';
                    elseif p < 0.01
                        pvals{pairi} = '**';
                    elseif p < 0.05
                        pvals{pairi} = '*';
                    else
                        pvals{pairi} = 'n.s.';
                    end
                    if h; lab_pair_sig(lab_pair(1), lab_pair(2)) = p; end

                    % To print statistic and p value
                    fprintf([LABEL_temp{vLABELi(lab_pair(1))} ' vs ' LABEL_temp{vLABELi(lab_pair(2))} ...
                        ': t(%d)=%.3f\t p=%.3e (' pvals{pairi} ')\n'], stats.df, stats.tstat, p)
                end
                disp(' ')

            end

        end







        % % ===== BOOTSTRAPPING/SUBSAMPLING FOR RESULT ROBUSTNESS =====
        % % subsetSize = 20;    % size of each subject subset
        % % numBoot    = 1000;  % number of bootstrap iterations
        % % jackknife
        % subsetSize = 19;    % size of each subject subset
        % numBoot    = 20;  % number of bootstrap iterations
        % 
        % nSbjv = length(vIDi);
        % nLabv = length(vLABELi);
        % nROIv = length(vROIi);
        % 
        % % visDiff = [];
        % visDiff = {1, [2 3]};
        % 
        % % --- preallocate effect-size distribution ---
        % if ~isempty(visDiff)
        %     nLabv = length(visDiff{1}) * length(visDiff{2});
        %     vLABELi = 1:nLabv;
        % end
        % matBoot = nan(numBoot, nLabv, nROIv);
        % sigBoot = nan(numBoot, nLabv, nROIv);
        % labPair = nchoosek(1:nLabv, 2);
        % nPair = size(labPair, 1);
        % pairSig = false(numBoot, nPair, nROIv);
        % 
        % for b = 1:numBoot
        % 
        %     rng(b);
        % 
        %     % --- sample a random subset of subjects without replacement ---
        %     if subsetSize == 19 && numBoot == 20
        %         idxSub = 1:numBoot;
        %         idxSub(b) = [];
        %     elseif subsetSize == 20
        %         idxSub = randsample(nSbjv, subsetSize, true);
        %     else
        %         idxSub = randsample(nSbjv, subsetSize);
        %     end
        % 
        %     % --- compute mean difference and pooled std for each ROI ---
        %     % single-label case: msrmat(:, vroii) is [nSbjv x 1]
        %     bootData = msrmat(idxSub, :, :);          % [subsetSize x N_label x N_roi]
        %     bootMean = squeeze(mean(bootData, 1));    % [N_label x N_roi]
        %     chanceVec = chance_lvl(vLABELi)';         % [nLabv × 1]
        %     if ~isempty(visDiff)
        %         bootDiff  = bootMean(visDiff{1}, :) - bootMean(visDiff{2}, :); % (nDiff, nROI)
        %         diffs = bootData(:, visDiff{1}, :) - bootData(:, visDiff{2}, :); % (subsetSize, nDiff, nROI)
        %     else
        %         bootDiff  = bootMean - chanceVec;         % [nLabv × nROI]
        %         chanceRep = reshape(chanceVec, [1, length(vLABELi), 1]);
        %         diffs     = bootData - chanceRep;         % [subsetSize × nLabv × nROI]
        %     end
        % 
        %     for labi = 1:nLabv
        %         for roii = 1:nROIv
        %             data_bs = squeeze(diffs(:, labi, roii));  % [subsetSize × 1]
        %             sigBoot(b, labi, roii) = ttest(data_bs, 0, 'Alpha', sgnf_lvl);
        %         end
        %     end
        % 
        %     % % ---------------------------------------------
        %     % % Pairwise t-tests between labels
        %     % % ---------------------------------------------
        %     % for roii = 1:nROIv
        %     %     for i = 1:nPair
        %     %         pair1 = labPair(i, 1);
        %     %         pair2 = labPair(i, 2);
        %     %         msrVec1 = bootData(:, pair1,   roii);  % subsetSize x 1
        %     %         msrVec2 = bootData(:, pair2,   roii);  % subsetSize x 1
        %     %         [h, pVal] = ttest(msrVec1, msrVec2, 'Alpha', sgnf_lvl);
        %     %         pairSig(b, i, roii) = h;  % 1 if significant
        %     %     end
        %     % end
        % 
        %     % --- then your accuracy-chance for each label × ROI ---
        %     matBoot(b,:,:) = bootDiff;      % store into a 3-D array now
        % 
        %     % % --- then your Cohen's d for each label × ROI ---
        %     % effDist(b,:,:) = bootDiff ./ bootStd;      % store into a 3-D array now
        % end
        % 
        % % --- modified: statistical test per label and ROI ---
        % 
        % % h_bs   = false(nLabv, nROIv);
        % % p_bs   = nan(nLabv, nROIv);
        % % 
        % % for labi = 1:nLabv
        % %     for roii = 1:nROIv
        % %         % extract bootstrap distribution for this label & ROI
        % %         data_bs = squeeze(meanBoot(:, labi, roii));  % [numBoot × 1]
        % %         [h_bs(labi, roii), p_bs(labi, roii)] = ttest(data_bs, 0, 'Alpha', sgnf_lvl);
        % %     end
        % % end
        % 
        % 
        % % ===== visualize effect‐size distributions per ROI in one figure =====
        % % figure('Units','normalized','Position',[.1 .1 .8 .8]);
        % figure;
        % for roii = 1:nROIv
        %     % subplot(2, nROIv/2, k);
        % 
        %     data_bs = squeeze(matBoot(:, :, roii)); % [numBoot × nLabv]
        %     sig_bs = logical(squeeze(sigBoot(:, :, roii)));
        % 
        %     nSig = sum(sig_bs, 1); % 1 x nLabv
        %     beeX1 = []; beeY1 = []; % for significant result
        %     beeX2 = []; beeY2 = []; % for nonsignificant result
        %     if ~isempty(visDiff)
        %         colors_temp = repmat([.85 .5 .5], nLabv, 1);
        %     else
        %         colors_temp = [.5 .5 .5; .9 .9 .9; 0 0 0];
        %     end
        %     for labi = 1:nLabv
        %         tempX1 = labi * ones(nSig(labi), 1);
        %         tempX2 = labi * ones(numBoot - nSig(labi), 1);
        %         tempY1 = data_bs(sig_bs(:, labi), labi);
        %         tempY2 = data_bs(~sig_bs(:, labi), labi);
        %         beeX1 = cat(1, beeX1, tempX1);
        %         beeX2 = cat(1, beeX2, tempX2);
        %         beeY1 = cat(1, beeY1, tempY1);
        %         beeY2 = cat(1, beeY2, tempY2);
        %     end
        %     colors_temp(nSig==0, :) = [];
        % 
        %     % plot over labels
        %     subplot(1, nROIv, roii);
        % 
        %     % siginificant results
        %     beeswarm(beeX1, beeY1, ...
        %         'colormap', repmat(colors_temp, sum(nSig>0), 1), ...
        %         'overlay_style', false, ...
        %         'corral_style', 'random', ...
        %         'MarkerFaceAlpha', .5);
        %     hold on
        %     % non-significant results
        %     beeswarm(beeX2, beeY2, ...
        %         'colormap', [.8 .8 .8; .2 .2 .2], ...
        %         'overlay_style', false, ...
        %         'corral_style', 'random');
        %     % beeswarm(beeX2, beeY2, ...
        %     %     'colormap', repmat([.5 .5 .5], sum((numBoot-nSig)>0), 1), ...
        %     %     'overlay_style', false, ...
        %     %     'corral_style', 'random');
        % 
        %     % all results
        %     if subsetSize == 19 && numBoot == 20
        %         beeswarm(reshape(repmat(1:nLabv, numBoot, 1), [], 1), ...
        %             reshape(data_bs, [], 1), ...
        %             'overlay_style', 'cik', ...
        %             'corral_style', 'gutter', ...
        %             'MarkerFaceColor', 'none');
        %     else
        %         beeswarm(reshape(repmat(1:nLabv, numBoot, 1), [], 1), ...
        %             reshape(data_bs, [], 1), ...
        %             'overlay_style', 'box2', ...
        %             'corral_style', 'gutter', ...
        %             'MarkerFaceColor', 'none');
        %     end
        %     hold on;
        % 
        %     % ---------------------------------------------
        %     % Visualize pairwise significance with lines
        %     % ---------------------------------------------
        %     % Recall: pairSig(b, i, roii) => 1 if label i vs label 4 is significant in iteration b
        %     % We'll loop over each iteration b, and for each i=1..3 that is significant,
        %     % draw a gray line from the (x=i, y=data_bs(b,i)) to (x=4, y=data_bs(b,4)).
        %     for b = 1:numBoot
        %         for i = 1:nPair
        %             if pairSig(b, i, roii) == 1
        %                 % Draw a line
        %                 xCoord = labPair(i, :);
        %                 yCoord = [data_bs(b, labPair(i, 1)), data_bs(b, labPair(i, 2))];
        %                 plot(xCoord, yCoord, ':', 'Color', [0.8 0.8 0.8], 'LineWidth', .3);
        %             end
        %         end
        %     end
        % 
        %     % overlay original mean accuracy - chance
        %     % msrmat_sbjmean: (nLabv, nROI)
        %     if ~isempty(visDiff)
        %         plot(1:nLabv, msrmat_sbjmean(visDiff{1}, roii) - msrmat_sbjmean(visDiff{2}, roii), ...
        %             'w_', 'LineWidth', 1.4, 'MarkerSize', 15);
        %     else
        %         plot(1:nLabv, msrmat_sbjmean(:, roii) - chance_lvl(vLABELi)', ...
        %             'w_', 'LineWidth', 1.4, 'MarkerSize', 15);
        %     end
        % 
        %     % % overlay mean Cohen's d per label
        %     % mns = mean(data_bs, 1);
        %     % plot(1:nLabv, mns, 'k*', 'LineWidth', 1.5, 'MarkerSize', 8);
        % 
        %     % % mark labels where the bootstrap test is significant
        %     % for lab = 1:nLabv
        %     %     if h_bs(lab, k)
        %     %         y = max(data_bs(:, lab));
        %     %         text(lab, y, '*', ...
        %     %             'FontSize', 16, ...
        %     %             'HorizontalAlignment', 'center', ...
        %     %             'Color', [0.7 0.2 0.3]);
        %     %     end
        %     % end
        % 
        %     line([0 nLabv+1], [0 0], ...
        %         'Color', [.5 .5 .5], ...
        %         'LineStyle', "--", ...
        %         'LineWidth', 1.5)
        % 
        %     title(ROI{vROIi(roii)}, 'Interpreter', 'none');
        %     xticks(1:length(vLABELi))
        %     if ~isempty(visDiff)
        %         ylabel('CCGP - SD');
        %         ylim([-.01 .01])
        %         xticklabels({'low', 'high'})
        %         xlabel('Uncertainty')
        %     else
        %         ylabel('CCGP - chance');
        %         switch subsetSize
        %             case 5
        %                 ylim([-.03 .08])
        %             case 10
        %                 ylim([-.02 .06])
        %             case 15
        %                 ylim([-.01 .06])
        %             case 19
        %                 ylim([-.01 .05])
        %             case 20
        %                 ylim([-.02 .06])
        %         end
        %         xlabel('CCGP', 'SD low', 'SD high')
        %     end
        %     hold off;
        % end
        % 
        % % add an overall title
        % if subsetSize == 19 && numBoot == 20
        %     sgtitle(sprintf('Jackknife resampling result (subset=%d, iterations=%d)', ...
        %         subsetSize, numBoot), 'FontSize', 14);
        % elseif subsetSize == 20
        %     sgtitle(sprintf('Bootstrapped effect sizes per ROI (subset=%d, iterations=%d)', ...
        %         subsetSize, numBoot), 'FontSize', 14);
        % else
        %     sgtitle(sprintf('Subsampled effect sizes per ROI (subset=%d, iterations=%d)', ...
        %         subsetSize, numBoot), 'FontSize', 14);
        % end







        
        
end

disp('step4 done')

%% 
disp('END:ccgp_result')

