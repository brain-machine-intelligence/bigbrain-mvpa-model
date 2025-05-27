
clear all

%% Set working directory

base_path = fullfile('..'); % If your current working directory is 'demo', use this. Otherwise, set your own path.
addpath(fullfile(base_path, 'functions', 'utils'))

% If your current working directory is 'demo', use this. Otherwise, set your own path that includes 'data'.
cd('..'); 


%% load setting 

% ===== experiment =====
Exp = 'Lee2014';  N_sbj = 20; Exp_sfx = [];
SepClass = [];


% ===== ROI mask for multi-voxel pattern (MVP) =====
% ROI = {'F3TL', 'F3TR', 'DLPFCL', 'DLPFCR', 'OFC', 'V1L', 'V1R', 'HIPPOL', 'HIPPOR', 'VentStrL', 'VentStrR'};
% ROI = {'ACCL', 'ACCR', 'SMAL', 'SMAR'};
ROI = {'F3TL', 'F3TR', 'DLPFCL', 'DLPFCR', 'OFC', 'ACCL', 'ACCR', 'preSMAL', 'preSMAR', 'V1L', 'V1R', 'HIPPOL', 'HIPPOR', 'VentStrL', 'VentStrR'};
% ROI = {'F3TL', 'F3TR'};
% ROI = {'DLPFCL', 'DLPFCR'};
% ROI = {'ACCL', 'ACCR'};
% ROI = {'SMAL', 'SMAR'};
% ROI = {'V1L', 'V1R'};
% ROI = {'HIPPOL', 'HIPPOR'};
% ROI = {'VentStrL', 'VentStrR'};


% ====================== labelling by task variable =======================
LABEL = {'Goal(6,7,8)', 'S2', 'S3'}; N_allclass = [3 4 4];


%% load setting 2


% ========================= MVP setting =========================
% (original main set)
% ROI_types = {'AAL3', 'AAL3', 'AAL', 'AAL', 'AAL', 'AAL3', 'AAL3', 'AAL3', 'AAL3', 'AAL3', 'AAL3'}; % R/L integrated OFC
% (full set)
ROI_types = {'AAL3', 'AAL3', 'AAL', 'AAL', 'AAL', 'AAL3', 'AAL3', 'JuBrain', 'JuBrain', 'AAL3', 'AAL3', 'AAL3', 'AAL3', 'AAL3', 'AAL3'};
% (others)
% ROI_types = {'AAL3', 'AAL3', 'AAL3', 'AAL3'};
% ROI_types = {'AAL3', 'AAL3', 'AAL', 'AAL', 'AAL', 'AAL3', 'AAL3', 'AAL3', 'AAL3', 'AAL3', 'AAL3', 'AAL3', 'AAL3', 'AAL3', 'AAL3'};
% ROI_types = {'AAL3', 'AAL3', 'AAL3', 'AAL3', 'AAL3', 'AAL3', 'AAL3', 'AAL3', 'AAL3', 'AAL3', 'AAL3', 'AAL3'};
% ROI_types = {'AAL3', 'AAL3'};
% ROI_types = {'AAL', 'AAL'};


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

% Single measure
measure = 'accs'; measure_cell = {}; CClabel_name = []; % decoding accuracy
% measure = 'accs_r'; measure_cell = {}; CClabel_name = []; % decoding accuracy
% measure = 'acc_mean'; measure_cell = {}; CClabel_name = []; % separability dimension
% measure = 'acc_mean'; measure_cell = {}; CClabel_name = 'UncCond'; % cross-condition generalization performance
% measure = 'cosine'; measure_cell = {}; CClabel_name = []; % cosine similarity


% ========================= save setting =========================

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

% --- decoding ---
save_path = 'data/decoding_result';

% save_name = ['fitrlinear_corrcoef_' Exp_sfx];
% save_name = ['fitrlinear_corrcoef_' Exp_sfx 'nullS3_'];

save_name = ['fitcecoc_CVacc_' Exp_sfx];

% save_name = ['fitcecoc_CCGP_' Exp_sfx trained_mdl_info];
% % % % save_name = ['fitcecoc_CCGP_Kim_wa_' trained_mdl_info];

% save_name = {['fitcecoc_CVacc_' Exp_sfx], ['fitcecoc_CCGP_' Exp_sfx trained_mdl_info]};
% % % % save_name = {'fitcecoc_CVacc_Kim_wa_', ['fitcecoc_CCGP_Kim_wa_' trained_mdl_info]};

% save_name = ['fitcecoc_weights_' Exp_sfx];


if contains(save_path, 'shattering')
    chance_lvl = 0.5 * ones(size(N_allclass));
elseif contains(save_path, 'decoding')
    if contains(save_name, 'fitrlinear')
        chance_lvl = 0./N_allclass;
    else
        chance_lvl = 1./N_allclass;
    end
end

if (contains(measure, 'diff') || contains(measure, 'cosine'))
    chance_lvl = zeros(size(N_allclass));
end

if (contains(measure, 'GI') || contains(measure, 'ratio'))
    chance_lvl = ones(size(N_allclass));
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
elseif contains(measure, 'diff') || contains(measure, 'ratio') 
    measure_ori = measure;
    msr_operation = measure;
elseif isempty(measure_cell)
    measure_cell = {measure};
    msr_operation = [];
elseif contains(measure, 'GI')
    measure_ori = measure;
    msr_operation = 'GI';
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
%     if contains(measure, 'cosine'); vec_len = 5000; end
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
                    
                    % ============ Decoding analysis ============
                    try
                        data_load = load([save_path '/' roi_type '/' ROI{roii} '/' , ...
                            LABEL{labi} '/' num2str(ev_name) '/', ...
                            save_name num2str(id) '.mat']);

                        if isfield(data_load.save_info, 'CCGP') && strcmp(measure, 'CCGP')
                            mean_msr_box(id, labi, roii, evi) = data_load.save_info.CCGP;
                            %                             if evi == 6
                            %                                 fprintf(data_load.save_info.CCGP)
                            %                             end
                        end

                        if isfield(data_load.save_info, 'accs') && strcmp(measure, 'accs')
                            mean_msr_box(id, labi, roii, evi) = mean(data_load.save_info.accs);
                        end

                        if isfield(data_load.save_info, 'beta_foldmeans_r') && strcmp(measure, 'cosine')
                            % 2021-07-23 temporary
                            data_load1 = load([save_path '/' roi_type '/' ROI{roii} '/' , ...
                                LABEL{1} '/' num2str(ev_name) '/', ...
                                save_name num2str(id) '.mat']);
                            data_load2 = load([save_path '/' roi_type '/' ROI{roii} '/' , ...
                                LABEL{2} '/' num2str(ev_name) '/', ...
                                save_name num2str(id) '.mat']);
                            weights1 = mean(data_load1.save_info.beta_foldmeans_r, 2);
                            weights2 = mean(data_load2.save_info.beta_foldmeans_r, 2);
                            norm_weights1 = weights1/vecnorm(weights1);
                            norm_weights2 = weights2/vecnorm(weights2);
                            mean_msr_box(id, labi, roii, evi) = norm_weights1' * norm_weights2;
                        end

                    catch Err
                        disp([num2str(Err.stack.line) ' ' Err.message])
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
        mean_msr_box = mean_msr_cell{1} ./ mean_msr_cell{2};
        measure = measure_ori;
    elseif contains(msr_operation, 'GI')
%         mean_msr_box = (mean_msr_cell{1} - .5) ./ (mean_msr_cell{2} - .5); % CCGP first, SD second
        mean_msr_box = (mean_msr_cell{1} - .4) ./ (mean_msr_cell{2} - .4); % CCGP first, SD second
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
elseif length(original_ROI) == 4
    roii1 = [1 2];
    roii2 = [3 4];
    bilateral_mode = 5;
elseif length(original_ROI) == 15
    roii1 = [1 2];
    roii2 = [3 4];
    roii3 = 5;
    roii4 = [6 7];
    roii5 = [8 9];
    roii6 = [10 11];
    roii7 = [12 13];
    roii8 = [14 15];
    bilateral_mode = 6;
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
            mean_msr_box(:,:,1,:) = mean(temp_msr_box(:, :, roii1, :), 3);
            mean_msr_box(:,:,2,:) = mean(temp_msr_box(:, :, roii2, :), 3);
            ROI = {'ACC', 'SMA'};
            
            indiv_msr_vec(:,:,1,:,:) = mean(temp_msr_vecs(:, :, roii1, :, :), 3);
            indiv_msr_vec(:,:,2,:,:) = mean(temp_msr_vecs(:, :, roii2, :, :), 3);

        case 6
            mean_msr_box(:,:,1,:) = mean(temp_msr_box(:, :, roii1, :), 3);
            mean_msr_box(:,:,2,:) = mean(temp_msr_box(:, :, roii2, :), 3);
            mean_msr_box(:,:,3,:) = mean(temp_msr_box(:, :, roii3, :), 3);
            mean_msr_box(:,:,4,:) = mean(temp_msr_box(:, :, roii4, :), 3);
            mean_msr_box(:,:,5,:) = mean(temp_msr_box(:, :, roii5, :), 3);
            mean_msr_box(:,:,6,:) = mean(temp_msr_box(:, :, roii6, :), 3);
            mean_msr_box(:,:,7,:) = mean(temp_msr_box(:, :, roii7, :), 3);
            mean_msr_box(:,:,8,:) = mean(temp_msr_box(:, :, roii8, :), 3);
            ROI = {'vlPFC', 'dlPFC', 'OFC', 'ACC', 'preSMA', 'V1', 'HPC', 'vStr'};

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
if bilateral_integ && length(ROI) == 4
    region_color_key = [4, ...  % ilPFC
        9, ...  % V1
        16,...  % HIPPO
        19];    % VentStr
else
    region_color_key = [];
end
% axis square; box off;
% set(0,'DefaultAxesFontSize',16,'DefaultTextFontSize',16);
set(0,'DefaultAxesFontSize',14,'DefaultTextFontSize',14);
% set(0,'DefaultAxesFontSize',10,'DefaultTextFontSize',10);
% set(0,'DefaultAxesFontSize',10,'DefaultTextFontSize',8);

% ===== visualization mode =====
% vis_mode = 'regional_info';
vis_mode = 'temporal_info'; temp_mode = 'group';     set(0,'DefaultAxesFontSize',12,'DefaultTextFontSize',12);
% vis_mode = 'temporal_info'; temp_mode = 'individual';     set(0,'DefaultAxesFontSize',12,'DefaultTextFontSize',12);
% vis_mode = 'corr_behav';  corr_mode = 'original';  set(0,'DefaultAxesFontSize',10,'DefaultTextFontSize',10);
% vis_mode = 'corr_behav';  corr_mode = 'ROI effects';  set(0,'DefaultAxesFontSize',10,'DefaultTextFontSize',10);
% vis_mode = 'corr_behav';  corr_mode = 'dichotomy effects';  set(0,'DefaultAxesFontSize',10,'DefaultTextFontSize',10);
% vis_mode = 'exhaust_classif';       set(0,'DefaultAxesFontSize',14,'DefaultTextFontSize',14);
% vis_mode = 'corr_pca'; pca_seed = 2021; set(0,'DefaultAxesFontSize',10,'DefaultTextFontSize',10);

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
% vLABELi = 1:length(LABEL);
vLABELi = 1;
% vLABELi = 2;
% vLABELi = [1 2];
% vLABELi = [2 3];

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
    vEVENTi = 2:7;    
    % vEVENTi = 2:6;
    % vEVENTi = 2:5;
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
%     y_upper = chance_lvl + 0.08; % shattering
    y_upper = chance_lvl + 0.1; % decoding
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
    
    case 'temporal_info'
        
        % ===== group result plot =====
        msrmat = mean_msr_box(vIDi, vLABELi, vROIi, vEVENTi);
        msrmat_size = size(msrmat);
        msrmat_sbjmean = reshape(mean(msrmat, 1), msrmat_size(2:end)); % [N_label x N_roi x N_event]
        msrmat_sbjste = reshape(std(msrmat, 0, 1) ... 
            ./sqrt(sum(~isnan(msrmat), 1)), msrmat_size(2:end)); % [N_label x N_roi x N_event]
        
        % --- A. subplot: roi, legend: label ---
        figure('Position', [0, 1080/2-80, 1920, 1080/2])
        temp_fig_title_sfx = [];
        for vroii = 1:length(vROIi)
            subplot(1, length(vROIi), vroii)
            
            % significance check (only for 2 labels)
            if length(vLABELi) == 2
                msrmat1 = squeeze(msrmat(:,1,vroii,:));
                msrmat2 = squeeze(msrmat(:,2,vroii,:));
                
                % record significant event
                sig_event = nan * ones(1, length(vEVENTi));
                sig_text = cell(1, length(vEVENTi));
                for vevi = 1:length(vEVENTi)                
                    [h, p] = ttest(msrmat1(:,vevi)', msrmat2(:,vevi)', ... 
                        'Alpha', sgnf_lvl);
                    if h; sig_event(vevi) = h; end
                    if p < 0.001
                        sig_text{vevi} = '***'; 
                    elseif p < 0.01
                        sig_text{vevi} = '**';
                    elseif p < 0.05
                        sig_text{vevi} = '*';
                    end
                end
                
                % Significance asterisk
                mean_msr_vec = squeeze(max(msrmat_sbjmean(:, vroii, :), [], 1));
                ste_msr_vec = squeeze(max(msrmat_sbjste(:, vroii, :), [], 1));
                text(find(~isnan(sig_event)), ...
                    mean_msr_vec(~isnan(sig_event))+mean(ste_msr_vec)*1.5, ...
                    sig_text(~isnan(sig_event)), 'FontSize', 20, ...
                    'HorizontalAlignment', 'center')
                
                hold on
            elseif isscalar(vLABELi)
                
                msrmat1 = squeeze(msrmat(:,:,vroii,:)); % N_sbj x N_event
                msrmat2 = chance_lvl(vLABELi) * ones(size(msrmat1));
                
                if isempty(y_upper)
%                     plot(msrmat1', 'LineWidth', 0.01, 'Marker', 'o', 'Color', [.9 .9 .9])
                    plot(msrmat1', 'LineWidth', 0.01, 'Color', [.9 .9 .9]); hold on
                    plot(msrmat1', 'LineStyle', 'none', ... 
                        'Marker', 'o', ... 
                        'MarkerSize', 4, ...
                        'Color', [.7 .7 .7])
                end
                hold on

                % record significant event
                sig_event = nan * ones(1, length(vEVENTi));
                sig_text = cell(1, length(vEVENTi));
                disp(ROI{vROIi(vroii)})
                for vevi = 1:length(vEVENTi)                
                    [h, p, ci, stats] = ttest(msrmat1(:,vevi)', msrmat2(:,vevi)', ...
                        'Alpha', sgnf_lvl);
                    if h; sig_event(vevi) = h; end
                    if p < 0.001
                        sig_text{vevi} = '***'; 
                        sig_event(vevi) = 1;
                    elseif p < 0.01
                        sig_text{vevi} = '**';
                        sig_event(vevi) = 1;
                    elseif p < 0.05
                        sig_text{vevi} = '*';
                        sig_event(vevi) = 1;
                    else
                        sig_text{vevi} = 'n.s.';
                    end

                    % To print statistic and p value
                    fprintf([EVENT_NAME{vevi} ': t(%d)=%.3f\t p=%.3e (' sig_text{vevi} ')\n'], stats.df, stats.tstat, p)

                end
                disp(' ')
                
                % Significance asterisk
                mean_msr_vec = squeeze(msrmat_sbjmean(:, vroii, :));
                ste_msr_vec = squeeze(msrmat_sbjste(:, vroii, :));
                text(find(~isnan(sig_event)), ...
                    mean_msr_vec(~isnan(sig_event))+mean(ste_msr_vec)*2, ...
                    sig_text(~isnan(sig_event)), 'FontSize', 14, ...
                    'HorizontalAlignment', 'center')
                
                temp_fig_title_sfx = [' ' LABEL{vLABELi}];
            end
            
            % errorbar plot
            EB = errorbar(repmat(1:length(vEVENTi), length(vLABELi), 1)', ...
                squeeze(msrmat_sbjmean(:, vroii, :))', ...
                squeeze(msrmat_sbjste(:, vroii, :))', ...
                'LineWidth', 3); % LineWidth = 3 for RLDM 2022 (original setting = 1.5)
%             for vlabi = 1:length(vLABELi)
%                 EB(vlabi).Color = colors{vlabi};
%             end
            EB(1).Color = [77 161 169]/255;
            
            % Draw the chance level line
            if isscalar(vLABELi)
                line([0 length(vEVENTi)], [chance_lvl(vLABELi) chance_lvl(vLABELi)], ...
                    'Color', [.5 .5 .5], ...
                    'LineStyle', "--", ...
                    'LineWidth', 1.5)
            end
            
            axis square; box off
            hold off
            
            % figure info
%             title({ROI{vROIi(vroii)}, temp_fig_title_sfx})
%             title({ROI{vROIi(vroii)}, temp_fig_title_sfx}, 'FontSize', 15)
            title({temp_fig_title_sfx, ROI{vROIi(vroii)}, ...
                }, 'FontSize', 15, 'FontWeight', 'normal')
%             title(ROI{vROIi(vroii)}, 'FontSize', 15)
            xticks(1:length(vEVENTi))
            xticklabels(ev_names); 
            xlim([0 length(vEVENTi)])
%             ax = gca; ax.FontSize = 9;
            if ~isempty(y_lower) && ~isempty(y_upper); ylim([min(y_lower(vLABELi)) max(y_upper(vLABELi))]); end
%             ylabel(['Mean ' measure], 'Interpreter', 'none')
%             ylabel(msr_name)
            
        end
        if length(vLABELi) > 1
            legend(EB(:), LABEL{vLABELi})
        end
        
        
end

disp('step4 done')


%% 
disp('END:shattering_result')

