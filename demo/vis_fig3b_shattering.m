
clear all

%% Set working directory

base_path = fullfile('..'); % If your current working directory is 'demo', use this. Otherwise, set your own path.
addpath(fullfile(base_path, 'functions', 'utils'))

% If your current working directory is 'demo', use this. Otherwise, set your own path that includes 'data'.
cd('..'); 


%% load setting 

% ===== experiment =====
Exp = 'Lee2014'; N_sbj = 20; Exp_sfx = []; 


% ===== ROI mask for multi-voxel pattern (MVP) =====
ROI = {'F3TL', 'F3TR', 'DLPFCL', 'DLPFCR', 'OFC', 'ACCL', 'ACCR', 'preSMAL', 'preSMAR', 'V1L', 'V1R', 'HIPPOL', 'HIPPOR', 'VentStrL', 'VentStrR'};
% ROI = {'F3TL', 'F3TR', 'DLPFCL', 'DLPFCR', 'OFC', 'V1L', 'V1R', 'HIPPOL', 'HIPPOR', 'VentStrL', 'VentStrR'};


% ====================== labelling by task variable =======================
LABEL = {'Goal(6,7,8)xUC LOROV goal', ...
    'Goal(6,7,8)xUC LOROV unc', ... 
    'Goal(6,7,8)xUC LOROV linear', ... 
    'Goal(6,7,8)xUC LOROV nonlinear'}; N_allclass = 6 * ones(1, length(LABEL)); SepClass = {' goal', ' unc', ' linear', ' nonlinear'};



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
% measure = 'diff';
% measure = 'ratio';
% measure = 'cosine';
% measure = 'corr_sig';   % significance of weight vector correlation
% measure = 'GI';

% multiple measures (mainly for 'mean_msr_box' only, not for 'indiv_msr_vec')
% measure_cell = {'CCGP', 'acc_mean', 'acc_mean', 'acc_mean'}; 
% measure_cell = {'CCGP', 'acc_mean', 'acc_mean'}; 
measure_cell = {};
% measure_cell = {'CCGP', 'acc_mean'}; % CCGP first, SD second for GI/ratio
% measure_cell = {'acc_mean', 'CCGP'}; % SD first, CCGP second for diff

% For shattering CCGP
CClabel_name = [];
% CClabel_name = 'CmplxCond';
% CClabel_name = 'UncCond';
% CClabel_name = 'Goal(6,7,8)';

% --- decoding ---
% measure = 'CCGP';
% measure = 'accs';
% measure = 'diff accs CCGP';
% measure = 'cosine';



% ========================= save setting =========================
% --- shattering ---
save_path = 'data/shattering_result';

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
% save_path = 'C:\Users\User\Desktop\Research\Dimensionality_control_for_prefrontal_meta_RL\temp_results\fitcecoc_decoding';

% save_name = ['fitrlinear_corrcoef_' Exp_sfx];
% save_name = ['fitrlinear_corrcoef_' Exp_sfx 'nullS3_'];

% save_name = ['fitcecoc_CVacc_' Exp_sfx];

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


% dichotomy map
diMap = [9, ... % 1) 1 4 vs 2 3 5 6 (goal 6) [2;3;5;6]
    18, ... % 2) 2 5 vs 1 3 4 6 (goal 7) [1;3;4;6]
    27, ... % 3) 3 6 vs 1 2 4 5 (goal 8) [3;6]
    7, ... % 4) 1 2 3 vs 4 5 6 (uncertainty) [4;5;6]
    1, ... % 5) 1 vs 2 3 4 5 6 (linear) [2;3;4;5;6]
    2, ... % 6) 2 vs 1 3 4 5 6 (linear) [1;3;4;5;6]
    4, ... % 7) 3 vs 1 2 4 5 6 (linear) [1;2;4;5;6]
    8, ... % 8) 4 vs 1 2 3 5 6 (linear) [1;2;3;5;6]
    16, ... % 9) 5 vs 1 2 3 4 6 (linear) [1;2;3;4;6]
    31, ... % 10) 6 vs 1 2 3 4 5 (linear) 6
    3, ... % 11) 1 2 vs 3 4 5 6 (linear) [3;4;5;6]
    5, ... % 12) 1 3 vs 2 4 5 6 (linear) [2;4;5;6]
    6, ... % 13) 2 3 vs 1 4 5 6 (linear) [1;4;5;6]
    24, ... % 14) 4 5 vs 1 2 3 6 (linear) [1;2;3;6]
    23, ... % 15) 4 6 vs 1 2 3 5 (linear) [4;6]
    15, ... % 16) 5 6 vs 1 2 3 4 (linear) [5;6]
    11, ... % 17) 1 2 4 vs 3 5 6 (linear) [3;5;6]
    19, ... % 18) 1 2 5 vs 3 4 6 (linear) [3;4;6]
    13, ... % 19) 1 3 4 vs 2 5 6 (linear) [2;5;6]
    26, ... % 20) 1 3 6 vs 2 4 5 (linear) [1;3;6]
    25, ... % 21) 1 4 5 vs 2 3 6 (linear) [2;3;6]
    22, ... % 22) 1 4 6 vs 2 3 5 (linear) [1;4;6]
    17, ... % 23) 1 5 vs 2 3 4 6 (nonlinear) [2;3;4;6]
    30, ... % 24) 1 6 vs 2 3 4 5 (nonlinear) [1;6]
    10, ... % 25) 2 4 vs 1 3 5 6 (nonlinear) [1;3;5;6]
    29, ... % 26) 2 6 vs 1 3 4 5 (nonlinear) [2;6]
    12, ... % 27) 3 4 vs 1 2 5 6 (nonlinear) [1;2;5;6]
    20, ... % 28) 3 5 vs 1 2 4 6 (nonlinear) [1;2;4;6]
    28, ... % 29) 1 2 6 vs 3 4 5 (nonlinear) [1;2;6]
    21, ... % 30) 1 3 5 vs 2 4 6 (nonlinear) [2;4;6]
    14]; % 31) 1 5 6 vs 2 3 4 (nonlinear) [1;5;6]


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

                    
                    try

                        if isempty(SepClass)
                            data_load = load([save_path '/' roi_type '/' ROI{roii} '/' , ...
                                LABEL{labi} '/' num2str(ev_name) '/', ...
                                save_name num2str(id) '.mat']);
                        else
                            data_load = load([save_path '/' roi_type '/' ROI{roii} '/' , ...
                                erase(LABEL{labi}, num2str(SepClass{labi})) '/' num2str(ev_name) '/', ...
                                save_name num2str(id) '.mat']);
                        end

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

                                if contains(LABEL{labi}, 'Goal(6,7,8) x Unc') || contains(LABEL{labi}, 'Goal(6,7,8)xUC')

                                    % goal main effect: 9 8 27
                                    % UC main effect: 7
                                    % nonlinear interaction: 10 12 17 20 29 30 (2 vs 4) 14 21 28 (3 vs 3)
                                    % full dichotomies: goal(3) + UC(1) + linear_it(18) + nonlinear(9)

                                    vec_len_temp = vec_len;

                                    switch SepClass{labi}
                                        case ' main'
                                            selec_dich = [9 18 27 7];
                                        case ' goal'
                                            selec_dich = [9 18 27]; % including Goal(6,7,8) dichotomies only
                                        case ' unc'
                                            selec_dich = 7; % including the UC_G dichotomy only
                                        case ' linear'
                                            selec_dich = 1:vec_len;
                                            selec_dich([9 18 27 7 10 12 17 20 29 30 14 21 28]) = []; % including linear 'interactions' only
                                        case ' nonlinear'
                                            selec_dich = [10 12 17 20 29 30 14 21 28]; % including nonlinear interactions(dichotomies) only
                                        otherwise
                                            selec_dich = str2double(SepClass{labi});
                                    end

                                    vec_len_temp = length(selec_dich);
                                    msr_vec = squeeze(mean(msr_vec(:, selec_dich, :, :), [3 4])); % (msr_vec: nchoosek(c,m), (2^m - 2)/2)

                                end

                                mean_msr_box(id, labi, roii, evi) = mean(msr_vec);
                                indiv_msr_vec(id, labi, roii, evi, 1:vec_len_temp) = msr_vec;

                            else
                                % ===== load shattering measure =====
                                % msr_mat: acc. or sep. of all possible classification
                                % size: (N_allclass choose N_shatteredclass) ...
                                % x (2^N_shatteredclass - 2)
                                msr_vec = data_load.save_info.acc_box. ...
                                    ShatteredClassesNumber{N_allclass(labi)}.(measure);

                                % ***** original *****
                                if isempty(SepClass)
                                    mean_msr_box(id, labi, roii, evi) = mean(msr_vec);

                                else

                                    if contains(LABEL{labi}, 'Goal(6,7,8) x Unc') || contains(LABEL{labi}, 'Goal(6,7,8)xUC')

                                        % goal main effect: 9 8 27
                                        % UC main effect: 7
                                        % nonlinear interaction: 10 12 17 20 29 30 (2 vs 4) 14 21 28 (3 vs 3)
                                        % full dichotomies: goal(3) + UC(1) + linear_it(18) + nonlinear(9)

                                        switch SepClass{labi}
                                            case ' main'
                                                selec_dich = [9 18 27 7];
                                            case ' goal'
                                                selec_dich = diMap(1:3);
                                            case ' unc'
                                                selec_dich = diMap(4);
                                            case ' linear'
                                                selec_dich = diMap(5:22);
                                                % selec_dich = 1:vec_len;
                                                % selec_dich([9 18 27 7 10 12 17 20 29 30 14 21 28]) = []; % including linear 'interactions' only
                                                % %                                                 case 'linear'
                                                % %                                                     selec_dich = 1:vec_len;
                                                % %                                                     selec_dich([10 12 17 20 29 30 14 21 28]) = []; % including linear dichotomies only
                                            case ' linear1'
                                                selec_dich = diMap(5:10);
                                            case ' linear2'
                                                selec_dich = diMap(11:16);
                                            case ' linear3'
                                                selec_dich = diMap(17:22);
                                            case ' nonlinear'
                                                selec_dich = diMap(23:31);
                                            case ' nonlinear1'
                                                selec_dich = diMap(23:28);
                                            case ' nonlinear2'
                                                selec_dich = diMap(29:31);
                                            otherwise
                                                selec_dich = str2double(SepClass{labi});
                                        end

                                        mean_msr_box(id, labi, roii, evi) = mean(msr_vec(selec_dich));

                                    end

                                end

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
% disp(data_load.run_info.fitcecoc_setting)
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
    roii5 = [8 9]; % preSMA
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
vis_mode = 'regional_info';
% vis_mode = 'temporal_info'; temp_mode = 'group';     set(0,'DefaultAxesFontSize',12,'DefaultTextFontSize',12);
% vis_mode = 'temporal_info'; temp_mode = 'individual';     set(0,'DefaultAxesFontSize',12,'DefaultTextFontSize',12);
% vis_mode = 'corr_behav';  corr_mode = 'original';  set(0,'DefaultAxesFontSize',10,'DefaultTextFontSize',10);
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
% BHV_load = 1;
BHV_load = 0;
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
% vLABELi = 2;
% vLABELi = 3;
% vLABELi = 4;
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
    % y_lower = chance_lvl - 0.025;
    y_lower = chance_lvl - 0.01;
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
    y_upper = chance_lvl + 0.055;
    % y_upper = chance_lvl + 0.06;
%     y_upper = chance_lvl + 0.075;
    % y_upper = chance_lvl + 0.08;
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
    case 'regional_info'



        % ===== Group result plot =====

        nSbj = length(vIDi);
        nLab = length(vLABELi);
        nROI = length(vROIi);

        msrmat = squeeze(mean(mean_msr_box(vIDi, vLABELi, vROIi, vEVENTi), 4)); % (nSbj, nLab, nROI)
        pMat = nan(nLab, nROI);
        for li = 1:nLab
            for ri = 1:nROI
                msrVec = squeeze(msrmat(:, li, ri));
                [~, p] = ttest(msrVec - chance_lvl(li));
                pMat(li, ri) = p;
            end
        end
        

        
        % averageing event
        msrmat_sbjmean = squeeze(mean(msrmat, 1)); % LEN_LABEL x LEN_ROI
        msrmat_sbjste = squeeze(std(msrmat, 0, 1) ... 
            ./sqrt(sum(~isnan(msrmat), 1))); % LEN_LABEL x LEN_ROI
        





        % #################################################################
        % Single label vs chance level
        % #################################################################

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

                % #####################################################
                % Permutation test for label vs. chance
                % #####################################################
                % We do a sign-flipping test. Under the null (accuracy == chance),
                % the difference (msrmat(:,vroii) - chance_lvl) could be +/- with equal prob.

                nPerm = 5000; % # of permutations
                realDiff = mean(msrmat(:, vroii) - chance_lvl(vLABELi));
                permDist = zeros(nPerm, 1);

                for perm_i = 1:nPerm
                    rng(perm_i)
                    % randomly flip sign of difference
                    flipSign = randsample([1, -1], length(msrmat(:, vroii)), true);
                    tmp = (msrmat(:, vroii) - chance_lvl(vLABELi)) .* flipSign;
                    permDist(perm_i) = mean(tmp);
                end

                % two-tailed p-value
                p_perm = mean(abs(permDist) >= abs(realDiff));
                % #####################################################
                % End of added permutation code
                % #####################################################

                if h
                    sig_roi(vroii) = h; 
                    fprintf(ROI{vROIi(vroii)})
                    fprintf(' p=%.2e \n', p)
                    disp(stats)
                end
                if p < 0.001
                    sig_text{vroii} = '***';
                elseif p < 0.01
                    sig_text{vroii} = '**';
                elseif p < 0.05
                    sig_text{vroii} = '*';
                else
                    sig_text{vroii} = ''; % no mark
                end

                % ----------------------------------------------------
                % Visualization of permutation distribution
                % ----------------------------------------------------
                figure('Name', sprintf('Permutation Dist: %s (ROI:%s)', LABEL{vLABELi}, ROI{vROIi(vroii)}))
                histogram(permDist, 30, 'FaceColor', [0.6 0.6 0.9], 'EdgeColor', 'none');
                hold on
                line([realDiff, realDiff], ylim, 'Color','r','LineWidth',2)
                xlabel('Mean Difference (msrmat - chance)')
                ylabel('Frequency')
                title({[ROI{vROIi(vroii)}, ' label vs chance: Perm Dist'], ...
                    ['Observed = ', num2str(realDiff, '%.3f'), ...
                    ', p_{perm} = ', num2str(p_perm, '%.4f')]})
                hold off
                % ----------------------------------------------------

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
        % #################################################################

        




        % #################################################################
        % Main visualization
        % #################################################################

        figure
        % subplot(2, length(vROIi), 1:length(vROIi))
        B = bar(msrmat_sbjmean', 'EdgeColor', 'none'); axis square; box off; % (nLab, nROI)
        % 4 colors
        diColor1(1,:) = [77 161 169]/255; % goal
        diColor1(2,:) = [226 16 103]/255; % uncertainty
        diColor1(3,:) = [155 80 148]/255; % linear
        diColor1(4,:) = [94 35 157]/255;  % nonlinear
        % 7 colors
        diColor2(1,:) = [77 161 169]/255;       % goal
        diColor2(2,:) = [226 16 103]/255;       % uncertainty
        diColor2(3,:) = 1.2*[155 80 148]/255;   % linear1
        diColor2(4,:) = .6*[226 16 103]/255;    % linear2 (close to uncertainty)
        diColor2(5,:) = .6*[77 161 169]/255;    % linear3 (close to goal)
        diColor2(6,:) = [94 35 157]/255;        % nonlinear 1
        diColor2(7,:) = [94 35 157]/255;        % nonlinear 2
        if length(vLABELi) == 4
            diff_erbr = .275;
            B(1).FaceColor = diColor1(1,:);
            B(2).FaceColor = diColor1(2,:);
            B(3).FaceColor = diColor1(3,:);
            B(4).FaceColor = diColor1(4,:);
        elseif length(vLABELi) == 7
            diff_erbr = .35;
            B(1).FaceColor = diColor2(1,:);
            B(2).FaceColor = diColor2(2,:);
            B(3).FaceColor = diColor2(3,:);
            B(4).FaceColor = diColor2(4,:);
            B(5).FaceColor = diColor2(5,:);
            B(6).FaceColor = diColor2(6,:);
            B(7).FaceColor = diColor2(7,:);
        end

        hold on
                

        % Make non-significant result gray
        for li = 1:nLab
            B(li).FaceColor = 'flat';

            for ri = 1:nROI
                if pMat(li, ri) >= .05                    
                    B(li).CData(ri,:) = [.8 .8 .8];
                else
                    if nLab == 4
                        B(li).CData(ri,:) = diColor1(li,:);
                    elseif nLab == 7
                        B(li).CData(ri,:) = diColor2(li,:);
                    end
                end

            end
        end

        
        % Individual dot plot, on the error bar plot
        for roii = 1:length(vROIi)
            if length(vLABELi) == 2
                plot([roii - diff_erbr, roii + diff_erbr], ...
                    squeeze(msrmat(:,:,roii))', ... % msrmat: [N_sbj x N_label x N_roi]
                    'LineWidth', 0.01, 'Marker', 'o', 'Color', [.8 .8 .8]);
            elseif length(vLABELi) == 1
                plot(roii, ...
                    squeeze(msrmat(:, roii))', ... % msrmat: [N_sbj x N_roi]
                    'LineWidth', 0.01, 'Marker', 'o', 'Color', [.8 .8 .8]);
            end
        end

        xdiff_erbr = linspace(-diff_erbr, diff_erbr, length(vLABELi));
        xtick_erbr = reshape((xdiff_erbr)' + (1:length(vROIi)), 1, []);     % LAB x REG -- 1 * (LAB * REG)
        if strcmp(measure, 'accs_rr')
            E = errorbar(xtick_erbr, reshape(msrmat_sbjmean,1,[]), ...
                reshape(msrmat_sbjste,1,[]), 'LineStyle', 'none', 'LineWidth', 1.5);
        else
            E = errorbar(xtick_erbr, reshape(msrmat_sbjmean,1,[]), ...
                reshape(msrmat_sbjste,1,[]), 'LineStyle', 'none', 'LineWidth', 3);
        end
        % errorbar LineWidth: original - 1.2, bold (BIG BRAIN 3rd seminar & COSYNE 2023 fig 3 bar) - 3
        E.Color = colors(21, :);

        % Draw the chance level line
        line([0 length(vROIi)+1], [chance_lvl(1) chance_lvl(1)], ...
            'Color', [.5 .5 .5], ...
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
%         ylim([0.35 0.7]); % To see UC G vs UC H individual plot
%         ylabel(['Mean ' measure], 'Interpreter', 'none')
        ylabel(msr_name, 'Interpreter', 'none')
        try
            legend(B(:), LABEL{vLABELi})
        catch
            legend(LABEL{vLABELi})
        end









        % ===== significance check B =====
        if length(vLABELi) > 1

            nSbjv = length(vIDi);
            nROIv = length(vROIi);
            nLABv = length(vLABELi);
            lab_pairs = nchoosek(1:nLABv, 2);
            nCk = size(lab_pairs, 1);

            diffArr = nan(nROIv, nSbjv, nCk);
            sigPairs = nan(nROIv, nCk);

            % To visualize pairwise comparison significance
            % figure
            % figure('Name','Pairwise Permutation Results - Ver1')
            % nColors = nCk;
            % pairColors = jet(nColors); % or lines(nColors), or any colormap
            pairColors = cat(1, .8 * ones(3, 3), repmat([94 35 157]/255, 3, 1));
            pairColors = pairColors([1 2 4 3 5 6], :);
            % figure('Name','Pairwise Permutation Results - Ver2')

            for vroii = 1:nROIv

                lab_pair_sig = ones(nLABv, nLABv);
                msrmat_roi = squeeze(msrmat(:, :, vroii));      % nSbjv x nLABv
                %                 pvals = nan * ones(size(lab_pairs, 1), 1);
                pvals = cell(nCk, 1);

                % subplot(ceil(nROIv/2), 2, vroii) % or subplot(1,nROIv,vroii) if you prefer 1 row
                hold on
                legendEntries = cell(nCk,1);

                disp(ROI{vROIi(vroii)})
                for pairi = 1:nCk
                    lab_pair = lab_pairs(pairi, :);
                    msrVec1 = msrmat_roi(:, lab_pair(1));
                    msrVec2 = msrmat_roi(:, lab_pair(2));

                    % Pairwise t-test (original)
                    [h, p, ~, stats] = ttest(msrVec1, msrVec2, 'Alpha', sgnf_lvl);

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
                    fprintf([erase(SepClass{vLABELi(lab_pair(1))}, ' ') ' vs' SepClass{vLABELi(lab_pair(2))} ...
                        ': t(%d)=%.3f\t p=%.3e (' pvals{pairi} ')\n'], stats.df, stats.tstat, p)
                    
                    
                end
                disp(' ')

            end




        end   % if length(vLABELi) > 1



















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
        % visDiff = {1:3, 4};
        % 
        % % --- preallocate effect-size distribution ---
        % if ~isempty(visDiff)
        %     nLabv = length(visDiff{1}) * length(visDiff{2});
        %     vLABELi = 1:nLabv;
        % end
        % meanBoot = nan(numBoot, nLabv, nROIv);
        % sigBoot = nan(numBoot, nLabv, nROIv);
        % pairSig = false(numBoot, 3, nROIv);
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
        %     % ---------------------------------------------
        %     % Pairwise t-tests: label i (1..3) vs label 4
        %     % ---------------------------------------------
        %     for roii = 1:nROIv
        %         for i = 1:3  % compare label i to label 4
        %             msrVec_i = bootData(:, i,   roii);  % subsetSize x 1
        %             msrVec_4 = bootData(:, 4,   roii);  % subsetSize x 1
        %             [h, pVal] = ttest(msrVec_i, msrVec_4, 'Alpha', sgnf_lvl);
        %             pairSig(b, i, roii) = h;  % 1 if significant
        %         end
        %     end
        % 
        %     % --- then your accuracy-chance for each label × ROI ---
        %     meanBoot(b,:,:) = bootDiff;      % store into a 3-D array now
        % 
        %     % % --- then your Cohen's d for each label × ROI ---
        %     % effDist(b,:,:) = bootDiff ./ bootStd;      % store into a 3-D array now
        % end
        % 
        % % --- modified: statistical test per label and ROI ---
        % nLabv = length(vLABELi);
        % nROIv  = length(vROIi);
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
        %     data_bs = squeeze(meanBoot(:, :, roii)); % [numBoot × nLabv]
        %     sig_bs = logical(squeeze(sigBoot(:, :, roii)));
        % 
        %     nSig = sum(sig_bs, 1); % 1 x nLabv
        %     beeX1 = []; beeY1 = []; % for significant result
        %     beeX2 = []; beeY2 = []; % for nonsignificant result
        %     if ~isempty(visDiff)
        %         colors_temp = [77 161 169; 226 16 103; 155 80 148] / 255;
        %     else
        %         colors_temp = [77 161 169; 226 16 103; 155 80 148; 94 35 157] / 255;
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
        %         'colormap', repmat([.5 .5 .5], sum((numBoot-nSig)>0), 1), ...
        %         'overlay_style', false, ...
        %         'corral_style', 'random');
        % 
        %     % all results
        %     if subsetSize == 19 && numBoot == 20
        %         beeswarm(reshape(repmat(1:nLabv, numBoot, 1), [], 1), ...
        %             reshape(data_bs, [], 1), ...
        %             'overlay_style', 'ci', ...
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
        %     % % ---------------------------------------------
        %     % % Visualize pairwise significance with lines
        %     % % ---------------------------------------------
        %     % % Recall: pairSig(b, i, roii) => 1 if label i vs label 4 is significant in iteration b
        %     % % We'll loop over each iteration b, and for each i=1..3 that is significant,
        %     % % draw a gray line from the (x=i, y=data_bs(b,i)) to (x=4, y=data_bs(b,4)).
        %     % for b = 1:numBoot
        %     %     for i = 1:3
        %     %         if pairSig(b, i, roii) == 1
        %     %             % Draw a line
        %     %             xCoord = [i, 4];
        %     %             yCoord = [data_bs(b, i), data_bs(b, 4)];
        %     %             plot(xCoord, yCoord, ':', 'Color', [0.8 0.8 0.8], 'LineWidth', .3);
        %     %         end
        %     %     end
        %     % end
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
        %         ylabel('SD difference');
        %         % ylim([-.01 .01])
        %         xticklabels({'goal-nonlin.', 'unc.-nonlin.', 'linear-nonlin.'})
        %         ylim([-.001 .035])
        %     else
        %         ylabel('SD - chance');
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
        %         xticklabels(SepClass(vLABELi))
        %     end
        % 
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
disp('END:shattering_result')

