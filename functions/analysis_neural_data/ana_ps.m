
function [] = ana_ps(id, varargin)

% comment out or set your own path
% addpath(genpath('/home/ydsung/A_Research/princeton-mvpa-toolbox-master'));

%%
% ========================= Experiment =========================
Exp = 'Lee2014'; Exp_sfx = []; prefix = 'wra';
% Exp = 'Lee2014pp2'; Exp_sfx = '_Lee2014pp2'; prefix = 'wra';
% Exp = 'Kim2019'; Exp_sfx = '_Kim2019_wa'; prefix = 'wa';
% Exp = 'Heo2021'; Exp_sfx = '_Heo2021_wra'; prefix = 'wra';
% Exp = 'Heo2021'; Exp_sfx = '_Heo2021_swra'; prefix = 'swra';
% Exp = 'Shin2025'; Exp_sfx = '_Shin2025_swar'; prefix = 'swar';


% ================ ROI mask for multi-voxel pattern (MVP) =================
% main ROI set
% ROI = {'F3TL', 'F3TR', 'DLPFCL', 'DLPFCR', 'OFC', 'V1L', 'V1R', 'HIPPOL', 'HIPPOR', 'VentStrL', 'VentStrR'};

% additional ROIs
% ROI = {'ACCL', 'ACCR', 'SMAL', 'SMAR'};
% ROI = {'ACCL', 'ACCR', 'preSMAL', 'preSMAR'};

% full ROI set
ROI = {'F3TL', 'F3TR', 'DLPFCL', 'DLPFCR', 'OFC', 'ACCL', 'ACCR', 'preSMAL', 'preSMAR', 'V1L', 'V1R', 'HIPPOL', 'HIPPOR', 'VentStrL', 'VentStrR'};

% AAL only
% ROI = {'DLPFCL', 'DLPFCR', 'OFC'};
% ROI = {'OFC'};

% test
% ROI = {'HIPPOL', 'HIPPOR'};



% ====================== labelling by task variable =======================

% LABEL_NULL = {'S3'};
LABEL_NULL = {};


% LABEL = {'UncCond'}; N_allclass = 2; trialSelect = {[]}; lab_shuffle = 0; % for fig.2c & d
% LABEL = {'UncCond', 'UncCond'}; N_allclass = [2 2]; trialSelect = {'G', 'H'}; lab_shuffle = 0; % for fig.2c & d
% LABEL = {'UncCond', 'UncCond'}; N_allclass = [2 2]; trialSelect = {'cL', 'cH'}; lab_shuffle = 0;
% LABEL = {'UncCond', 'UncCond'}; N_allclass = [2 2]; trialSelect = {'binPMB1', 'binPMB0'}; lab_shuffle = 0;
% LABEL = {'UncCondd1', 'UncCondd2', 'UncCondd1', 'UncCondd2', ...
%     'UncConda1', 'UncConda2', 'UncConda1', 'UncConda2'}; N_allclass = 2 * ones(size(LABEL)); trialSelect = {'G', 'G', 'H', 'H', 'G', 'G', 'H', 'H'}; lab_shuffle = 0;
% LABEL = {'UncCond delay1', 'UncCond delay2', 'UncCond delay1', 'UncCond delay2', ...
%     'UncCond advance1', 'UncCond advance2', 'UncCond advance1', 'UncCond advance2'}; N_allclass = 2 * ones(size(LABEL)); trialSelect = {'G', 'G', 'H', 'H', 'G', 'G', 'H', 'H'}; lab_shuffle = 0;
% LABEL = {'UncCond delay1', 'UncCond delay2', 'UncCond delay1', 'UncCond delay2'}; N_allclass = [2 2 2 2]; trialSelect = {'G', 'G', 'H', 'H'}; lab_shuffle = 0;
% LABEL = {'UncCond advance1', 'UncCond advance2', 'UncCond advance1', 'UncCond advance2'}; N_allclass = [2 2 2 2]; trialSelect = {'G', 'G', 'H', 'H'}; lab_shuffle = 0;
% LABEL = {'UncCond delay1', 'UncCond delay2', 'UncCond delay1', 'UncCond delay2'}; N_allclass = [2 2 2 2]; trialSelect = {'cL', 'cL', 'cH', 'cH'}; lab_shuffle = 0;
% LABEL = {'UncCond shuff', 'UncCond shuff'}; N_allclass = [2 2]; trialSelect = {'G', 'H'}; lab_shuffle = 1;
% LABEL = {'UncCond perm', 'UncCond perm'}; N_allclass = [2 2]; trialSelect = {'G', 'H'}; lab_shuffle = 0;
LABEL = {'Goal(6,7,8)', 'S2', 'S3'}; N_allclass = [3 4 4]; trialSelect = {[], 'G', 'G'}; lab_shuffle = 0; % for fig.2a & b


CClabel_name = 'UncCond'; ccSelect = [];



% to categorize continuous variables into several classes
categorize = cell(1, length(LABEL));  % empty cell for no categorization
% categorize = num2cell(2 * ones(1, length(LABEL)));
% categorize = num2cell(4 * ones(1, length(LABEL)));
% categorize = num2cell(5 * ones(1, length(LABEL)));



% ========================= MVP setting =========================
% prefix = 'wra*';  % preprocessing prefix
% prefix = 'swa*';

% SignalType = 'percent';
SignalType = 'zscore'; % default setting

% ROI_types = {};
% (main ROI set)
% ROI_types = {'AAL3', 'AAL3', 'AAL', 'AAL', 'AAL', 'AAL3', 'AAL3', 'AAL3', 'AAL3', 'AAL3', 'AAL3'}; % +DLPFC & R/L integrated OFC
% (full ROI set)
ROI_types = {'AAL3', 'AAL3', 'AAL', 'AAL', 'AAL', 'AAL3', 'AAL3', 'JuBrain', 'JuBrain', 'AAL3', 'AAL3', 'AAL3', 'AAL3', 'AAL3', 'AAL3'};
% (additional ROIs)
% ROI_types = {'AAL3', 'AAL3', 'JuBrain', 'JuBrain'};

% roi_type = 'AAL3';
% roi_type = 'AAL';
% roi_type = 'JuBrain';
roi_type = [];

MaskNum_map = containers.Map({'F3TL', 'F3TR', 'V1L', 'V1R', 'HIPPOL', 'HIPPOR', 'VentStrL', 'VentStrR', ...
    'DLPFCL', 'DLPFCR', 'OFC', 'ACCL', 'ACCR', 'SMAL', 'SMAR', 'preSMAL', 'preSMAR'}, ...
    {9, ... % F3TL
    10, ... % F3TR
    47, ... % V1L
    48, ... % V1R
    41, ... % HIPPOL
    42, ... % HIPPOR
    157, ... % VentStrL
    158, ... % VentStrR
    7, ... % DLPFCL
    8, ... % DLPFCR
    [5 9 15 27 6 10 16 28], ... % OFC
    [153 155], ... % ACCL
    [154 156], ... % ACCR
    15, ... % SMAL
    16, ... % SMAR
    nan, ... % preSMAL
    nan}); % preSMAR

% EVENT = {[1 2], [4 5], [7 8]};
EVENT = [1 2 4 5 7 8 9 10];      % mainly for 'strict' MVP
% EVENT = [9 10];      % mainly for 'strict' MVP
% EVENT = 4;      % 2021-07-07 test for weight map
% EVENT = [1 2 5 7 8 9 10];      % 2021-07-07 test for weight map
% EVENT = -3:8;
% EVENT = -3:12;
% EVENT = 1:8;
% EVENT = {[2 3], [5 6], 8};

strict = 0;
% strict = 1;
interpolation = 0;



% ========================= Undersampling setting =========================
% Repeat = 1;    % # of undersample
% Repeat = 5;    % # of undersample
% Repeat = 10;    % # of undersample
% Repeat = 20;    % # of undersample
Repeat = 100;    % # of undersample

rRepeat = 0;
% rRepeat = 5;
% rRepeat = 10;
% rRepeat = 100;


% ========================= save setting =========================
save_path = 'ps_result';
save_name = ['ps_' CClabel_name ccSelect Exp_sfx];






% main ====================================================================


fprintf('%d ================================================ \n', id)


% #########################################################################
% =========================== loop1: label load ===========================
% #########################################################################


labels = [];

for labi = 1:length(LABEL)
    
    if contains(LABEL{labi}, 'delay')
        var_name = split(LABEL{labi}, ' delay');
        var_name = var_name{1};
    elseif contains(LABEL{labi}, 'advance')
        var_name = split(LABEL{labi}, ' advance');
        var_name = var_name{1};
    else
        var_name = LABEL{labi};
    end    
    var_name = erase(var_name, ' shuff');
    var_name = erase(var_name, ' LOROV');
    var_name = erase(var_name, ' 5fold');
    var_name = erase(var_name, ' perm');
    
    if contains(var_name, ' run')
        runs_temp = var_name((strfind(var_name, ' run')+3):end);    % 'run keyword must be the last
        runs = nan(1,length(runs_temp));
        for i=1:length(runs); runs(i) = str2double(runs_temp(i)); end
        session = arbMBMF_load_var(Exp, 'Session', id, []);
        selector = 1 * ismember(session, runs);
        var_name = erase(var_name, var_name(strfind(var_name, ' run'):end));
    elseif ~isempty(trialSelect{labi})
        selector = trialSelect{labi};
        trial_select_sfx = [' ' trialSelect{labi}];
        trial_select_sfx_cell{labi} = trial_select_sfx;
    else
        selector = [];
        trial_select_sfx = [];
        trial_select_sfx_cell{labi} = trial_select_sfx;
    end
   
    label = arbMBMF_load_var(Exp, var_name, id, selector);

    if contains(LABEL{labi}, ' delay')
        tokens = regexp(LABEL{labi}, 'delay\s*(\d+)', 'tokens');
        nDelay = str2double(tokens{1}{1}); % # of trials to delay the label
        % labelOri = label; % only used for test
        label = [nan(1, nDelay), label( 1:(end-nDelay) )];
        % imagesc([labelOri; label]); % only used for test
    end

    if contains(LABEL{labi}, ' advance')
        tokens = regexp(LABEL{labi}, 'advance\s*(\d+)', 'tokens');
        nAdvance = str2double(tokens{1}{1}); % # of trials to advance the label
        % labelOri = label; % only used for test
        label = [label( (1+nAdvance):end ), nan(1, nAdvance)];
        % imagesc([labelOri; label]); % only used for test
    end
    
    if lab_shuffle || contains(LABEL{labi}, 'shuff')
        valid_lab_idx = ~isnan(label);
        valid_lab_shuffled = label(valid_lab_idx);
        valid_lab_shuffled = valid_lab_shuffled(randperm(length(valid_lab_shuffled)));
        label(valid_lab_idx) = valid_lab_shuffled;
    end
    
    if ~isempty(categorize{labi})
        label = regr2class(label, categorize{labi});
    end

    labels = cat(1, labels, label);

end % labi

% load label to be nullified
if ~isempty(LABEL_NULL)
    label_nulls = [];
    for ni = 1:length(LABEL_NULL)
        %     label_null = arbMBMF_load_var(Exp, var_name, id, selector);
        label_null = arbMBMF_load_var(Exp, LABEL_NULL{ni}, id, []);
        label_nulls = [label_nulls; label_null];
    end
    save_name = [save_name, 'null', LABEL_NULL{:} '_'];
end


% For CCGP: cross-condition labeling for each trial
if ~isempty(CClabel_name)
    CClabel = arbMBMF_load_var(Exp, CClabel_name, id, ccSelect);
else 
    CClabel = [];
end


% #########################################################################
% ============== loop2: decoding label for each ROI & event ===============
% #########################################################################


for roii = 1:length(ROI)
    roi = ROI{roii};
    fprintf([roi ' ================================ \n'])
    
    if ~isempty(ROI_types)
        roi_type = ROI_types{roii};
    end


    % multi-voxel pattern MVP: [N_voxel x N_event x N_trial]
    MVPfull = arbMBMF_boldpat(Exp, roi, id, ...
        'SignalType', SignalType, ...
        'roi_type', roi_type, ...
        'MaskNum', MaskNum_map(roi), ...
        'Event', EVENT, ...
        'strict', strict, ...
        'interpolation', interpolation, ...
        'prefix', prefix);


    % removing mean neural acvitity evoked by null-label
    if ~isempty(LABEL_NULL)
        for ni = 1:length(LABEL_NULL)
            MVPfull = lab_pat_nullify(label_nulls(ni, :), MVPfull, ...
                'NullDim', 3);
        end
    end
    
    for labi = 1:length(LABEL)
        
        fprintf([LABEL{labi} ' ================ \n'])
        
        for evi = 1:length(EVENT)
            
            % loaded label 2021-03-21 - under the event loop for some reason
            label = labels(labi, :);
            
            % event specific MVP
            if ndims(MVPfull) > 2
                MVP = squeeze(MVPfull(:,evi,:));
            else
                MVP = MVPfull;
            end
            % event name for saving
            if iscell(EVENT)
                ev_name = EVENT{evi};
            else
                ev_name = EVENT(evi);
            end
                        
            % sanity check
            trial_len = size(label, 2);
            if size(MVP,2)~=trial_len
                disp('label-pattern size unmatch'); 
                continue
            end
            
            if ~contains(save_name, 'fitrlinear')
                unqSet = unq_elms(label);
                if length(unqSet)~=N_allclass(labi)
                    disp('(arbMBMF_decoding) label element error');
                    continue
                end
            end
            
            % running information
            run_info.Exp = Exp;
            run_info.ROI = ROI;
            run_info.LABEL = LABEL;
            run_info.N_Allclass = N_allclass;
            run_info.categorize = categorize;
            run_info.lab_shuffle = lab_shuffle;
            run_info.MVP_setting.SignalType = SignalType;
            run_info.MVP_setting.roi_type = roi_type;
            run_info.MVP_setting.EVENT = EVENT;
            run_info.MVP_setting.strict = strict;
            run_info.MVP_setting.interpolation = interpolation;
            run_info.MVP_setting.prefix = prefix;
            run_info.stat_setting.Repeat = Repeat;
            run_info.stat_setting.rRepeat = rRepeat;
            
            % save information
            save_info.id = id;
            save_info.roi = ROI{roii};
            save_info.roi_type = roi_type;
            save_info.label_name = LABEL{labi};
            save_info.event = ev_name;
            save_info.save_path = save_path;
            save_info.save_name = save_name;

            
            % ================ main coding distance analysis ================
            MVP_ori   = MVP;      % original MVP timecourse
            label_ori = label;    % original label timecourse
            
            runs = arbMBMF_load_var(Exp, 'Session', id, []);
            [setRun, lenRuns] = unq_elms(runs);
            [setCC, lenCC] = unq_elms(CClabel);

            nRun = length(setRun);
            nCC = length(setCC);
            nVox = size(MVP_ori, 1);
            nClass = N_allclass(labi);

            % class pairs to calculate coding directions
            cPairs = nchoosek(1:nClass, 2);
            nPair = size(cPairs, 1);

            % Prepare to store parallelism score results
            ps  = nan(nPair, Repeat);
            
            for rpi = 1:Repeat

                rng(rpi); % Keep the same randomization approach if needed

                ps_run = nan(nPair, nRun);

                for runi = 1:nRun
                    thisRun = setRun(runi);  % e.g., runSet might be [1, 2, 3, ...]

                    cDirec = nan(nCC, nPair, nVox); % coding direction

                    for cci = 1:nCC
                        thisCC = setCC(cci);

                        % idx extraction
                        tempIdx = (runs == thisRun) & (CClabel == thisCC);

                        % === Extract data for this run ===
                        [label_run, MVP_run] = lab_pat_undersample( ...
                            label_ori(tempIdx), MVP_ori(:, tempIdx) ...
                            );
                        
                        % === Compute mean pattern per class ===
                        [~, ~, idxClass] = unique(label_run);

                        meanPatterns = zeros(nVox, nClass);
                        for c = 1:nClass
                            meanPatterns(:, c) = mean(MVP_run(:, idxClass == c), 2);  % (nVoxel, nClass)
                        end                      

                        for pi = 1:nPair
                            pair = cPairs(pi, :);
                            pat1 = meanPatterns(:, pair(1));
                            pat2 = meanPatterns(:, pair(2));
                            cDirec(cci, pi, :) = pat1 - pat2;
                        end

                    end

                    ccPairs = nchoosek(1:nCC, 2); 
                    nccPair = size(ccPairs, 1); % if nCC==2, then nccPair = 1
                    ps_cc = nan(nPair, nccPair);
                    for pi = 1:nPair
                        for cpi = 1:nccPair
                            ccPair = ccPairs(cpi, :);
                            dir1 = squeeze(cDirec(ccPair(1), pi, :)); % (nVox, 1)
                            dir2 = squeeze(cDirec(ccPair(2), pi, :)); % (nVox, 1)
                            nd1 = dir1 / norm(dir1); % normalized column vector
                            nd2 = dir2 / norm(dir2);
                            cs = nd1' * nd2; % cosine similarity
                            ps_cc(pi, cpi) = cs;
                        end
                    end

                    ps_run(:, runi) = mean(ps_cc, 2); % RHS: (nPair, 1), ccPairs-averaged

                end % for runi (loop over runs)

                % Weighted average of distances across runs
                noNanRun = all(~isnan(ps_run), 1);
                ps(:, rpi) = ps_run(:, noNanRun) * lenRuns(noNanRun) ...
                    / sum(lenRuns(noNanRun)); % RHS: (nPair, nRun) * (nRun, 1) / scalar

            end % for r = 1:Repeat

            % === Save info in the same structure you used before ===
            save_info.ps = ps; % (nPair, Repeat)
            save_info.psMean = mean(ps(:));  


            % save
            save_dir = [save_path '/' roi_type '/' ROI{roii} '/' , ...
                LABEL{labi} trial_select_sfx_cell{labi} '/' num2str(ev_name) '/'];

            try
                save([save_dir, ...
                    save_name num2str(id) '.mat'], ...
                    'run_info', 'save_info')
            catch
                mkdir(save_dir)
                save([save_dir, ...
                    save_name num2str(id) '.mat'], ...
                    'run_info', 'save_info')
            end

            disp([save_dir, save_name num2str(id) ' saved'])

        end % evi

    end % labi

end % roii


