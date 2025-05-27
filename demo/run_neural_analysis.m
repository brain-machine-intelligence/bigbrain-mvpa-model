

%% Add path

base_path = fullfile('..'); % If your current working directory is 'demo', use this. Otherwise, set your own path.
addpath(fullfile(base_path, 'data', 'subj_masked_EPI'))
addpath(fullfile(base_path, 'data', 'human_behavior'))
addpath(fullfile(base_path, 'data', 'fmri_supp'))
addpath(fullfile(base_path, 'functions', 'analysis_neural_data'))
addpath(fullfile(base_path, 'functions', 'utils'))

addpath(genpath(fullfile(base_path, 'packages')))

% If your current working directory is 'demo', use this. Otherwise, set your own path that includes 'data'.
cd('..'); 


%% Run decoding analysis
% Run this section in a directory including 'data'.

disp('start')

parfor id = 1:20
    
    try
        ana_decoding(id)
    catch Err
        disp(Err.message)
    end
        
end

disp('end');

%% Run shattering analysis
% Run this section in a directory including 'data'.

disp('start')

parfor id = 1:20
    
    try
        ana_shattering(id)
    catch Err
        disp(Err.message)
    end
        
end

disp('end');

%% Run CCGP analysis
% Run this section in a directory including 'data'.

disp('start')

parfor id = 1:20
    
    try
        ana_ccgp(id)
    catch Err
        disp(Err.message)
    end
        
end

disp('end');


%%






