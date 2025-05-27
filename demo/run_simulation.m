
%% Add path

base_path = fullfile('..'); % If your current working directory is 'demo', use this. Otherwise, set your own path.
addpath(fullfile(base_path, 'data', 'human_behavior'))
addpath(fullfile(base_path, 'functions', 'analysis_simulation'))
addpath(fullfile(base_path, 'functions', 'analysis_simulation', 'behavior_record'))
addpath(fullfile(base_path, 'functions', 'utils'))

% addpath(genpath(fullfile(base_path, 'packages')))

% If your current working directory is 'demo', use this. Otherwise, set your own path that includes 'functions/analysis_simulation/behavior_record'.
cd('..'); 

%% Run MB agent simulation
% Run this section in a directory including 'behavior_record'.

disp('start')

saveDir = 'mb_virtual_episode';
nSimul = 100; % number of iteration
humanSimul = [0 0]; % fitting to human data (1 1) vs. virtual simulation (0 0)
mbmf = 1; % model-based learning agent

for id = 1:20
    
    try
        BATCH_ORI_FULL(id, nSimul, humanSimul, mbmf, saveDir);
    catch Err
        disp(Err.message)
    end
        
end

disp('end')


%% Run MF agent simulation
% Run this section in a directory including 'behavior_record'.

disp('start')

saveDir = 'mf_virtual_episode';
nSimul = 100; % number of iteration
humanSimul = [0 0]; % fitting to human data (1 1) vs. virtual simulation (0 0)
mbmf = 2; % model-free learning agent

for id = 1:20
        
    try
        BATCH_ORI_FULL(id, nSimul, humanSimul, mbmf, saveDir);
    catch Err
        disp(Err.message)
    end
        
end

disp('end')




