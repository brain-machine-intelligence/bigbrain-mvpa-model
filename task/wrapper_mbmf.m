function [] = wrapper_mbmf( sub_code, WAVE )
%%% FUNCTION TO CONTROLL THE PARADIGM

% test flag: if activated experiment will run with 4 blocks only
% test_switch = 'test';
test_switch = 'real';

%%% %%% %%%

% set session to save files
SESSION = ['wave_' num2str(WAVE) '_' num2str(sub_code)];

% set color randomisation
color_case = ceil(rand(1)*6);

% deactivate annoying PTB debugging
Screen('Preference', 'VisualDebugLevel', 3);
Screen('Preference', 'SuppressAllWarnings', 1);

% randomize stuff a little bit
srand = RandStream('mt19937ar', 'Seed', sub_code*WAVE); RandStream.setGlobalStream(srand);
% srand = RandStream('mt19937ar', 'Seed', sum(clock)); RandStream.setGlobalStream(srand);

% delete old logifles
save_folder = 'C:\ExpFiles\Sebastian\final_edit\result_save';
delete(fullfile(save_folder, '*'));

%%%%% START EXPERINMENT CONTROL %%%%%

disp(['this is participant number ' num2str(sub_code)]);
disp('press key to continue with training!');
pause;


%%% START TRAINING
SIMUL_arbitration_fmri3(sub_code, 1, 'pre', color_case, test_switch);

disp('training finished');
disp(' ');
disp('please wait for instructions!');

%%% WAIT TOGETHER TO START SESSION 1
% KbName('keynameswindows')
continue_key = 70; % this is key 'F'
press = 0;
while press == 0;
    [~, ~, kb_keycode] = KbCheck;
    if find(kb_keycode)==continue_key;
        press = 1;
    end
end

%%% START SESSION 1
SIMUL_arbitration_fmri3(sub_code, 1, 'fmri', color_case, test_switch);

disp('session 1 finished');
disp(' ');
disp('please wait for instructions!');
pause;


%%% WAIT TOGETHER TO START SESSION 2
continue_key = 71; % this is key 'G'
press = 0;
while press == 0;
    [~, ~, kb_keycode] = KbCheck;
    if find(kb_keycode)==continue_key;
        press = 1;
    end
end


%%% START SESSION 2
SIMUL_arbitration_fmri3(sub_code, 2, 'fmri', color_case, test_switch);

disp(' ');
disp('please wait...');

%%% COPY FILES TO THE MOTHER
source_path = save_folder;
target_path = 'N:\client_write\Sebastian\results';
mkdir(target_path, SESSION);
copyfile(fullfile(source_path, '*'), fullfile(target_path, SESSION));

%%% STUFF IS DONE
disp(' ');
disp('THANK YOU - EXPERIMENT FINISHED');

% end of function
end
