function [D] = n_back(t,n,init,num_targets)
% one block of n-back task
% t: number of trials for which a correct/incorrect action is defined
% n: n-back


stim_root =  pwd ;
if init
    config_display(0, 5, [0 0 0], [1 1 1], 'Arial', 60, 4);
    config_keyboard(100,5,'exclusive');
    start_cogent;
    cgloadbmp(3, fullfile(stim_root,'instructions_n_back.bmp') )
    cgdrawsprite(3,0,0)
    cgflip(3)
    cgpencol(1,1,1)
    waitkeydown(inf)
end

for ii = fliplr(1:5)
    cgtext(['Report ' int2str(n) '-back targets: start in ' int2str(ii) ' secs'],0,0); cgflip(0,0,0); pause(1)
end
pt  = 0.5;      % from  jaeggi, studer-luethi et al (2010)
iti = 2.5;        % was 2.5 in the above;

r   = get_random_sequence(t,n,100000000,num_targets);
for trl = 1:t
    t0  = time;
    clearkeys;
    if n == 1;
        cgrect(0,512,1280,30);
    end
    cgtext(r(trl),0,0); cgflip(0,0,0); pause(pt);
    if n == 1;
        cgrect(0,512,1280,30);
    end
    cgflip(0,0,0);
    pause(iti);
    
    [r_ rt_]    = response(t0); % dependent measures
    D.r(trl)    = r_;
    D.rt(trl)   = rt_;
    if trl > n
        D.target(trl) = strcmp(r(trl),r(trl-n)); % indicates a correct n-back target (against which to assess subjects' performance).
    end
end
D.TP = D.r( D.target); % total true positives
D.FP = D.r(~D.target); % total false positive
D.TN = ~D.r(~D.target); % total true negatives
D.FN = ~D.r(D.target); % total false negatives

D.RT = D.rt(D.target); % RT for true positives

stop_cogent


function [r] = get_random_sequence(t,n,II,num_targets)
% t: pick numbers for the total number of trials t
% n as in n-back
% II = 1 returns a random string sequence
% II = big uses brute force to return a random string subject to constraint that there are exactly 4 targets

for ii = 1:II;
    set     = char(['A':'H']);          % jaeggi, studer-luethi et al used 8 images not letters
    nset    = length(set);
    i       = ceil(nset*rand(1,t)) ;    % with repeat
    r       = set(i);
    
    for trl = 1:t %
        if trl > n
            target(trl) = strcmp(r(trl),r(trl-n)); % indicates a correct n-back target (against which to assess subjects' performance).
        end
    end
    if sum(target) == num_targets; % want e.g. 6 targets in the block?
        break
    end
end

function [r rt] = response(t0)
readkeys
[ key, t, n ] = getkeydown;
if ~n  % no keypress
    r  = 0; rt = 0;
else
    r  = 1; rt = t(1) - t0;
end
