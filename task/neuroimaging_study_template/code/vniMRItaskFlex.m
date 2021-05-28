function [subjdata] = vniMRItaskFlex(subjID, reStart, isreal)
% VNI neuroimaging script updated on 6/8/20 to be more flexible in the
% event the study needs to be stopped part way through.
% Hayley Brooks, University of Denver

% subjID: '###'
% reStart: 0 = starting new study; 1 = reStarting with block 1;
% 2=reStarting with block 2; 3 = reStarting with block 3
% isreal: 0 = test; 1 = real;

% each block is 12.025 minutes ~ 12 minutes and 1.5 seconds (not including
% the 15 seconds of fixation before the block starts)

try
    homepath = 'C:\Users\inc\Documents\Sokol-Hessner Lab\';
    cd(homepath);
catch
    homepath = pwd;
    cd(homepath);
end

clc % clear command window

if nargin < 3
    isreal = 1; %assume real (if its a test, the trials will be shorter)
    if nargin <2
        reStart = 0; % assume we are not restarting a run
        if nargin < 1 % if number of input arguments for this function is less than one,
            subjID = '000'; % set the default subject ID in case we don't include one.
        end
    end
end

if length(subjID) >3
    error('Subject ID too long')
elseif length(subjID) < 3
    error('Subject ID too short')
end

fname = sprintf('VNI_behavior_VNI%s_%.4f.mat',subjID,now);
% this will create a unique data file name so that data will not be
% overwrote each iteration

%----Prep----%
%Initialize rand to a different state each time:
if verLessThan('matlab','7.12')
    rand('twister',sum(100*clock))
    fprintf('Set rng using oldest method\n');
elseif verLessThan('matlab','7.14')
    RandStream.setDefaultStream ...
        (RandStream('mt19937ar','seed',sum(100*clock))); % 'shuffle' makes the seed based on the current time
    fprintf('Set rng using intermediate method\n')
else
    rng('shuffle'); % 'shuffle' makes the seed based on the current time
    fprintf('Set rng using newest method\n')
end

if isreal
    HideCursor; % remove the cursor so it doesn't just show up on our screen
end

Screen('Preference', 'SkipSyncTests', 1);
% problems with display refresh rate synchronization tests, so set to 1 to
% skip tests and avoid errors, these tests are more important for tasks
% where timing needs to be very precise

[wind,rect] = Screen('OpenWindow', max(Screen('Screens'))); % open window on the max number of screens available (if two screens,max = 1)
Priority(MaxPriority(wind)); % make sure the computer will default to updating the window you care about
h = rect(4); w = rect(3); % define screens height and width


Screen('TextFont',wind,'Courier'); % set the font to something non-offensive & portable
Screen('TextSize',wind,40); % set the font size to something reasonable

blk = BlackIndex(wind); % automatically pull out the black, white, and gray CLUTs
wht = WhiteIndex(wind);
gry = GrayIndex(wind);

Screen('FillRect', wind, blk); % set background
Screen('Flip',wind); %THIS YOU DON'T HAVE TO RE-SET EACH FLIP


DrawFormattedText(wind,'Setting up...','center','center',wht,40); % put up something on the screen while we finish prep
Screen('Flip',wind);
WaitSecs(2); % display content for two seconds

%keys
KbName('UnifyKeyNames'); % help PTB deal with different keyboards.
esc_key_code = KbName('ESCAPE'); % this is the key we'll use to exit the experiment if we want to
continue_key_code = KbName('space'); % this is key to continue
scanner_trigger = KbName('5%');
resp_keys = {'1!','2@'}; % 1 for selecting left,2 for selecting right
resp_key_codes = KbName(resp_keys); % get the key codes for the response keys

if isreal
    nT = 219; % number of trials for real task
    nB = 3;
    triBlockStart = [1,74,147];
else
    nT = 12;% number of trials for testing
    nB = 3;
    triBlockStart = [1,5,9];
    testTriBlock = [1:4,1:4,1:4];
end

% set the timing for each trial event (% ISI and ITI will be defined
% differently depending on reStart input)
prestime = 2; % number of seconds the image is presented
choicetime = 2; % number of seconds I have to respond
feedbacktime = 1; % number of seconds for feedback

% Big diff is going to be if we have to reStart blocks 1, 2, and/or 3
if reStart == 0 % we are starting fresh!
    cs = vniMRIchoiceset; %generate new choiceset
    
    feedbackdelay = cs.isi; % get isis from new choice set.
    ititime= cs.iti; % get itis from new choice set
    
    % Prep the struct which we'll eventually save out with all the study data
    subjdata.subjID = subjID; % subject ID
    subjdata.nT = nT; % number of trials
    subjdata.cs = cs; % store the choiceset
    subjdata.cs.choice = nan(size(cs,1),1); %used size of cs instead of nT for testing purposes when nT = 6
    subjdata.cs.outcome = nan(size(cs,1),1);
    subjdata.cs.loc = nan(size(cs,1),1);% what side of the screen were gamble and alt presented;loc = 1: gamble on left, alternative or right, loc = 2 alt on left, gamble on the right
    subjdata.cs.response = cell(size(cs,1),1);
    subjdata.cs.subjectIndex = repmat(subjID,size(cs,1),1); % put subject index in every row - will want this for analysis
    subjdata.cs.RTs = nan(size(cs,1),1);
    subjdata.cs.trial = (1:size(cs,1))'; % add trial number
    
    subjdata.ts = struct(); %create a struct to store variable timing events
    subjdata.ts.blockStart = nan(3,1);
    
    subjdata.params = struct(); % create a struct to store the preset timing events
    subjdata.params.prestime = prestime;
    subjdata.params.choicetime = choicetime;
    subjdata.params.feedbacktime = feedbacktime;
    subjdata.params.feedbackdelay = feedbackdelay;
    subjdata.params.ititime = ititime;
    
    nTp = 3; %3 practice trials only happening for reStarting
    subjdata.practice = struct();
    subjdata.practice.riskyGain = [27.63,33.15,29.87];
    subjdata.practice.riskyLoss = [0,0,0];
    subjdata.practice.alternative = [10.15,18.35,15.23];
    subjdata.practice.choice = nan(nTp,1); %used size of cs instead of nT for testing purposes when nT = 6
    subjdata.practice.outcome = nan(nTp,1);
    subjdata.practice.loc = nan(nTp,1);% what side of the screen were gamble and alt presented;loc = 1: gamble on left, alternative or right, loc = 2 alt on left, gamble on the right
    subjdata.practice.response = cell(nTp,1);
    subjdata.practice.RTs = nan(nTp,1);
    
    
else  %if we are reStarting a block
    
    subjdata = pullData(subjID);
    subjdata = subjdata.subjdata;
    
    feedbackdelay = nan(nT,1);
    feedbackdelay = subjdata.cs.isi;
    
    ititime = nan(nT,1);
    ititime = subjdata.cs.iti;
    
    
    %testText = class(subjdata); % this works!
    %testText = subjdata.subjID; % this also works!
    %testText = num2str(size(subjdata.cs)); % this works: '219 22'
    %testText = num2str(size(subjdata.cs.iti));% this also works: '219 1'
    %testText = num2str(feedbackdelay(1)); %this works when feedbackdelay
    %is a vector of NaNs and when assigning subjdata.cs.isi to
    %feedbackdelay
    %testText = num2str(ititime(1)); % this works when ititime is a vector
    % of NaNs and when assigning subjdata.cs.iti to ititime
    %DrawFormattedText(wind,testText,'center','center',wht,40);
    %Screen('Flip',wind);
    %WaitSecs(2);
    
    
end % if reStart == 0;

fid = fopen(sprintf('VNI_behavior_VNI%s_%.4f.txt',subjID,now),'w'); % open a new text file, set it to "write" status



% setting up some stimulus stuff %
% rect dimensions for placing shapes and text
% rect dims: rect(1)=left, rect(2)=top, rect(3)=right, rect(4)=bottom
xCenter = rect(3)/2; % center of x axis
yCenter = rect(4)/2; % center of y axis
radius = round(w*.15); % define radius

lCent = w/4; % center of the left side
rCent = 3*w/4; % center of the right side

%dimensions for left rectangle and right rectangle
lRect = [lCent-radius h/2-radius lCent+radius h/2+radius];
rRect = [rCent-radius h/2-radius rCent+radius h/2+radius];


% The start of the script will be different depending on reStart input

try
    if reStart == 0 % if we are starting a new session, do intro and practice
        %-------------Brief INTRO-----------%
        DrawFormattedText(wind,'As discussed in the instructions, this task consists of 219 trials. You will choose between a gamble and a guaranteed alternative.\n\n\n Press the left button to choose the option on the left and the right button to choose the option on the right.\n\n Press either button to continue ','center','center',wht,40);
        Screen('Flip',wind);
        WaitSecs(2);
        while 1
            [keyIsDown,~,keyCode] = KbCheck;
            if (keyIsDown && size(find(keyCode),2)==1)
                if keyCode(esc_key_code)
                    error('Experiment aborted by user'); % allow aborting the study here
                elseif any(keyCode(resp_key_codes)) % participant can start study when ready
                    break % change screen as soon as they respond
                end
            end
        end
        
        %---------------------Practice-----------------------%
        DrawFormattedText(wind,'We will now begin the practice task. \n\n\n When you''re ready to begin, press the left or right button. ','center','center',wht,40);
        Screen('Flip',wind);
        WaitSecs(2);
        while 1
            [keyIsDown,~,keyCode] = KbCheck;
            if (keyIsDown && size(find(keyCode),2)==1)
                if keyCode(esc_key_code)
                    error('Experiment aborted by user'); % allow aborting the study here
                elseif any(keyCode(resp_key_codes)) % participant can start study when ready
                    break % change screen as soon as they respond
                end
            end
        end
        
        %--Start practice--%
        subjdata.practice.studystart = GetSecs; % log the study start time if it's the first trial
        for t = 1:nTp
            subjdata.practice.trialStart(t) = GetSecs; %when does the trial start? there will be a slight delay between this and stimStart(t)
            loc = randperm(2,1); %randomly select location of options on the screen. loc = 1: gamble on left, alternative or right, loc = 2 alt on left, gamble on the right
            subjdata.practice.loc(t) = loc; %store location
            %------stimulus presentation---------%
            %prepare stimulus%
            Screen('FillOval',wind,wht,[lRect' rRect']); %draw two ovals on the left and right side of the screen
            if loc == 1
                Screen('DrawLine', wind, blk, lRect(1),yCenter,lRect(3), yCenter, 5); %draw line in the middle of gamble on left side of screen
                %place text
                DrawFormattedText(wind,sprintf('$%.2f',subjdata.practice.riskyGain(t)),'center','center',blk,[],[],[],[],[],[lRect(1) lRect(2) lRect(3) lRect(4)- radius]); %gain text
                DrawFormattedText(wind,sprintf('$%.2f',subjdata.practice.riskyLoss(t)),'center','center',blk,[],[],[],[],[],[lRect(1) lRect(2)+radius lRect(3) lRect(4)]); %loss text
                DrawFormattedText(wind,sprintf('$%.2f',subjdata.practice.alternative(t)),'center','center',blk,[],[],[],[],[],rRect); % alt text
            elseif loc == 2
                Screen('DrawLine', wind, blk, rRect(1),yCenter, rRect(3), yCenter, 5); % draw line in middle of gamble on right side of screen
                %place text
                DrawFormattedText(wind,sprintf('$%.2f',subjdata.practice.alternative(t)),'center','center',blk,[],[],[],[],[],lRect); %alt text
                DrawFormattedText(wind,sprintf('$%.2f',subjdata.practice.riskyGain(t)),'center','center',blk,[],[],[],[],[],[rRect(1) rRect(2) rRect(3) rRect(4)-radius]); %gain text
                DrawFormattedText(wind,sprintf('$%.2f',subjdata.practice.riskyLoss(t)),'center','center',blk,[],[],[],[],[],[rRect(1) rRect(2)+radius rRect(3) rRect(4)]); % loss text
            end
            DrawFormattedText(wind,'OR','center','center',wht); % place OR between both ovals
            
            
            subjdata.practice.stimStart(t) = Screen('Flip',wind,[],1); %show the stimuli
            
            % ready the response period stuff that will show up in the response
            % collection period
            DrawFormattedText(wind, 'Left','center', 'center',wht,[],[],[],[],[],[lRect(1),lRect(4) lRect(3) lRect(4)+lRect(2)]);
            DrawFormattedText(wind, 'Right','center', 'center',wht,[],[],[],[],[],[rRect(1),rRect(4) rRect(3) rRect(4)+rRect(2)]);
            
            while GetSecs - subjdata.practice.studystart < t*prestime + (t-1)*(choicetime+feedbacktime) + sum(feedbackdelay(1:(t-1)))+ sum(ititime(1:(t-1)))
                [keyIsDown,~,keyCode] = KbCheck;
                %if keyIsDown
                if (keyIsDown && size(find(keyCode),2) ==1)
                    if keyCode(esc_key_code)
                        error('Experiment aborted by user'); % allow aborting during stimulus presentation
                    end
                end
            end
            
            %----------response collection-----------%
            %display keys for choices on screen
            subjdata.practice.choiceStart(t) = Screen('Flip',wind); % show the v and n on the screen, nows ps can respond
            
            while GetSecs - subjdata.practice.studystart < t*(prestime + choicetime) + (t-1)*(feedbacktime) + sum(feedbackdelay(1:(t-1))) + sum(ititime(1:(t-1)))
                [keyIsDown,~,keyCode] = KbCheck;
                %outp(adOut,1); %send marking pulse
                %if keyIsDown
                if (keyIsDown && size(find(keyCode),2) ==1)
                    if keyCode(esc_key_code)
                        error('Experiment aborted by user'); % allow aborting the study here
                    elseif any(keyCode(resp_key_codes)) % IF the pressed key matches a response key...
                        subjdata.practice.RTs(t) = GetSecs - subjdata.practice.choiceStart(t); % record RT
                        subjdata.practice.extraTime(t) = choicetime - subjdata.practice.RTs(t); % store the extra time left over once p makes a choice
                        %subjdata.ts.feedbackdelayStart(t) = GetSecs; %now the isi starts
                        %outp(adOut,0); % turn off marker - tell biopac that event is over
                        break % change screen as soon as they respond
                    end
                end
            end
            
            
            %------ISI------%
            Screen('gluDisk',wind,wht,xCenter,yCenter,3); % make a white dot that will be displayed during ISI
            subjdata.practice.isiStart(t) = Screen('Flip',wind); %now isi starts
            
            if any(keyCode(resp_key_codes))
                subjdata.practice.response{t} = KbName(keyCode); % record response
                if subjdata.practice.loc(t) == 1 && strcat(subjdata.practice.response{t}(1)) == '1' %if gamble is on the left and p selected v
                    subjdata.practice.choice(t) = 1; % then they gambled
                elseif subjdata.practice.loc(t)==2 && strcat(subjdata.practice.response{t}(1)) == '2' % if gamble is on the right and p selected n
                    subjdata.practice.choice(t) = 1; %then they gambled
                else
                    subjdata.practice.choice(t) = 0; % otherwise they did not gamble
                end
            else
                %subjdata.cs.response{t} = NaN; %did not respond in time
                subjdata.practice.choice(t) = NaN; %did not respond in time
            end
            
            while GetSecs - subjdata.practice.studystart < t*prestime + min([subjdata.practice.RTs(t) choicetime]) + (t-1)*(choicetime+feedbacktime) + sum(feedbackdelay(1:t))+ sum(ititime(1:(t-1)));
                [keyIsDown,~,keyCode] = KbCheck;
                if (keyIsDown && size(find(keyCode),2) ==1)
                    %if keyIsDown
                    if keyCode(esc_key_code)
                        error('Experiment aborted by user'); % allow aborting during stimulus presentation
                    end
                end
            end
            
            %-----FEEDBACK-----%
            
            % prepare feedback stimulus during ISI
            if ~isnan(subjdata.practice.choice(t))
                if subjdata.practice.choice(t) == 1 %if gamble taken
                    OC = randi(2,1); % generate a random number, 1 or 2; if 1, outcome is gain. if 2, outcome is loss
                    if OC == 1 %if OC is 1, gain is outcome
                        subjdata.practice.outcome(t) = subjdata.practice.riskyGain(t); %store gain amount
                        if loc == 1 %if gamble is on left
                            Screen('FillArc',wind,wht,lRect,270,180); %display gain amount on left
                            DrawFormattedText(wind,sprintf('$%.2f',subjdata.practice.riskyGain(t)),'center','center',blk,[],[],[],[],[],[lRect(1) lRect(2) lRect(3) lRect(4)-radius]);
                        elseif loc ==2 % if gamble is on the right
                            Screen('FillArc',wind,wht,rRect,270,180); %display gain amount on right
                            DrawFormattedText(wind,sprintf('$%.2f',subjdata.practice.riskyGain(t)),'center','center',blk,[],[],[],[],[],[rRect(1) rRect(2) rRect(3) rRect(4)-radius]);
                        end
                    elseif OC ==2  %if OC is 2, loss is the outcome
                        subjdata.practice.outcome(t) = subjdata.practice.riskyLoss(t); % store loss amount
                        if loc ==1 %if gamble is on the right
                            subjdata.practice.outcome(t) = subjdata.practice.riskyLoss(t); % store loss amount
                            Screen('FillArc',wind,wht,lRect,90,180); %display outcome
                            DrawFormattedText(wind,sprintf('$%.2f',subjdata.practice.riskyLoss(t)),'center','center',blk,[],[],[],[],[],[lRect(1) lRect(2)+radius lRect(3) lRect(4)]);
                        elseif loc == 2 % if
                            Screen('FillArc',wind,wht,rRect,90,180);
                            DrawFormattedText(wind,sprintf('$%.2f',subjdata.practice.riskyLoss(t)),'center','center',blk,[],[],[],[],[],[rRect(1) rRect(2)+radius rRect(3) rRect(4)]);
                        end
                    end
                else
                    subjdata.practice.outcome(t) = subjdata.practice.alternative(t); %store alt amount
                    if loc == 1 %if alt is on the right side of screen
                        Screen('FillOval', wind, wht,rRect);
                        DrawFormattedText(wind,sprintf('$%.2f',subjdata.practice.alternative(t)),'center','center',blk,[],[],[],[],[],rRect);
                    else
                        
                        Screen('FillOval', wind, wht,lRect);
                        DrawFormattedText(wind,sprintf('$%.2f',subjdata.practice.alternative(t)),'center','center',blk,[],[],[],[],[],lRect);
                    end
                end
            else
                DrawFormattedText(wind,'You did not respond in time','center','center', wht, 40);
                subjdata.practice.outcome(t) = NaN;
                
            end % end this feedback for loop; now feedback has been prepped and is ready to be displayed when isi is over
            
            
            subjdata.practice.otcStart(t) =  Screen('Flip',wind); % show feedback
            
            
            Screen('gluDisk',wind,wht,xCenter,yCenter,3); % make a white dot that will be displayed during iti
            
            
            while GetSecs - subjdata.practice.studystart < t*(prestime+feedbacktime) + (t-1)*(choicetime) + sum(feedbackdelay(1:t)) +sum(ititime(1:(t-1))) + min([subjdata.practice.RTs(t) choicetime])
                [keyIsDown,~,keyCode] = KbCheck;
                if (keyIsDown && size(find(keyCode),2) ==1)
                    %if keyIsDown
                    if keyCode(esc_key_code)
                        error('Experiment aborted by user'); % allow experiment aborting here
                    end
                end
            end
            
            
            %subjdata.ts.otcEnd(t) = GetSecs; % end of outcome display
            
            % ------ ITI ------ %
            
            subjdata.practice.itiStart(t) = Screen('Flip',wind); %show white dot on the screen
            
            %fprintf(fid,'riskyGain %0.2f, safe, %0.2f, loc, %g, choice, %g, outcome, %0.2f, RT, %0.2f, ISI, %g, ITI, %g\n', subjdata.practice.riskyGain(t),subjdata.practice.alternative(t),subjdata.practice.loc(t), subjdata.practice.choice(t), subjdata.practice.outcome(t), subjdata.practice.RTs(t)); %save txt file of what we have
            save(fname,'subjdata'); % save what data we have
            
            %wait a certain amount of time until iti ends
            while GetSecs - subjdata.practice.studystart < t*(prestime+choicetime+feedbacktime) + sum(feedbackdelay(1:t)) + sum(ititime(1:t))  % control the entire experiment length on a per-trial basis
                [keyIsDown,~,keyCode] = KbCheck;
                if (keyIsDown && size(find(keyCode),2) ==1)
                    %if keyIsDown
                    if keyCode(esc_key_code)
                        error('Experiment aborted by user'); % allow experiment aborting here
                    end
                end
            end
            
        end %ends the loop for the practice trials.
        
        subjdata.practice.studystop = GetSecs; % mark the end-time
        
        
        DrawFormattedText(wind,'Practice Complete!\n\n\n If you have any questions, please hold them and the experimenter will be with you in a few minutes. \n\n\n Press either button to continue.','center','center', wht, 40);
        Screen('Flip',wind);
        WaitSecs(2);
        while 1
            [keyIsDown,~,keyCode] = KbCheck;
            if (keyIsDown && size(find(keyCode),2) ==1)
                if keyCode(esc_key_code)
                    error('Experiment aborted by user'); % allow aborting the study here
                elseif any(keyCode(resp_key_codes)) % participant can start study when ready
                    break % change screen as soon as they respond
                end
            end
        end
    end %reStart ==0 // practice and intro complete
    

    % Once we are ready to start a block, show waiting for experimenter screen:
    DrawFormattedText(wind, 'Waiting for experimenter.', 'center', 'center', wht,40);
    Screen('Flip', wind); % update display to show break_message
    WaitSecs(1);
    while 1
        [keyIsDown,~,keyCode] = KbCheck;
        if (keyIsDown && size(find(keyCode),2) ==1)
            if keyCode(esc_key_code)
                error('Experiment aborted by user'); % allow aborting the study here
            elseif keyCode(continue_key_code) %
                break % change screen as soon as they respond
            end
        end
    end
    DrawFormattedText(wind, 'Waiting for scanner.', 'center', 'center', wht,40);
    Screen('Flip', wind); % update display to show break_message
    WaitSecs(1);
    % wait for scanner pulse trigger to start
    while keyCode(scanner_trigger) == 0 % loop continues until scanner trigger
        
        [~, ~, keyCode] = KbCheck; % check keyboard for a response, ~'s used as dummy placeholders
        
        % wait 1 ms before checking the keyboard again to prevent
        % overload of the machine at elevated Priority()
        WaitSecs(0.001);
    end
    DrawFormattedText(wind, '+', 'center', 'center', wht,40);
    Screen('Flip', wind); % update display to show break_message
    WaitSecs(15.0);
    
    %Now the real task can begin!
    
    if reStart == 0
        subjdata.ts.studystart = GetSecs; % log the study start time
        b = 1;
        triStart = 1;
    else
        b = reStart; % block is indicated by the reStart value
        triStart = triBlockStart(b); % start on the first trial of current block
    end
    
    
    for t = triStart:nT
        
        if isreal
            triBlock = subjdata.cs.triBlock(t);
        else
            triBlock = testTriBlock(t);
        end
        
        %----------interblock break----------%
        if t > ((nT/nB)*b)
            b=b+1;
            DrawFormattedText(wind, 'There will now be a brief break before continuing the rest of the task. \n\n Waiting for experimenter.', 'center', 'center', wht,40);
            Screen('Flip', wind); % update display to show break_message
            subjdata.ts.blockBreakStart(b) = GetSecs;
            WaitSecs(1);
            while keyCode(continue_key_code) == 0
                [keyIsDown,~,keyCode] = KbCheck;
                if (keyIsDown && size(find(keyCode),2) ==1)
                    if keyCode(esc_key_code)
                        error('Experiment aborted by user'); % allow aborting the study here
                    elseif keyCode(continue_key_code) %
                        break % change screen as soon as they respond
                    end
                end
                % wait 1 ms before checking the keyboard again to prevent
                % overload of the machine at elevated Priority()
            end
            
            DrawFormattedText(wind, 'Waiting for scanner.', 'center', 'center', wht,40);
            Screen('Flip', wind); % update display to show break_message
            WaitSecs(1);
            % wait for scanner pulse trigger to start
            while keyCode(scanner_trigger) == 0 % loop continues until scanner trigger
                
                [~, ~, keyCode] = KbCheck; % check keyboard for a response, ~'s used as dummy placeholders
                
                % wait 1 ms before checking the keyboard again to prevent
                % overload of the machine at elevated Priority()
                WaitSecs(0.001);
            end
            
            DrawFormattedText(wind, '+', 'center', 'center', wht,40);
            Screen('Flip', wind); % update display to show break_message
            WaitSecs(15.0);
        end %if t > ((nT/nB)*b)
        
        subjdata.ts.trialStart(t) = GetSecs;
        
        %store the start of each block as when the first trial starts
        if t == triBlockStart(1)
            subjdata.ts.blockStart(1) = subjdata.ts.trialStart(t);
        elseif t == triBlockStart(2)
            subjdata.ts.blockStart(2) = subjdata.ts.trialStart(t);
        elseif t == triBlockStart(3)
            subjdata.ts.blockStart(3) = subjdata.ts.trialStart(t);
        end
        
        loc = randperm(2,1); %randomly select location of options on the screen. loc = 1: gamble on left, alternative or right, loc = 2 alt on left, gamble on the right
        subjdata.cs.loc(t) = loc; %store location
        %------stimulus presentation---------%
        %prepare stimulus%
        Screen('FillOval',wind,wht,[lRect' rRect']); %draw two ovals on the left and right side of the screen
        if loc == 1
            Screen('DrawLine', wind, blk, lRect(1),yCenter,lRect(3), yCenter, 5); %draw line in the middle of gamble on left side of screen
            %place text
            DrawFormattedText(wind,sprintf('$%.2f',subjdata.cs.riskyGain(t)),'center','center',blk,[],[],[],[],[],[lRect(1) lRect(2) lRect(3) lRect(4)- radius]); %gain text
            DrawFormattedText(wind,sprintf('$%.2f',subjdata.cs.riskyLoss(t)),'center','center',blk,[],[],[],[],[],[lRect(1) lRect(2)+radius lRect(3) lRect(4)]); %loss text
            DrawFormattedText(wind,sprintf('$%.2f',subjdata.cs.alternative(t)),'center','center',blk,[],[],[],[],[],rRect); % alt text
        elseif loc == 2
            Screen('DrawLine', wind, blk, rRect(1),yCenter, rRect(3), yCenter, 5); % draw line in middle of gamble on right side of screen
            %place text
            DrawFormattedText(wind,sprintf('$%.2f',subjdata.cs.alternative(t)),'center','center',blk,[],[],[],[],[],lRect); %alt text
            DrawFormattedText(wind,sprintf('$%.2f',subjdata.cs.riskyGain(t)),'center','center',blk,[],[],[],[],[],[rRect(1) rRect(2) rRect(3) rRect(4)-radius]); %gain text
            DrawFormattedText(wind,sprintf('$%.2f',subjdata.cs.riskyLoss(t)),'center','center',blk,[],[],[],[],[],[rRect(1) rRect(2)+radius rRect(3) rRect(4)]); % loss text
        end
        DrawFormattedText(wind,'OR','center','center',wht); % place OR between both ovals
        
        subjdata.ts.stimStart(t) = Screen('Flip',wind,[],1); %show the stimuli
        
        % ready the response period stuff that will show up in the response
        % collection period
        DrawFormattedText(wind, 'Left','center', 'center',wht,[],[],[],[],[],[lRect(1),lRect(4) lRect(3) lRect(4)+lRect(2)]);
        DrawFormattedText(wind, 'Right','center', 'center',wht,[],[],[],[],[],[rRect(1),rRect(4) rRect(3) rRect(4)+rRect(2)]);
        
        while GetSecs - subjdata.ts.blockStart(b) < triBlock*prestime + (triBlock-1)*(choicetime+feedbacktime) + sum(feedbackdelay(triBlockStart(b):(t-1))) + sum(ititime(triBlockStart(b):(t-1)))
            [keyIsDown,~,keyCode] = KbCheck;
            if (keyIsDown && size(find(keyCode),2) ==1)
                if keyCode(esc_key_code)
                    error('Experiment aborted by user'); % allow aborting during stimulus presentation
                end
            end
        end
        
        %----------response collection-----------%
        %display keys for choices on screen
        subjdata.ts.choiceStart(t) = Screen('Flip',wind); % show the v and n on the screen, nows ps can respond
        
        while GetSecs - subjdata.ts.blockStart(b) < triBlock*(prestime + choicetime) + (triBlock-1)*(feedbacktime) + sum(feedbackdelay(triBlockStart(b):(t-1))) + sum(ititime(triBlockStart(b):(t-1)))
            [keyIsDown,~,keyCode] = KbCheck;
            %outp(adOut,1); %send marking pulse
            if (keyIsDown && size(find(keyCode),2) ==1)
                if keyCode(esc_key_code)
                    error('Experiment aborted by user'); % allow aborting the study here
                elseif any(keyCode(resp_key_codes)) % IF the pressed key matches a response key...
                    subjdata.cs.RTs(t) = GetSecs - subjdata.ts.choiceStart(t); % record RT
                    subjdata.ts.extraTime(t) = choicetime - subjdata.cs.RTs(t); % store the extra time left over once p makes a choice
                    %outp(adOut,0); % turn off marker - tell biopac that event is over
                    break % change screen as soon as they respond
                end
            end
        end
        
        
        %------ISI------%
        Screen('gluDisk',wind,wht,xCenter,yCenter,3); % make a white dot that will be displayed during ISI
        subjdata.ts.isiStart(t) = Screen('Flip',wind); %now isi starts
        
        if any(keyCode(resp_key_codes))
            subjdata.cs.response{t} = KbName(keyCode); % record response
            if subjdata.cs.loc(t) == 1 && strcat(subjdata.cs.response{t}(1)) == '1' %if gamble is on the left and p selected v
                subjdata.cs.choice(t) = 1; % then they gambled
            elseif subjdata.cs.loc(t)==2 && strcat(subjdata.cs.response{t}(1)) == '2' % if gamble is on the right and p selected n
                subjdata.cs.choice(t) = 1; %then they gambled
            else
                subjdata.cs.choice(t) = 0; % otherwise they did not gamble
            end
        else
            %subjdata.cs.response{t} = NaN; %did not respond in time
            subjdata.cs.choice(t) = NaN; %did not respond in time
        end
        
        
        while GetSecs - subjdata.ts.blockStart(b) < triBlock*(prestime) + min([subjdata.cs.RTs(t) choicetime]) + (triBlock-1)*(choicetime+feedbacktime) + sum(feedbackdelay(triBlockStart(b):t)) + sum(ititime(triBlockStart(b):(t-1)))
            [keyIsDown,~,keyCode] = KbCheck;
            if (keyIsDown && size(find(keyCode),2) ==1)
                %if keyIsDown
                if keyCode(esc_key_code)
                    error('Experiment aborted by user'); % allow aborting during stimulus presentation
                end
            end
        end
        
        
        %-----FEEDBACK-----%
        
        % prepare feedback stimulus during ISI
        if ~isnan(subjdata.cs.choice(t))
            if subjdata.cs.choice(t) == 1 %if gamble taken
                OC = randi(2,1); % generate a random number, 1 or 2; if 1, outcome is gain. if 2, outcome is loss
                if OC == 1 %if OC is 1, gain is outcome
                    subjdata.cs.outcome(t) = subjdata.cs.riskyGain(t); %store gain amount
                    if loc == 1 %if gamble is on left
                        Screen('FillArc',wind,wht,lRect,270,180); %display gain amount on left
                        DrawFormattedText(wind,sprintf('$%.2f',subjdata.cs.riskyGain(t)),'center','center',blk,[],[],[],[],[],[lRect(1) lRect(2) lRect(3) lRect(4)-radius]);
                    elseif loc ==2 % if gamble is on the right
                        Screen('FillArc',wind,wht,rRect,270,180); %display gain amount on right
                        DrawFormattedText(wind,sprintf('$%.2f',subjdata.cs.riskyGain(t)),'center','center',blk,[],[],[],[],[],[rRect(1) rRect(2) rRect(3) rRect(4)-radius]);
                    end
                elseif OC ==2  %if OC is 2, loss is the outcome
                    subjdata.cs.outcome(t) = subjdata.cs.riskyLoss(t); % store loss amount
                    if loc ==1 %if gamble is on the right
                        subjdata.cs.outcome(t) = subjdata.cs.riskyLoss(t); % store loss amount
                        Screen('FillArc',wind,wht,lRect,90,180); %display outcome
                        DrawFormattedText(wind,sprintf('$%.2f',subjdata.cs.riskyLoss(t)),'center','center',blk,[],[],[],[],[],[lRect(1) lRect(2)+radius lRect(3) lRect(4)]);
                    elseif loc == 2 % if
                        Screen('FillArc',wind,wht,rRect,90,180);
                        DrawFormattedText(wind,sprintf('$%.2f',subjdata.cs.riskyLoss(t)),'center','center',blk,[],[],[],[],[],[rRect(1) rRect(2)+radius rRect(3) rRect(4)]);
                    end
                end
            else
                subjdata.cs.outcome(t) = subjdata.cs.alternative(t); %store alt amount
                if loc == 1 %if alt is on the right side of screen
                    Screen('FillOval', wind, wht,rRect);
                    DrawFormattedText(wind,sprintf('$%.2f',subjdata.cs.alternative(t)),'center','center',blk,[],[],[],[],[],rRect);
                else
                    
                    Screen('FillOval', wind, wht,lRect);
                    DrawFormattedText(wind,sprintf('$%.2f',subjdata.cs.alternative(t)),'center','center',blk,[],[],[],[],[],lRect);
                end
            end
        else
            DrawFormattedText(wind,'You did not respond in time','center','center', wht, 40);
            subjdata.cs.outcome(t) = NaN;
            
        end % end this feedback for loop; now feedback has been prepped and is ready to be displayed when isi is over
        
        
        subjdata.ts.otcStart(t) =  Screen('Flip',wind); % show feedback
        
        Screen('gluDisk',wind,wht,xCenter,yCenter,3); % make a white dot that will be displayed during iti
        
        
        while GetSecs - subjdata.ts.blockStart(b) < triBlock*(prestime+feedbacktime) + (triBlock-1)*(choicetime) + sum(feedbackdelay(triBlockStart(b):t)) + sum(ititime(triBlockStart(b):(t-1))) + min([subjdata.cs.RTs(t) choicetime])
            [keyIsDown,~,keyCode] = KbCheck;
            if (keyIsDown && size(find(keyCode),2) ==1)
                %if keyIsDown
                if keyCode(esc_key_code)
                    error('Experiment aborted by user'); % allow experiment aborting here
                end
            end
        end
        
        
        %subjdata.ts.otcEnd(t) = GetSecs; % end of outcome display
        
        % ------ ITI ------ %
        
        subjdata.ts.itiStart(t) = Screen('Flip',wind); %show white dot on the screen
        
        
        fprintf(fid,'riskyGain %0.2f, safe, %0.2f, loc, %g, choice, %g, outcome, %0.2f, RT, %0.2f, ISI, %g, ITI, %g, trial, %g, block, %g, triBlock, %g\n', subjdata.cs.riskyGain(t),subjdata.cs.alternative(t),subjdata.cs.loc(t), subjdata.cs.choice(t), subjdata.cs.outcome(t), subjdata.cs.RTs(t), subjdata.params.feedbackdelay(t), subjdata.params.ititime(t),t,b,triBlock); %save txt file of what we have
        save(fname,'subjdata'); % save what data we have
        
        %wait a certain amount of time until iti ends
        while GetSecs - subjdata.ts.blockStart(b) < triBlock*(prestime+choicetime+feedbacktime) + sum(feedbackdelay(triBlockStart(b):t)) + sum(ititime(triBlockStart(b):t))
            [keyIsDown,~,keyCode] = KbCheck;
            if (keyIsDown && size(find(keyCode),2) ==1)
                if keyCode(esc_key_code)
                    error('Experiment aborted by user'); % allow experiment aborting here
                end
            end
        end
        
    end %for t = triStart:nT
    
    subjdata.ts.studystop = GetSecs; % mark the end-time
    
    %---Outcome selection and wrap up---%
    DrawFormattedText(wind,'Task Complete!','center','center',wht,40); % put up a screen placeholder
    Screen('Flip',wind);
    WaitSecs(2);
    
    DrawFormattedText(wind, 'Randomly selecting 1 outcome...', 'center', 'center', wht, 40); %display text while we randomly select choice
    Screen('Flip', wind);
    WaitSecs(2);
    
    
    subjdata.realOutcomesInds = find(isfinite(subjdata.cs.outcome)); %find in outcomes the rows with finite numbers, ignore NaNs
    subjdata.payoutTrial = subjdata.realOutcomesInds(randi(length(subjdata.realOutcomesInds),1)); % randomly select one of the trials to pay out
    subjdata.fullOutcome = subjdata.cs.outcome(subjdata.payoutTrial); % pull out the amount of the outcome
    subjdata.payout = subjdata.fullOutcome * .25; % for VIC sliding scale
    
    
    %display results of random generator
    DrawFormattedText(wind,sprintf('Randomly selected trial: #%d.\n\nThe outcome of this trial is $%.2f.\n\n As a result, you will receive $%.2f.\n\n Please wait for further instructions.', subjdata.payoutTrial, subjdata.fullOutcome, subjdata.payout), 'center', 'center', wht,40);
    Screen('Flip', wind);
    WaitSecs(1);
    while keyCode(continue_key_code) == 0 % wait for experimenter to press space bar
        [keyIsDown,~,keyCode] = KbCheck;
        if (keyIsDown && size(find(keyCode),2) ==1)
            if keyCode(esc_key_code)
                error('Experiment aborted by user'); % allow aborting the study here
            elseif keyCode(continue_key_code) %
                break % change screen as soon as they respond
            end
        end
        % wait 1 ms before checking the keyboard again to prevent
        % overload of the machine at elevated Priority()
    end
    
    fclose(fid); % close the text file
    save(fname,'subjdata'); % save what we have
    
    
    sca % clear the screen (note: this also clears things like "HideCursor"...)
    Priority(0); % reset screen priority
    
    %fprintf('\nExperiment took %0.5f seconds when it should have taken %0.1f seconds.\n',subjdata.ts.studystop-subjdata.ts.studystart, nT*(prestime+choicetime+feedbacktime) + sum(feedbackdelay) + sum(ititime))
    
catch ME % what to do if something breaks in the "try" statment, ME will show you the error
    
    fclose(fid); % close the text file
    save(fname,'subjdata'); % save what data we have
    sca % clear the screen (note: this also clears things like "HideCursor"...)
    Priority(0); % return screen priority to default
    rethrow(ME); %show the error
end
end % end this function