function [subjdata] = vniPREmriTask(subjID, isreal)
%Pre-scan task (brief overview of instructions and practice trials)


%cd S:\Projects\VIC\task\

%cd /Volumes/shlab/Projects/VNI/task/ % set working directory

try
    homepath = 'C:\Users\sokolhessnerlab\Desktop';
    cd(homepath);
catch
    homepath = pwd;
    cd(homepath);
end

clc % clear command window

if nargin < 2
    isreal = 1; % assume real (if its a test, the trials will be shorter)
    if nargin < 1 % if number of input arguments for this function is less than one,
        subjID = '000'; % set the default subject ID in case we don't include one.
    end
end


if length(subjID) >3
    error('Subject ID too long')
elseif length(subjID) < 3
    error('Subject ID too short')
end

fname = sprintf('VNI_preMriTask_VNI%s_%.4f.mat',subjID,now);


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
Screen('Preference', 'SkipSyncTests', 1); % for testing purposes when sync is out of wack

[wind,rect] = Screen('OpenWindow', max(Screen('Screens'))); % open window on the max number of screens available (if two screens,max = 1)
Priority(MaxPriority(wind)); % make sure the computer will default to updating the window you care about
h = rect(4); w = rect(3); % define screens height and width

Screen('TextFont',wind,'Courier'); % set the font to something non-offensive & portable
Screen('TextSize',wind,40); % set the font size to something reasonable

blk = BlackIndex(wind); % automatically pull out the black, white, and gray CLUTs
wht = WhiteIndex(wind);
gry = GrayIndex(wind);

Screen('FillRect', wind, blk); % set background <--
%THIS YOU DON'T HAVE TO RE-SET EACH FLIP
Screen('Flip',wind);

DrawFormattedText(wind,'Setting up...','center','center',wht,40); % put up something on the screen while we finish prep
Screen('Flip',wind);
WaitSecs(2); % display content for two seconds

%keys
KbName('UnifyKeyNames'); % help PTB deal with different keyboards.
esc_key_code = KbName('ESCAPE'); % this is the key we'll use to exit the experiment if we want to
%trig_key_code = KbName('SPACE'); % this is trigger key
%exp_key_code = KbName('Return'); % this is experimenter key
exp_key_code = KbName('Space'); % this is experimenter key
resp_keys = {'1!' '2@'}; % x for selecting left,n for selecting right
resp_key_codes = KbName(resp_keys); % get the key codes for the response keys

%-------- study set up ------------%
%addpath('C:\Users\sokolhessnerlab\Desktop\VNI task\');
%THIS WILL BE WHEREVER THE TASK IS LOCATED DURING THE ACTUAL EXPERIMENT


nT=5;

% set the timing for each trial event
prestime = 2; % number of seconds the image is presented
choicetime = 2; % number of seconds I have to respond
feedbacktime = 1; % number of seconds for feedback
feedbackdelay = [1,2,1,1,3]'; % get isis from choice set.
ititime= [4,1,2,2,1]'; % get itis from choice set

fid = fopen(sprintf('VNI_preMriTask_VNI%s_%.4f.txt',subjID,now),'w'); % open a new text file, set it to "write" status

% Prep the struct which we'll eventually save out with all the study data
subjdata.subjID = subjID; % subject ID
subjdata.nT = nT; % number of trials
subjdata.cs = struct; % store the choiceset
subjdata.cs.riskyGain = [8.64, .64, 48.60, 21.17, 39.01];
subjdata.cs.riskyLoss = [0 0 0 0 0];
subjdata.cs.alternative = [5.53, 8.80, 14.96, 18.76, 29.10];
subjdata.cs.choice = nan(nT,1); %used size of cs instead of nT for testing purposes when nT = 6
subjdata.cs.outcome = nan(nT,1);
subjdata.cs.loc = nan(nT,1);% what side of the screen were gamble and alt presented;loc = 1: gamble on left, alternative or right, loc = 2 alt on left, gamble on the right
subjdata.cs.response = cell(nT,1);
subjdata.cs.subjectIndex = repmat(subjID,nT,1); % put subject index in every row - will want this for analysis
subjdata.cs.RTs = nan(nT,1);
subjdata.cs.trial = (1:nT)'; % add trial number


subjdata.ts = struct(); %create a struct to store variable timing events

subjdata.params = struct(); % create a struct to store the preset timing events
subjdata.params.prestime = prestime;
subjdata.params.choicetime = choicetime;
subjdata.params.feedbacktime = feedbacktime;
subjdata.params.feedbackdelay = feedbackdelay;
subjdata.params.ititime = ititime;

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


% OK, we're ready - give control to the experimenter
DrawFormattedText(wind,'Waiting for the experimenter','center','center',wht,40); %draw text
Screen('Flip',wind);
WaitSecs(2);
while 1
    [keyIsDown,~,keyCode] = KbCheck;
    if (keyIsDown && size(find(keyCode),2) ==1)
        if keyCode(esc_key_code)
            error('Experiment aborted by user'); % allow aborting the study here
        elseif any(keyCode(exp_key_code)) % experimenter key required
            break % change screen as soon as they respond
        end
    end
end

try % try running the study
    %-------------Brief INTRO-----------%
    DrawFormattedText(wind,'As discussed in the instructions, this task consists of 219 trials. You will choose between a gamble and a guaranteed alternative.\n\n\n Press the blue key to choose the left and yellow key to choose the right.\n\n\n If you have any questions, please ask the experimenter now.','center','center',wht,40);
    Screen('Flip',wind);
    WaitSecs(2);
    while 1
        [keyIsDown,~,keyCode] = KbCheck;
        if (keyIsDown && size(find(keyCode),2) ==1) % making sure that when participants press more than 1 button, the screen does not change
            if keyCode(esc_key_code)
                error('Experiment aborted by user'); % allow aborting the study here
            elseif any(keyCode(exp_key_code)) % experimenter key required
                break % change screen as soon as they respond
            end
        end
    end
    
    %---------------------Practice-----------------------%
    DrawFormattedText(wind,'We will now begin the practice task. \n\n\n When you''re ready to begin, press the blue or yellow key. ','center','center',wht,40);
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
    
    
    
    %--Start practice--%
    
    
    subjdata.ts.studystart = GetSecs; % log the study start time if it's the first trial
    for t = 1:nT
        subjdata.ts.trialStart(t) = GetSecs; %when does the trial start? there will be a slight delay between this and stimStart(t)
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
        
        while GetSecs - subjdata.ts.studystart < t*prestime + (t-1)*(choicetime+feedbacktime) + sum(feedbackdelay(1:(t-1)))+ sum(ititime(1:(t-1)))
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
        subjdata.ts.choiceStart(t) = Screen('Flip',wind); % show the v and n on the screen, nows ps can respond
        
        while GetSecs - subjdata.ts.studystart < t*(prestime + choicetime) + (t-1)*(feedbacktime) + sum(feedbackdelay(1:(t-1))) + sum(ititime(1:(t-1)))
            [keyIsDown,~,keyCode] = KbCheck;
            %outp(adOut,1); %send marking pulse
            %if keyIsDown
            if (keyIsDown && size(find(keyCode),2) ==1)
                if keyCode(esc_key_code)
                    error('Experiment aborted by user'); % allow aborting the study here
                elseif any(keyCode(resp_key_codes)) % IF the pressed key matches a response key...
                    subjdata.cs.RTs(t) = GetSecs - subjdata.ts.choiceStart(t); % record RT
                    subjdata.ts.extraTime(t) = choicetime - subjdata.cs.RTs(t); % store the extra time left over once p makes a choice
                    %subjdata.ts.feedbackdelayStart(t) = GetSecs; %now the isi starts
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
        
        while GetSecs - subjdata.ts.studystart < t*prestime + min([subjdata.cs.RTs(t) choicetime]) + (t-1)*(choicetime+feedbacktime) + sum(feedbackdelay(1:t))+ sum(ititime(1:(t-1)));
            [keyIsDown,~,keyCode] = KbCheck;
            if (keyIsDown && size(find(keyCode),2) ==1)
                %if keyIsDown
                if keyCode(esc_key_code)
                    error('Experiment aborted by user'); % allow aborting during stimulus presentation
                end
            end
        end
        
        %
        %         if doSCR ==1
        %             outp(adOut,0); %tell biopac to turn marker off
        %         end
        %
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
        
        %         if doSCR == 1
        %             outp(adOut,1); % turn on marker - tell biopac that event is starting
        %         end
        
        Screen('gluDisk',wind,wht,xCenter,yCenter,3); % make a white dot that will be displayed during iti
        
        
        while GetSecs - subjdata.ts.studystart < t*(prestime+feedbacktime) + (t-1)*(choicetime) + sum(feedbackdelay(1:t)) +sum(ititime(1:(t-1))) + min([subjdata.cs.RTs(t) choicetime])
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
        
        %         if doSCR == 1
        %             outp(adOut,0); % turn off marker - tell biopac that event is ending
        %         end
        
        fprintf(fid,'riskyGain %0.2f, safe, %0.2f, loc, %g, choice, %g, outcome, %0.2f, RT, %0.2f, ISI, %g, ITI, %g\n', subjdata.cs.riskyGain(t),subjdata.cs.alternative(t),subjdata.cs.loc(t), subjdata.cs.choice(t), subjdata.cs.outcome(t), subjdata.cs.RTs(t), subjdata.params.feedbackdelay(t), subjdata.params.ititime(t)); %save txt file of what we have
        save(fname,'subjdata'); % save what data we have
        
        %wait a certain amount of time until iti ends
        while GetSecs - subjdata.ts.studystart < t*(prestime+choicetime+feedbacktime) + sum(feedbackdelay(1:t)) + sum(ititime(1:t))  % control the entire experiment length on a per-trial basis
            %while GetSecs - subjdata.ts.stimStart(1) < (t*(prestime+choicetime+feedbacktime+feedbackdelay))+sum(ititime(1:t)) % control the entire experiment length on a per-trial basis
            [keyIsDown,~,keyCode] = KbCheck;
            if (keyIsDown && size(find(keyCode),2) ==1)
                %if keyIsDown
                if keyCode(esc_key_code)
                    error('Experiment aborted by user'); % allow experiment aborting here
                end
            end
        end
        
    end %ends the loop for the task.
    
    subjdata.ts.studystop = GetSecs; % mark the end-time
    
    
    %experimenter take back control of task
    DrawFormattedText(wind,'Practice Complete!\n\n\n If you have any questions, ask the experimenter now. \n\n\n\n\n Waiting for experimenter.','center','center', wht, 40);
    Screen('Flip',wind);
    WaitSecs(2);
    while 1
        [keyIsDown,~,keyCode] = KbCheck;
        if (keyIsDown && size(find(keyCode),2) ==1)
            if keyCode(esc_key_code)
                error('Experiment aborted by user'); % allow aborting the study here
            elseif any(keyCode(exp_key_code)) % participant can start study when ready
                break % change screen as soon as they respond
            end
        end
    end
    
    fclose(fid); % close the text file
    
    %save(sprintf('vic_behavior_VIC%s_%.4f.mat',subjID,now),'subjdata'); % save out that subject data struct only with a non-overwriting name
    save(fname,'subjdata'); % save what we have
    
    while 1
        [keyIsDown,~,keyCode] = KbCheck;
        if (keyIsDown && size(find(keyCode),2) == 1)
            %if keyIsDown
            if keyCode(esc_key_code)
                error('Experiment aborted by user'); % if they pressed the escape key, end the study
            elseif any(keyCode(exp_key_code))
                break % otherwise leave the while loop and finish!
            end
        end
    end
    
    sca % clear the screen (note: this also clears things like "HideCursor"...)
    Priority(0); % reset screen priority
    
    fprintf('\nExperiment took %0.5f seconds when it should have taken %0.5f seconds.\n',subjdata.ts.studystop-subjdata.ts.studystart, nT*(prestime+choicetime+feedbacktime) + sum(feedbackdelay) + sum(ititime))
    
catch ME % what to do if something breaks in the "try" statment, ME will show you the error
    
    fclose(fid); % close the text file
    save(fname,'subjdata'); % save what data we have
    sca % clear the screen (note: this also clears things like "HideCursor"...)
    Priority(0); % return screen priority to default
    rethrow(ME); %show the error
    
end
end

