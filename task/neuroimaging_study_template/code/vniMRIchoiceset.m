function [ choiceset ] = vniMRIchoiceset()
% value in context neuroimaging study
%  THIS SCRIPT CALLS THE FUNCTION: 'resample'
% updated 7/6/20

% create a distribution-like sampling around each level of mean EV
% where 2/3 of the trials within each level of mean EV will be sampled
% +/-$2 around indifference and 1/3 will be sampled  +/- >$2 around
% indifference

% task timing: each block is 12.025 minutes (36.08 minutes)

lnslope = -2; %the slope of the diagonals where we get safe and gain values
nT = 210; % there are 210 trials right now, 9 get added on at the end.
EVs = [5:5:25];% levels of EVs
evLev = size(EVs,2); % how many levels?

% what are the points on each diagonal of the EVs 5-25?
nPts = nT/evLev;

evPointsX = NaN(nPts,evLev);
evPointsY = NaN(nPts,evLev);

% setting up the "distribution":
evPointsX(:,1) = [linspace(0, 3,7) linspace(3.1, 6.9,28) linspace(7, 10,7)]; %ev 5
evPointsX(:,2) = [linspace(5, 8,7) linspace(8.1, 11.9,28) linspace(12, 15,7)]; %ev 10
evPointsX(:,3) = [linspace(10,13,7) linspace(13.1,16.9,28) linspace(17,20,7)]; %ev 15
evPointsX(:,4) = [linspace(15, 18,7) linspace(18.1, 21.9,28) linspace(22, 25,7)]; %ev 20
evPointsX(:,5) = [linspace(20, 23,7) linspace(23.1, 26.9,28) linspace(27, 30,7)]; %ev 25

evPointsY = flip(evPointsX)*-lnslope; % reverse X matrix * slope of the diagonals


% add noise
% create array of noise that will be applied to both X and Y coordinates

addNoise = -.5+1*rand(size(evPointsX)); % generate random numbers between -.5 and .5

evPointsXnoise = round(evPointsX + addNoise,2);
evPointsYnoise = round(evPointsY + 2*addNoise,2);

% put the two matrices into one
evXYcoords = [evPointsXnoise(:,1) evPointsYnoise(:,1) evPointsXnoise(:,2) evPointsYnoise(:,2) evPointsXnoise(:,3) evPointsYnoise(:,3) evPointsXnoise(:,4) evPointsYnoise(:,4) evPointsXnoise(:,5) evPointsYnoise(:,5)];
evXYcoords(evXYcoords < 0) = 0;  % any negatives, make them 0.
evPointsXnoise(evPointsXnoise < 0) = 0;
evPointsYnoise(evPointsYnoise < 0) = 0;


% Use a matrix design to ensure that we cover all possible shifts in a
% consistent and even way across all paticipants
n =0;

while 1
    n = n+1;
    mtx = ones(5,5);
    
    for r = 1:5
        mtx(r,:) = r;
        mtx(r,r) = NaN;
    end
    
    evs = NaN();
    for i = 1:20
        
        
        if i == 1
            evs(i) = randsample(5,1); %pick a column on first iteration
        else
            evs(i) = evnxt;  %otherwise ev becomes evnxt(previous row becomes column)
        end
        
        if i <20
            if ~any(find(~isnan(mtx(:,evs(i))))) % if there are no finite numbers in this column, start back at the top of while loop
                break 
            end
            
            evnxt = resample(find(~isnan(mtx(:,evs(i)))),1);
            mtx(evnxt, evs(i)) = NaN;
        end
        
        if all(isnan(mtx(:,evs(1))))
            mtx(evs(1),:) = NaN;
        end
    end
    
    if i == 20
        break % exit the while loop
    end


end % end while loop

evs = [evs, evs(1)];






rho = [2,1.5,.75,.5]; % a range of rho values
yvals = [0:5:55]; % gain values
xvalsSeekingH = (.5*(yvals.^rho(1))).^(1/rho(1)); % safe values raised to rho
xvalsSeekingL = (.5*(yvals.^rho(2))).^(1/rho(2));
xvalsAverseH = (.5*(yvals.^rho(3))).^(1/rho(3));
xvalsAverseL = (.5*(yvals.^rho(4))).^(1/rho(4));
xvalsNeutral = yvals*.5;




%{
 Each level of EV is made up of 4 runs (6,9,12,15 trials). 
 There are 42 possible trials in each run, taking care of 20 runs
 We need a 21st run because we will be going back up or down to the first run, 
 so we will add 9 more trials (nTri = 219) evenly spaced on the 
 final EV level [1:5:42) =  1  6 11 16 21 26 31 36 41.
 The only different across people for the LAST run is that they will be at different 
 EV levels but run length and index on the EV level will be the same.
 %} 
  
  %translate EVs (1-5) into the actual EVs (5-25)
  realEvs = evs;
  
  realEvs(realEvs==5) = EVs(5);
  realEvs(realEvs==4) = EVs(4);
  realEvs(realEvs==3) = EVs(3);
  realEvs(realEvs==2) = EVs(2);
  realEvs(realEvs==1) = EVs(1);

  %run lengths and blocks (for the scanner) each level of EV needs to have the same number of runs with the same amount of run lenghths (4 runs, lengths = 6,9,12,15), with the last run being 9 trials (we add that on at the end)
  %we also need to split the task into 3 blocks and ensure that a block does not end right before or after a shift.

  rL = [6,9,12,15]; % possible run lengths
  
  prestime = 2;
  choicetime = 2;
  feedbacktime=1;
  nBlock = 3;
  nTri = 219;
  nTb =nTri/nBlock; %73 trials per 3 blocks
  nShifts = 20;

  %isiTimes = 2:6
  isiTimes = 1.75:5.75;
  isiDistTimes = isiTimes.^-1.05;
  isiDistTimes = round(nTb*isiDistTimes/sum(isiDistTimes), 0); %number of times we see 2,3,4,5,6 (4 2 2 1 1)

  ISIs = NaN(nTb,1);

  ISIs = [repmat(isiTimes(1), isiDistTimes(1),1);
      repmat(isiTimes(2), isiDistTimes(2),1);
      repmat(isiTimes(3), isiDistTimes(3),1);
      repmat(isiTimes(4), isiDistTimes(4),1);      
      repmat(isiTimes(5), isiDistTimes(5),1)
      ];
  
  
  %itiTimes = 1:5;
  itiTimes = .75:4.75;
  itiDistTimes = itiTimes.^-1.05;
  itiDistTimes = round(nTb*itiDistTimes/sum(itiDistTimes), 0); 

  ITIs = NaN(nTb,1);

  ITIs = [repmat(itiTimes(1), itiDistTimes(1),1);
      repmat(itiTimes(2), itiDistTimes(2),1);
      repmat(itiTimes(3), itiDistTimes(3),1);
      repmat(itiTimes(4), itiDistTimes(4),1);      
      repmat(itiTimes(5), itiDistTimes(5),1)
      ];

  
n=0;
while 1
  n = n+1;
  runLength= NaN();

  runmtx = repmat(rL,evLev,1);
  
  for i = 1:20 % for the first 20 EVs (not 21, because we already know the last run is 9 trials)
    
    selectRow = evs(i); % row is the first EV
    selectCol = resample(find(~isnan(runmtx(selectRow,:)))',1);
    runLength(i)= runmtx(selectRow, selectCol);
    runmtx(selectRow,selectCol) = NaN;
    
  end 
   
  runLength = [runLength,9]; % add the last run 
  nRuns = size(runLength,2);
  
  % now distribute across the 3 blocks
  timeArray = NaN(nTb,4,nBlock); % column: isi, iti, run length, ev level
  timeArray(:,1,1) = ISIs(randperm(length(ISIs)));
  timeArray(:,2,1) = ITIs(randperm(length(ITIs)));
  timeArray(:,1,2) = ISIs(randperm(length(ISIs)));
  timeArray(:,2,2) = ITIs(randperm(length(ITIs)));
  timeArray(:,1,3) = ISIs(randperm(length(ISIs)));
  timeArray(:,2,3) = ITIs(randperm(length(ITIs)));
  
      
    % check that each block is same amt of time
   b1length = nTb*prestime + nTb*choicetime + nTb*feedbacktime + sum(timeArray(:,1,1) + sum(timeArray(:,2,1))); % block 1
   b2length = nTb*prestime + nTb*choicetime + nTb*feedbacktime + sum(timeArray(:,1,2) + sum(timeArray(:,2,2))); % block 2
   b3length = nTb*prestime + nTb*choicetime + nTb*feedbacktime + sum(timeArray(:,1,3) + sum(timeArray(:,2,3))); % block 3
  
   %length(unique([b1length,b2length,b3length]); % should all be the same

   
   addRuns = [];
   addShifts=[];
   
   for r=1:length(runLength)
       onerun = repmat(runLength(r),runLength(r),1);
       oneShift = repmat(realEvs(r), runLength(r),1);
       addRuns = [addRuns; onerun];
       addShifts = [addShifts;oneShift];
   end
  
    timeArray(:,3,1) = addRuns(1:length(timeArray));
    timeArray(:,3,2) = addRuns((length(timeArray)+1):(length(timeArray)*2));
    timeArray(:,3,3) = addRuns((length(timeArray)*2+1):length(addRuns));
    
    timeArray(:,4,1) = addShifts(1:length(timeArray));
    timeArray(:,4,2) = addShifts((length(timeArray)+1):(length(timeArray)*2));
    timeArray(:,4,3) = addShifts((length(timeArray)*2+1):length(addShifts));
    
   
   
    % we want the block to end somewhere in between a run 
    % ideally a block ends atleast 4 trials after a shift to see how
    % behavior might change after a shift and a block starts with atleast 4
    % trials left in the run (to allow people to adjust before a shift)
    
    if timeArray((nTb-4),4,1) == timeArray(4,4,2) && timeArray((nTb-4),4,2) == timeArray(4,4,3)
        break
    end
     
    
end %end while loop
   

  %disp(n); % how many times did it take?
  % check to see if our criteria are met
%   timeArray((nTb-4),4,1)
%   timeArray(4,4,2) % should be the same as above
%   timeArray((nTb-4),4,2)
%   timeArray(4,4,3) % should be the same


  % create a table of our choiceset
  choiceset = NaN(nTri,9);
  choiceset= array2table(choiceset, 'VariableNames',{'groundEV','evInd','runSize','block','isi','iti','riskyGain', 'riskyLoss', 'alternative'});

  count = 1; % start count at trial 1
  for i = 1:nRuns
    choiceset.groundEV(count:(count+runLength(i)-1)) = realEvs(i); % record the ground EV for every trial
    choiceset.evInd(count:(count+runLength(i)-1)) = evs(i);
    choiceset.runSize(count:(count+runLength(i)-1)) = runLength(i);
    count = count + runLength(i); % update the count
  end
  
  choiceset.isi = [timeArray(:,1,1); timeArray(:,1,2); timeArray(:,1,3)];
  choiceset.iti = [timeArray(:,2,1); timeArray(:,2,2); timeArray(:,2,3)];
  choiceset.block = [repmat(1,nTb,1); repmat(2,nTb,1); repmat(3,nTb,1)];
  choiceset.triBlock = [1:nTb 1:nTb 1:nTb]'; % trial within each block
  
  %Populate choice set with gains and losses such that we hit all the possible values for each level of EV given the number of trials in each run
 
  nTevLev = 1:42;
  valueInds = repmat(nTevLev',1,evLev);
  
  for n = 1:nT
      
      selectCol = choiceset.evInd(n); %column is level of EV for that trial, n
      
      selectRow = resample(find(~isnan(valueInds(:,selectCol))),1); % select a row
      
      choiceset.riskyGain(n) = evPointsYnoise(selectRow,selectCol);
      choiceset.alternative(n)= evPointsXnoise(selectRow,selectCol);
      choiceset.riskyLoss(n) = 0; %gain only task
      
      valueInds(selectRow,selectCol) = NaN ; % replace with NaN so we don't select this index again
      
  end
  
  %  sum(~isnan(valueInds)); % double check things went right
  
    
  % Add the last 9 trials: evenly spaced on the final EV level [1:5:42] =  1  6 11 16 21 26 31 36 41
  last9rowInd= 1:5:42;
  last9rowInd = last9rowInd(randperm(length(last9rowInd))); %shuffle the vector
  choiceset.riskyGain((nT+1):nTri) = evPointsYnoise(last9rowInd, evs(21));
  choiceset.alternative((nT+1):nTri) = evPointsXnoise(last9rowInd, evs(21));
  choiceset.riskyLoss((nT+1):nTri) = 0;
  
  % add some manipulation checks where safe is zero and the gamble is very
  % high
  ind = find(choiceset.riskyGain>50); % find trials where gain is > $50
  smpl = ind(randperm(size(ind,1),2)); % select two of those trials
  choiceset.alternative(smpl) = 0; % change the alternatives to 0
  % participants should 100% gamble here
  
  scaleBy = max(choiceset.riskyGain);
  choiceset.gainSC = round(choiceset.riskyGain/scaleBy,2);
  choiceset.lossSC = choiceset.riskyLoss/scaleBy;
  choiceset.altSC = round(choiceset.alternative/scaleBy,2);
  choiceset.EV = round((choiceset.riskyGain*.5) + (choiceset.riskyLoss*.5),2); % ev of the gamble
  choiceset.meanEV = round(((choiceset.riskyGain*.5) + choiceset.alternative)/2,2);
  

  
%   plot(choiceset.groundEV,'o')
%   hold on
%   plot(74,choiceset.groundEV(74),'r*')
%   plot(147,choiceset.groundEV(147),'r*')
  
  
end % end the function
