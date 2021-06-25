function DOSE = runDose(DataFile, utilityFunction, prior_typeRho, prior_typeLambda,prior_typeMu, Choice, outputFile, numQus, softmaxGroup, saveMatlab, header)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Prior Convergence  Experiment:
%%% Running DOSE procedure on Econographics Survey Respondents to test for
%%% convergence aggregating the posterior distribution of the parameters over
%%% all subjects and using this as the prior for next run of the experiment
%%%
% DataFile  = string of **csv** name in input_files which contains subjects data
%             or dataset/matrix already in memory
%             Column 1 must contain integer subject ID
%             Colums 2-5 must contain question list and user choice [z x y a]
% utilityFunction = specifies name of softmax file to be used
% Priors    = Int denoting the number of elements in the support of the priors
%             or a string 'custom' denoting custom priors are to be used. Files
%             for priors must be in the input_files folder and named
%             {rho, lam, mu}_parameters.csv with the first row containing
%             the parameter values and the second row containg the probabilities
% Choice    = String one of {'optimal', 'original'}. Optimal produce optimal
%             question ordering using the Kullback-Leibler divergence,
%             otherwise uses the original order and just calculates bayes stats
%
% outputFile = string of **csv** name to be placed in output_files with results
% numQus = number of questions to be run, set as `maximum; if want all questions
% softmaxGroup = whether to specify softmax for all subjects at once ('allSubjects', or by individual ('bySubject')
% 				Former is faster, but will not work on v.large datasets due to matrix size
% saveMatlab = 'yes' to save matlab file as well (this can be dataintensive)
% header = write 'noheader' if .csv input has no header row  or 'header'

%%%      By: Judah Okwuobi, Jonathan Chapman
%%% Updated: 5/17/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Subject data loaded as a table containing SubID, question set, and subj answer
try
    if ischar(DataFile)
		disp('reading csv file')
		if strcmp(header, 'noheader')
			disp('reading file with no header')
			SubjData = csvread(['input_files/' DataFile '.csv'],0,0);
			%SubjData = dataset('FILE', ['input_files/' DataFile '.csv'], 'Delimiter', ',');
		else
			disp('first row of .csv a header row')
			SubjData = csvread(['input_files/' DataFile '.csv'],1,0);
        end
    else
        SubjData = double(DataFile);
        SubjData = SubjData(:,[1:6 15 7:14]);
    end
catch
    warning(['The data set should come in as a csv or data object that is' ...
             ' convertible to double in Matlab'])
end

% Run Experiment for all Subjects
numSubs = numel(unique(SubjData(:,1)));
sub_ids = unique(SubjData(:,1));

% % Priors should come in as a matrix/double col1 = values, col2 = probabilites
% % These are uniform priors over a standard range for each parameter
if isnumeric(prior_typeRho)
	if strcmp(utilityFunction,'CRRA')
	% Priors are constructed by inputing prior ranges
		priors.rho = makePriors(0.2, 1.661589941, prior_typeRho);
		priors.lam = makePriors(0, 4.561383585, prior_typeLambda);
		priors.mu  = makePriors(0, 8, prior_typeMu);
	elseif strcmp(utilityFunction,'sShapedCARA') || strcmp(utilityFunction,'concaveCARA')
	% Priors are constructed by inputing prior ranges
		priors.rho = makePriors(-0.661589941,0.8, prior_typeRho);
		priors.lam = makePriors(0, 4.561383585, prior_typeLambda);
		priors.mu  = makePriors(0, 8, prior_typeMu);
	end	
elseif strcmp(prior_typeRho,'custom')
    priors.rho = csvread('input_files/priors/rho_parameters.csv');
    priors.lam = csvread('input_files/priors/lambda_parameters.csv');
    priors.mu  = csvread('input_files/priors/mu_parameters.csv');

else
    warning('Priors  must be one of {Int, "custom", "wep", "shfPosterior"}')
    return
end

numRhos = length(priors.rho);
numLambdas = length(priors.lam);
numMus = length(priors.mu);

% Cutt off conditon: In current paridigm we prespecify number of runs,
% however in future we may use a cuttoff conditions for convergence
maxruns = 1; % number of runs before procedure should halt
run  = 1;

while run <= maxruns
    run    %display run number
    if run == 2
        disp('Beginning Feedback Sequence')
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% SECTION 1: Pre-calculate the probability space for all combinations
    %%% of questions and parameter values with associated probabilities .
    %%% Parameter vals are given in first run, but are calculated from
    %%% aggregate posterior distribution in subsequent runs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    disp('Initial Run on Uniform Priors')
    % Run Dose for all Subjects in Parallel
    posteriors = cell(1,numSubs); % Storing posteriors for passing to TimePref
    estimates = cell(1,numSubs); %Storing aggregate  estimates
    params  = cell(1,numSubs); % Storing Parameter priors and posteriors for all subjects

    %probabilities for the prior distribution
    probTensor = zeros(numMus, numLambdas, numRhos);
    checkExist=exist('priorTensor');
    if checkExist == 1
        priorTensor = priorTensor;
    else
        disp 'priorTensor is not given. Object is being created...'
        priorTensor = kron(priors.rho(2,:), kron(priors.lam(2,:), priors.mu(2,:)));
        priorTensor = reshape(priorTensor(:), [numMus numLambdas numRhos]);
    end

    % If selected above, softmax Array is calculated 
	% once prior to each run using all the
    % unique questions across all subjects to obtain a master q.list
	% Faster than calculating by subject but may 
	% not work for v.large datasets with many questions
	% Since matrix will be too big
    % Function call returns struct with two fields containing the
    % softmaxArray, and the matrix of all unique questions from across
    % all subjects
    % Runs>1 use the aggregated posteriors from previously adjacent run
	if strcmp(softmaxGroup, 'allSubjects')
		if strcmp(utilityFunction,'CRRA')
			softmaxArray = softmaxCRRA(SubjData(:,2:4), priors, 0.5);
		elseif strcmp(utilityFunction,'concaveCARA')
			softmaxArray = softmaxConcaveCARA(SubjData(:,2:4), priors, 0.5);
		elseif	strcmp(utilityFunction,'sShapedCARA')
			softmaxArray = softmaxsShapedCARA(SubjData(:,2:4), priors, 0.5);
		else
		warning('Need to specify "CRRA", "concaveCARA" or "sShapedCARA" utility')
		return
		end
		disp('End of Softmax');
	end

    % parfor (i=1:numSubs, workers) %change upper limit to numSubs
    for i = 1:numSubs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% SECTION 2: Get the parameter estimates for each subject given their
    %%% data of answers to a question set. In the first run the priors are
    %%% given, however in subsequent runs they are constructed from
    %%% aggregated posteriors.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sub_id = sub_ids(i);
        disp(['Subject  ' num2str(sub_id)]);

        % Nx4 Matrix of subject data including question list and user
        % choice over N rounds of question from Survey
        subjdata = SubjData(ismember(SubjData(:,1), sub_id), [2:5 7]);
		
		% Calculate softmax for the individual if selected in the options
		% slower than calculating all together, but may be needed for 
		% very large datasets when re-estimating
		if strcmp(softmaxGroup, 'bySubject')
			if strcmp(utilityFunction,'CRRA')
				softmaxArray = softmaxCRRA(subjdata(:,1:3), priors, 0.5);
			elseif strcmp(utilityFunction,'concaveCARA')
				softmaxArray = softmaxConcaveCARA(subjdata(:,1:3), priors, 0.5);
			elseif	strcmp(utilityFunction,'sShapedCARA')
				softmaxArray = softmaxsShapedCARA(subjdata(:,1:3), priors, 0.5);
			else
			warning('Need to specify "CRRA", "concaveCARA" or "sShapedCARA" utility')
			return
			end
			disp('End of Softmax');
		end
		
		% Number of questions
		if isnumeric(numQus)
			numQu = numQus;
		elseif strcmp(numQus, 'maximum')	
			numQu = size(subjdata,1);
		else
            warning(['for numQu enter a number of questions or "maximum"'])
		end
        % Get estimates given question ordering
        if strcmp(Choice, 'optimal')
            disp('optimal ordering')
            results = shChoiceOptimal(subjdata, priors, priorTensor, sub_id, softmaxArray,numQu);
        elseif strcmp(Choice,'original')
            disp('original ordering')
            results = shChoiceNoKL(subjdata, priors, priorTensor, sub_id, softmaxArray,numQu);
        else
            warning(['question order type must be one of "optimal" or "original"'...
                     ' otherwise the procedure might produce nonsensical results'])
        end
        estimatesTime(i) = toc(t);
        disp(['Subject  ' num2str(sub_id)  '  Elapsed Time:  ' num2str(estimatesTime(i))])

        %%% Probability tensor is aggregated across subjects
        probTensor = probTensor + results.probTensor;

        %%% Each cell contains an individuals result array
        estimates{i} = results.subjestimates;

        %%% Store each subjects user ID and posterior and prior estimates
        %posteriors{i,1} = deal(sub_id, results.posteriors);
        tempPost = cell(1,2);
        tempPost{1} = sub_id;
        tempPost{2}= results.posteriors;
        posteriors{i} = tempPost;

    end


    % Converting all cell estimates to matrix for saving
    estimates = cell2mat(transpose(estimates));

    run = run + 1;
    %clearvars -except priors run maxruns SubjData numRhos prior_type t % results
end

probTensor = probTensor./numSubs;
 DOSE.estimates  = estimates;
 DOSE.probTensor = probTensor;
 DOSE.posteriors = reshape([posteriors{:}], [2,numSubs]);
 DOSE.priors     = priors;

% save filename with dates
% Order of columns is: 
% K-L number=1 if no K-L number chosen
% Answer of 999=missing
%header1={'WinAmount' 'LossAmount' 'SureAmount' 'KL_Number'  'Avg_Rho' 'SE_Rho' 'Avg_Lam' ...
           % 'SE_Lam' 'Avg_Mu' 'SE_Mu' 'Answer' 'Round' 'SubID' 'qID'};
FileName=['output_files\' outputFile ,datestr(now, 'dd-mmm-yyyy'),'.csv']
dlmwrite(FileName,estimates,'precision',11)

if strcmp(saveMatlab, 'yes')
	disp('saving matlab file')
	FileName=['output_files\' outputFile ,datestr(now, 'dd-mmm-yyyy'),'.mat']
	save(FileName)
else
	disp('not saving matlab file, set option to yes at start of file to save.')
end
 f= DOSE;