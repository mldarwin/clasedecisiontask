function f = shChoiceNoKL(subjdata, priors, priorTensor, sub_id, softmaxArray,numQ)
%%% Function for running the dose experiment design on subjects given a
%%% question list and subject''s answers to the questions. The optimal question
%%% order is obtained using the KL divergence metric, then posteriors are
%%% calculated, and estimates over all rounds of questions are obtained.
% subjdata       = 3 column  matrix  of questions from a subject [z x y]
% priors         = 1x1 struct with 3 fields containing 2xN matrices of parameter
%                  distributions priors.rho, priors.lam, priors.mu
% sub_id         = Subject ID number
% softmaxArray   = N-D array of probabilities used in calculating bayesian stats
%%%      By: Judah Okwoubi
%%% Updated: 7/24/2107

%%% The questions and the softmax for that question are not removed from
%%% the respective objects, as is done in the regular SH_Choice
t = tic;
% Hardcoded values
p = 0.5;

% Remove non-answers
%subjdata = subjdata(subjdata(:,4) ~=999,:);

% Num Questions
%numQ = size(subjdata,1);

% Load in prior distributions from "priors" argument
% Each array has two rows. The first holds the value of the parameter
% The second holds the probabilities associated with each value.

rhoArray = priors.rho ; lamArray = priors.lam ; muArray = priors.mu ;
numRhos = length(rhoArray); numLambdas = length(lamArray); numMus = length(muArray);

% Parameter probability tensor from Kronecker prod of intial prior probs
% generated and reshaped  with Mu, Lam, Rho, along dimensions 1,2,3 respectively
priorTensor = kron(priors.rho(2,:), kron(priors.lam(2,:), priors.mu(2,:)));
priorTensor = reshape(priorTensor(:), [numMus numLambdas numRhos]);

% Create an array that will hold our questions that we would ask, the KL number
% for the question, the average and standard error of lam, rho, and mu, as well
% as the question id
resultArray = ones(numQ, 14);

% Colums are as follows :
% |1. GainQ   |2. LossQ   |3. SureQ |4. KL# (always 1 here) |5-10. Estimates/SD |11. Choice
% |12. Round# |13. SubjID |14. qID |15. Identifier for Frydman if SHF dataset

% Question ordering is already determinined in BROAD procedure in Survey
% Question set, round#, and subject ID are inputed into the result array
subjdata = subjdata(1:numQ,:);
resultArray(1:numQ, [1:3 11:13 14]) = ...
    [subjdata(1:numQ,1:4) transpose(1:numQ) repmat(sub_id,numQ,1) subjdata(:,end)];
subjdata = subjdata(:,1:4);

% Loop over the number of questions and input results into array
for round = 1:numQ

    %Question Ordering with KL divergence is no longer needed

    user_answer = subjdata(round,4);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Softmx array is calculated from all questions used on all subjects so
    % the question index number must be obtained from a matching on the
    % master question list.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    question_num = find(ismember(softmaxArray.master_qlist, ...
                        subjdata(round, 1:3), 'rows'));
    resultArray(round, 14) = question_num;

    % Store the softmax for the question asked.
    % If the users answer is yes, then we use the softmax. If the user
    % answers no then we set softmax = 1 - sofmax
    if user_answer ==1
        tempsoftmaxArray = softmaxArray.array(:,:,:, question_num);
    elseif user_answer ==0
        tempsoftmaxArray = 1 - softmaxArray.array(:,:,:, question_num);
    else
        if round == 1 % Case where subject doesn't answer on the first of Q's
            tempEstimates = [];
            tempParams = {'rho', 'lam', 'mu'};
            for i = 1:3
                tempEstimates(2*i-1) = mean(priors.(tempParams{i})(1,:));
                tempEstimates(2*i) = std(priors.(tempParams{i})(1,:));
            end
            resultArray(round,5:10) = tempEstimates;

        else
            resultArray(round,5:10) = resultArray(round-1, 5:10);
        end
        continue
    end

    % Calculate the new array of probabilities for every combination of
    % lambda, rho and mu
    % P(rho,lam,mu | ans )  for all param values are now in tempProbMatrix
    tempProbMatrix = probUpdater(tempsoftmaxArray, priorTensor);
    tempProbMatrix = reshape(tempProbMatrix, [numMus numLambdas numRhos]);

    % Set the posterior dist ribution as the new prior for the next round
    priorTensor = tempProbMatrix;

    % Marginal distributions are calculated for obtaining point estimates
    tempRho = sum(sum(tempProbMatrix,1),2); rhoArray(2,:) = tempRho(:);
    tempLam = sum(sum(tempProbMatrix,3),1); lamArray(2,:) = tempLam(:);
    tempMu  = sum(sum(tempProbMatrix,3),2); muArray(2,:)  = tempMu(:);


    % calculate and store the average values and SDs for rho, lambda, and mu

    % Rho: Mean and Variance
    avgRho = sum(times(rhoArray(1,:), rhoArray(2,:)));
    sdRho  = sum(rhoArray(2,:).*(rhoArray(1,:)-avgRho).^2);
    resultArray(round,5) = avgRho;
    resultArray(round,6) = sqrt(sdRho);

    % Lambda: Mean and Variance
    avgLam = sum(times(lamArray(1,:), lamArray(2,:)));
    sdLam  = sum(lamArray(2,:).*(lamArray(1,:)-avgLam).^2);
    resultArray(round,7) = avgLam;
    resultArray(round,8) = sqrt(sdLam);

    % Mu: Mean and Variance
    avgMu = sum(times(muArray(1,:), muArray(2,:)));
    sdMu  = sum(muArray(2,:).*(muArray(1,:)-avgMu).^2);
    resultArray(round,9)  = avgMu;
    resultArray(round,10) = sqrt(sdMu);

    % Display the portion of the resultArray that has been updated
    round;
    sub_id;
    resultArray(round,:);

    endTime(round) = toc(t);
    if ismember(round, [numQ])
    disp(['Subject  ' num2str(sub_id)  '  Round:  ' num2str(round) ...
          ' Elapsed Time:  ' num2str(endTime(round))])
    end
end

% Output is the result array and the array of priors and posteriors
result.probTensor     = tempProbMatrix;
result.posteriors.rho = rhoArray;
result.posteriors.lam = lamArray;
result.posteriors.mu  = muArray;
result.subjestimates  = resultArray;
result.priors.rho     = priors.rho;
result.priors.lam     = priors.lam;
result.priors.mu      = priors.mu;

f = result;
