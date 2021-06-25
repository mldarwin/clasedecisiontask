function f = shChoiceOptimal(subjdata, priors, priorTensor, sub_id, softmaxArray, numQ)
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
%%% Updated: 8/05/2017

t= tic;


% Loop over all questions and remove the ones that the user did not answer
%subjdata = subjdata(subjdata(:,4) ~=999,:);

% Number of questions
%numQ = size(subjdata,1);

% Load in prior distributions from "priors" argument
% Each array has two rows. The first holds the value of the parameter
% The second holds the probabilities associated with each value.

rhoArray = priors.rho ; lamArray = priors.lam ; muArray = priors.mu ;
numRhos = length(rhoArray); numLambdas = length(lamArray); numMus = length(muArray);

% Parameter probability tensor from Kronecker prod of intial prior probs is
% generated and reshaped  with Mu, Lam, Rho, along dimensions 1,2,3 respectively
%priorTensor = kron(priors.rho(2,:), kron(priors.lam(2,:), priors.mu(2,:)));
%priorTensor = reshape(priorTensor(:), [numMus numLambdas numRhos]);
%priorTensor = priorTensor;


% Create an array that will hold the questions asked, the Kullback-Leibler div
% for our question, and the average and standard error of lam, rho, and mu
% Colums are as follows :
% |1. GainQ   |2. LossQ   |3. SureQ |4. KL# |5-10. Estimates/SD |11. Choice
% |12. Round# |13. SubjID |14. qID |15. Identifier for Frydman if SHF dataset
resultArray = ones(numQ, 14);
resultArray(1:numQ,13) = repmat(sub_id,numQ,1);
subjdata = subjdata(:,1:4);

% Avoid modifying main softmax array each round by copying it and keeping
% only information relevant to current subject.
q_index = ismember(softmaxArray.master_qlist, subjdata(:,1:3), 'rows');
klSoftmaxArray = softmaxArray.array(:,:,:,q_index);
klSoftmaxQlist = softmaxArray.master_qlist(q_index,:);

% Loop over the number of questions and input results into array
%endTime = ones(numQ,1); % For recording time of each round
for round = 1:numQ
    %round
    resultArray(round, 12) = round; % Store the round number in the array

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% KL Procedure for Optimal Q Order : Calculate the combination of
    %%% x, y, and z with the largest information number
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    klSoftmaxQlist_z = klSoftmaxQlist(:,1);
    klSoftmaxQlist_x = klSoftmaxQlist(:,2);
    klSoftmaxQlist_y = klSoftmaxQlist(:,3);
    tempKLArray = ones(size(subjdata,1),5);

    %priorTensorNoMu = sum(sum(priorTensor, 1),2); % Integrate over Mu
    parfor (question = 1:size(subjdata,1),24)
		% for (question = 1:size(subjdata,1))
        % Uncomment for integrating out the Mu parameter
        %klSoft = sum(sum(klSoftmaxArray(:,:,:,question), 1),2); % Integrate over Mu
        klSoft = klSoftmaxArray(:,:,:,question);
        tempKL = KLCalculator(priorTensor, klSoft);
        z = klSoftmaxQlist_z(question);
        x = klSoftmaxQlist_x(question);
        y = klSoftmaxQlist_y(question);
        answer = subjdata(ismember(subjdata(:,1:3), [z x y], 'rows'),4);
        tempKLArray(question,:) = [z x y abs(tempKL) answer];
    end

    % check that no elements are NaN
    if sum(isnan(tempKLArray(:,4)))>0
        disp('Check KL: Some elements of KL procedure produce NaN values')
        return
    end

    % Input the optimal question,its KL number and the subjects answer into
    % the result array
    if numel(tempKLArray(tempKLArray(:,4) == max(tempKLArray(:,4)))) > 1
        % Condition for ties for maximum KL due to numerical inaccuracy
        disp 'Randomized draw from Question set due to equal KL values'
        tempKLArray = tempKLArray(tempKLArray(:,4) == max(tempKLArray(:,4)),:);
        msize = numel(tempKLArray(:,4));
        index = randi(msize);
        resultArray(round,[1:4 11]) = tempKLArray(index, :);
    else
        %tempKLArray(ismember(tempKLArray(:,4), max(tempKLArray(:,4)), 'rows'), :)
        resultArray(round,[1:4 11]) = ...
            tempKLArray(ismember(tempKLArray(:,4), max(tempKLArray(:,4)), 'rows'), :);
    end

    vals = num2cell(resultArray(round,[1:3 11]));
    [z,x,y,user_answer] = deal(vals{:});

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Softmx array is calculated from all questions used on all subjects so
    % the question index number must be obtained from a matching on the
    % master question list.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % For the question_num/id to be the same across all subjects we have to get
    % the values from the master_qlist
    question_num = find(ismember(klSoftmaxQlist, [z x y], 'rows'));
    resultArray(round, 14) = find(ismember(softmaxArray.master_qlist, ...
                                            [z x y], 'rows'));


    % Obtain the softmax for the question asked.
    % If the users answer is yes, then we use the softmax. If the user
    % answers no then we set softmax = 1 - sofmax
    if user_answer == 1
        tempsoftmaxArray = klSoftmaxArray(:,:,:, question_num);
    elseif user_answer == 0
        tempsoftmaxArray = 1 - klSoftmaxArray(:,:,:, question_num);
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


        tempStruct = remove_data(subjdata, klSoftmaxQlist, klSoftmaxArray, ...
                              question_num,[z x y user_answer]);
        subjdata = tempStruct.subjdata;
        klSoftmaxQlist = tempStruct.klSoftmaxQlist;
        klSoftmaxArray = tempStruct.klSoftmaxArray;
        continue
    end

    % Calculate the new array of probabilities for every combination of
    % rho, lambda and mu
    % P(rho,lam,mu | ans )  for all param values are now in tempProbMatrix
    tempProbMatrix = probUpdater(tempsoftmaxArray, priorTensor);
    tempProbMatrix = reshape(tempProbMatrix(:), [numMus numLambdas numRhos]);
    priorTensor    = tempProbMatrix; % transpose(tempProbMatrix(:));



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

    % Remove optimal question from the subject data and softMax structure
    % before next iteration

    tempStruct = remove_data(subjdata, klSoftmaxQlist, klSoftmaxArray, question_num,[z x y user_answer]);
    subjdata = tempStruct.subjdata;
    klSoftmaxQlist = tempStruct.klSoftmaxQlist;
    klSoftmaxArray = tempStruct.klSoftmaxArray;



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
result.time           = endTime;
f = result;

function rm = remove_data(subjdata, klSoftmaxQlist, klSoftmaxArray, ...
                          question_num, question_vec)
    subjdata(ismember(subjdata, question_vec, 'rows'), :) = [];
    klSoftmaxQlist(question_num,:) = [];
    klSoftmaxArray(:,:,:,question_num) = [];
    rm.subjdata = subjdata; rm.klSoftmaxQlist = klSoftmaxQlist;
    rm.klSoftmaxArray = klSoftmaxArray;
end
end
