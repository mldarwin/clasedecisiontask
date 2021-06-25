function f = KLCalculator(priorTensor, softmaxArray)
% Calculate the information number of a question based on the Kullback-Leibler
% divergence
%      By: Judah Okwuobi
% Updated: 8/5/2017

% Calculate the likelihoods for each combination of rho, lambda, and mu
% for if the player accepts or rejects the risk
% The multArray is the product accumulation of softmaxArrays over several rounds

if numel(priorTensor) ~= numel(softmaxArray) % Cehck for equal number of elements
    disp('The likelihood tensor and the prior tensor have mismatched dimensions')
    return
end


likelihoodArrayYes = softmaxArray;
likelihoodArrayNo  = 1-softmaxArray;

% All model joint prior probabilities are stored in vector 'pTheta'
% The ordering of the parameters is important due to how softmaxArray is created

pTheta = priorTensor(:);
lkYes = likelihoodArrayYes(:);
lkNo  = likelihoodArrayNo(:);

% Ensure the values for the likelihood arent zero
% zs = numel(lkYes(lkYes==0)) + numel(lkNo(lkNo==0));
% disp(['Number of zeroes in likelihoodArray ' num2str(zs)])
lkYes(lkYes==0) = 10^(-300);
lkNo(lkNo ==0)  = 10^(-300);

% Denominator sum with all values
sumDenomYes = sum(times(pTheta,lkYes));
sumDenomNo  = sum(times(pTheta,lkNo));

% Likelihood divided by the normalization constant
fracYes = lkYes./sumDenomYes;
fracNo  = lkNo./sumDenomNo;

EntropyYes = sum(times(times(pTheta, lkYes), log(fracYes)));
EntropyNo  = sum(times(times(pTheta, lkNo), log(fracNo)));

KLnum = EntropyYes + EntropyNo;
f = KLnum;
