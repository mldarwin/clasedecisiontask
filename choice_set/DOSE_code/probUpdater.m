% Calculates the new probabilities for each combination of lambdas and rhos
% using the users answer and the current probabilities for each
% combination
function f = probUpdater(softmaxArray,priorTensor)

%denom = 0;

% Calculate the sume of all softmax probabilities but if the answer is no,
% then we want the sum of the 1 - softmax probabilities

denomTensor = times(softmaxArray(:), priorTensor(:));
sumDenom = sum(denomTensor(:));
postProbs = denomTensor/sumDenom;

% Reshape object back into MxNxP array
%postProbs = reshape(postProbs, [numMu numLamda numRho]);


f = postProbs;
