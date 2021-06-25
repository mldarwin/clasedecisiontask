function f = makePriors(lowerLim,upperLim, granularity)

% Create a uniform discrete distribution with reduced probabilities on the
% end mass points for use in the DOSE procedure.
% lowerLim      = Lower limit of the prior distribution
% upperLim      = Upper limit fo the prior distribution
% granularity   = No. of mass points in the distribution
%      By : Judah Okwuobi 
% Updates : 7/21/2017

% Get uniformily spaced values and set in result matrix 
priorValues = linspace(lowerLim,upperLim,granularity); 
priorMat(1,:) = priorValues;
priorMat(2,:) = (1/granularity);
% Old version: end points of the distribution have 1/10 the 
% probability mass as the intermediate points.
% prob = (10/(10*granularity-18));
%probVec = prob*ones(granularity-2,1);
%priorMat(2,:) = [prob/10  probVec' prob/10];

%Check that probability mass sums to 1
if sum(priorMat(2,:))-1 > 1e-10
    disp 'Probability mass does not sum to 1'
    sum(priorMat(2,:))
    return
end

f = priorMat; 