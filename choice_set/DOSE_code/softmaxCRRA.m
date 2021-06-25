function softmax = softmax(subjdata, priors, p)
% This is a script to calculate all combinations of the softmax given a
% master list of questions and an N point prior distribution of the the
% three parameters (rho, lam, mu)
%
% subjdata    = 3 column  matrix  of questions from all subjects [z x y]
% priors      = 1x1 struct with 3 fields containing 2xN matrices of parameter
%               distributions
% p           = Gamble probability
%
%      By: Judah Okwuobi
% Updated: 7/22/2017

t = tic;

%Master question list  of unique values for all subjects
master_qlist = unique(subjdata, 'rows'); % output should be sorted
qSet = sortrows(master_qlist);  % sort again just to make sure :)


rhoVals = priors.rho(1,:);
lamVals = priors.lam(1,:);
muVals  = priors.mu(1,:);
[T3,T2,T1] = ndgrid(muVals,lamVals,rhoVals);

paramSet = [T1(:) T2(:) T3(:)]; % PxQxR combinations of parameter values
qComb = kron(qSet,ones(size(paramSet,1),1));
prmComb = kron(ones(size(qSet,1),1),paramSet);
%softmaxComb = [prmComb qComb]; %PxQxRxS combinations of parameters & questions

z = qComb(:,1);
x = qComb(:,2);
y = qComb(:,3);
rho = prmComb(:,1);
lam = prmComb(:,2);
mu  = prmComb(:,3);

% tic;
% uz(z>=0) = z(z>=0).^rho(z>=0);
% uz(z<0)  = -lam(z<0).*(abs(z(z<0))).^rho(z<0);
%
% ux(x>=0) = x(x>=0).^rho(x>=0);
% ux(x<0)  = -lam(x<0).*(abs(x(x<0))).^rho(x<0);
%
% uy(y>=0) = y(y>=0).^rho(y>=0);
% uy(y<0)  = -lam(y<0).*(abs(y(y<0))).^rho(y<0);
% toc;

uz = (z.^rho).*(z>=0) +  (-lam.*(abs(z).^rho)).*(z<0);
ux = (x.^rho).*(x>=0) +  (-lam.*(abs(x).^rho)).*(x<0);
uy = (y.^rho).*(y>=0) +  (-lam.*(abs(y).^rho)).*(y<0);

U = -mu.*((1-p)*uz + p*ux - uy);

softmax = (1 + exp(U)).^(-1);

softmaxArray = reshape(softmax, [numel(muVals) numel(lamVals) numel(rhoVals) ...
                                 size(master_qlist,1)]);

% ****************************** IMPORTANT NOTE ******************************
% Mu param index is along the first dimension, lam second, rho third
% This is desired for when the array is flattened using array(:)
% rearraging the dimensions as below, produces softmax array identitcal to
% previous procedure that used nested forloops [softmax(rho_i,lam_j,mu_k,quest)]
% softmaxArray = permute(softmaxArray, [3 2 1 4]);
% In this case rho is along first dimension, lam is second, mu is last

result.array = softmaxArray;
result.master_qlist = qSet;

disp(['Softmax Elapsed Time:  ' num2str(toc(t))])
softmax = result;
