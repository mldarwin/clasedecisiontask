function newRow = resample(x,n)
% created for VNI choiceset matrix

 newRow= x(randsample(size(x,1),n));

end

