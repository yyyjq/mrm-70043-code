function [SIM_M,SIM_m] = funsampleIM(IM_iM,X)

SIM_M = imresize(IM_iM,[size(X,1) size(X,2)]);
SIM_M = SIM_M-min(min(SIM_M));
SIM_M = SIM_M./max(max(SIM_M));
SIM_m = reshape(SIM_M,numel(X),1);

end

