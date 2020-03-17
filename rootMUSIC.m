function [ DoAs ] = rootMUSIC( Y,K )
% Mark Wagner, Dec 31 2019
% Root MUSIC for ULAs using Vandermonde decomposition. 

% K         = number of sources
% Y         = MxL matrix, M = # sensors, L = # snapshots
% DoAs      = estimated DOAs in degrees


[U,~,~]     = svd(Y*Y');                    %find signal subspace
[M,~]       = size(Y);                      %number of sensors
En          = U(:,K+1:end);
G           = En*En';
p           = zeros(2*(M-1),1);
for i = -(M-1):(M-1)
    p(i+M) = sum(diag(G,i));
end
doa_phasors = roots(p);
mags        = abs(doa_phasors);

doa_phasors(abs(mags)>1) = [];
mags(abs(mags)>1) = [];
[~,inds]    = sort(mags,'descend');
DoAs        = doa_phasors(inds(1:K));
DoAs        = asin(-angle(DoAs)/pi)/pi*180;
DoAs        = sort(DoAs);

% [~,inds]    = sort(abs(mags-1),'ascend');
% DoAs        = doa_phasors(inds(1:2*K));
% DoAs        = unique(round(-angle(DoAs),4));
% DoAs        = sort(DoAs/pi*90);

end


