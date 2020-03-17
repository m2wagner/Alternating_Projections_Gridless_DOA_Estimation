function [ DoAs ] = irregular_rootMUSIC( Y,r,K )
%Mark Wagner, Dec 31 2019
% Root MUSIC for ULAs using Vandermonde decomposition. 

% K         = number of sources
% Y         = MxL matrix, M = # sensors, L = # snapshots
% DoAs      = estimated DOAs in degrees
r           = r(:);
M           = length(r);
if diag(squareform(pdist(r)),1) == ones(M-1,1)
    DoAs    = rootMUSIC(Y,K);
    return
else

[t_est,~]           = IVD(r,K,Y*Y');
DoAs                = asin(sort(t_est)*2)/pi*180;

end





