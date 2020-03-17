function [ root_locs, c1 ] = IVD( r, K, T )
% Mark Wagner, Dec 31, 2019
% GENERALIZED VANDERMONDE DECOMPOSITION: FINDS DOAS FROM COVARIANCE MATRIX FORMED FROM
% UNEVENLY SAMPLED SUM OF SIGNALS IN NOISE 
% INPUTS: 
% r             SENSOR LOCATIONS 
% K             NUMBER OF DOAS PRESENT
% T             COVARIANCE MATRIX
% OUTPUT:
% root_locs     DOAS GIVEN BETWEEN [-.5, .5)
% c1            SOURCE POWERS

%% Make sure inputs are correct size, find # sensors
r               = real(r(:)).';
M               = size(T,1);
P_roots         = MS_root( T, r, K, 22, M, M, 0 );
root_locs       = angle(P_roots)./(2*pi);
W_est           = exp(1i*angle(P_roots(:)).'.*r(:));%regenerate irregular Vandermonde Matrix
W_inv           = pinv(W_est);
c1              = diag(W_inv*T*W_inv');% X\Y can throw error on rank deficient matrix
% c1              = diag( ((W_est\T)*(W_est\T)') );

% figure(3)
% zplane(P_roots)
% pause(.01)
end