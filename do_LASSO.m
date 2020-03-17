function [DoAs] = do_LASSO( Y_obs, r, K, res, N )

[M,L]       = size(Y_obs);
thetas      = -90:res:90;
% N           = length(thetas);
r           = r(:);
A           = exp(1i*pi.*sind(thetas).*r);

cvx_begin
variable x_l1(length(thetas),L) complex;
cvx_precision high
cvx_quiet(true)
minimize( sum(norms(x_l1,2,2)) ) ;
subject to
norm(A*x_l1 - Y_obs,'fro') <= 1*norm(N,'fro');
cvx_end

[mu_val_LASSO,peak_LASSO] = findpeaks(norms(x_l1,2,2),'SORTSTR','descend','Npeaks', K);
DoAs    = sort(thetas(peak_LASSO));
% [RMSE_deg_LASSO,~]    = rmseCAL( theta , theta_found_LASSO , 1 )
% A_RESULTS.thetaLA{iii} = theta_found_LASSO;



end

