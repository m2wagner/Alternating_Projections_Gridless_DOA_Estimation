function [ Y, V, D ] = ITP(K, r, X )
%Extended Toeplitz Projection of matrix X

[t_est,c1]          = IVD(r, K, X );   %decompose
alpha_est           = exp(1i*2*pi*t_est);         %reconstruct alpha parameters
D                   = (diag(c1));                 %signal strengths
V                   = (alpha_est(:).').^r(:);     %reconstruct extended Vandermonde matrix
Y                   = (V*D*V');         %reconstruct Tu

% if sum(isnan(Y(:)))+ sum(isinf(Y(:))) > 0
%     Y = zeros(size(X));
% end

end

