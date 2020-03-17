function [ Y ] = PSD( X )
%Postive semi-definite (PSD) projection of a matrix
X(isinf(X)|isnan(X)) = 1/size(X,1);
% [VZ,DZ]             = eig((X+X')/2);
[VZ,DZ]             = eig(X);
dZ                  = real(diag(DZ));                   %eigenvalues must be real
idx                 = dZ>0;
Y                   = VZ(:,idx)*diag(dZ(idx))*VZ(:,idx)';%reconstruct
% Y                   = (Y+Y')/2;


end

