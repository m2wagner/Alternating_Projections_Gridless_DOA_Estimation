function [ T ] = Tproj( A )
%Symmetric toeplitz projection of matrix A
[M,~]       = size(A);
toe         = zeros(1,M);
for     i   = 1:M
    toe(i)  = mean( [diag(A,(i-1));conj(diag(A,-(i-1)))] );
end
toe(1)      = abs(toe(1));
T           = toeplitz(toe);
end

