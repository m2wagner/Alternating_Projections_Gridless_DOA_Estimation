function [thetas] = gen_DOAs(K, tol)
%% Generate values between [0,1] with at least tol spacing between them
n   = 1e6-round(tol*1e6);
b   = ceil(tol*n);
[as,is]=sort(randperm(n-(K-1)*(b-1),K)); % if it throws an error then b/k/n are not compatible
a = as + (0:K-1)*(b-1);
thetas = (a(is)/1e6) +(tol/2);


end