function [ spec ] = CBF_spec( Y, N, r)

r           = r(:);
DOAs        = linspace(-1,1,N);
Ryy         = (1/L)*Y*Y'; 
M           = size(Y,1);
spec        = zeros(length(DOAs),1);
for i = 1:N
    DOA     = DOAs(i);
    a_theta     = exp(1i*pi*DOA.*r);
    spec(i) = abs(a_theta'*Ryy*a_theta);
end
spec        = 20*log10(spec./max(spec));
 

end

