function [ peak_locs ] = find_CBF_peaks( Y, r, K, sep, tol )
%Find CBF peaks with specified separation to arbitrtrary precision

[M,L]           = size(Y);
r               = real(r(:)).';
Ryy             = (1/L)*(Y*Y');


%% Evaluate CBF spectrum
samples                 = 20*(M);
sep_samples             = round(samples*sep);
spacing                 = 1/samples;
f                       = linspace(-.5-spacing,.5+spacing,samples+2); %add 2 extra points incase peak is on extrema       
CBF_spectrum            = zeros(length(f),1);
for i = 1:length(f)
    a                   = exp(1i*2*pi*f(i).*r');        %steering vector
    CBF_spectrum(i)     = 20*log10(abs(a'*Ryy*a));              %CBF spectrum
end

%% visualize
% figure(4)
% clf
% plot(f*2,CBF_spectrum)
% title('CBF spectrum')
% xlabel('$\sin(\theta)$','interpreter','latex')
% ylabel('Magnitude')
% grid on
% pause(.1)

%% Find approximate peaks
[pks,inds]              = findpeaks(CBF_spectrum,(1:samples+2),'MinPeakDistance',sep_samples); %find the peaks
if sum(ismember(inds,2)) && sum(ismember(inds,samples))
    [~,id]              = find(inds == samples);
    pks(id)             = -M;            
end
[~,id]                  = sort(pks,'descend');
peak_locs               = f(inds(id(1:min([length(id),K]))));

%% Starting from approximate root locations, zoom into function until actual root is found
peak_locs_refined       = zeros(K,1);                               %allocate space
fun                     = @(f) -20*log10( abs( diag((exp(1i*2*pi*f.*r'))' * Ryy * (exp(1i*2*pi*f.*r')) )));     %generalized null spectrum 
options                 = optimset('Display','none','TolX',tol);                      %set options
for i = 1:length(peak_locs)                                                             %for each approximate root
    peak_locs_refined(i)= fminsearch(fun,peak_locs(i)-(abs(f(1)-f(2))/2),options);      %find local minima about generalized null spectrum
end      
peak_locs               = sort((peak_locs_refined)); %sorted DoAs on [-.5,.5)




end

