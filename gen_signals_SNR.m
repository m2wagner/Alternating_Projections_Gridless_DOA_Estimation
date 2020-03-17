function [ Y, Y2, noise, A, A2, X] = gen_signals_SNR( L,M,K,SNR,alpha,r, src_range )
% Mark Wagner, January 10, 2019
%GENERATE GAUSSIAN RANDOM SIGNAL ARRIVING AT SENSOR ARRAY WITH ELEMENT SPACING
%GIVEN BY 'R'.
%INPUTS:
% L         SNAPSHOTS
% M         SENSORS
% K         SOURCES
% SNR       SIGNAL TO NOISE RATIO
% alpha     DIRECTIONS OF ARRIVAL AS COMPLEX VALUES ON UNIT CIRCLE
% R         SENSOR POSITIONS FROM ARBITRARY STARTING POINT


%% Generate array steering matrices and channel matrix

%% Old
r_ULA       = (0:M-1)';
A           = (alpha(:).').^r_ULA;                 %steering matrix, ULA
A2          = (alpha(:).').^r(:);                  %steering matrix, NUA
X           = repmat((src_range.^rand(K,1)),1,L).*exp(1i*2*pi*rand(K,L)); %random source signals
W           = randn(M,L)+1i*randn(M,L);                 %noise
Y           = A*X;
Y2          = A2*X;                                     %measurements

Y           = Y./norm(Y,'fro');                         %scale
Y2          = Y2./norm(Y2,'fro');
W           = (1/(10^(SNR/10)))*(W./norm(W,'fro'));
noise       = W;

end

