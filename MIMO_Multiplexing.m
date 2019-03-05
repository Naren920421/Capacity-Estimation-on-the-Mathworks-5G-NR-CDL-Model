function Capacity = MIMO_Multiplexing(H,Ns,Noise_power,BW,TX_power_Microwave)
%%% This is the function to calculate channel capacity by multiplexing for microwave signals.
% Input
%           H                   : narrowband channel matrix, Nr x Nt
%           Ns                  : number of data streams, 1
%           Noise_power         : dBm/Hz
%           BW                  : Hz
%           TX_power_Microwave  : dBm
% Ouput
%           Capacity            : capacity in bits/s
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TX_power_Microwave = 10^(TX_power_Microwave/10)*1e-3; % watt
Noise_power = 10^(0.1*Noise_power)*1e-3; % watt/Hz

% Covariance matrix of transmitted symbols
CovTx = eye(Ns)*(1/Ns)*TX_power_Microwave;

% Noise power across the whole bandwidth
Rn = BW*Noise_power;

% Capacity calculation, bits/s
Capacity = BW*abs(log2(det(eye(Ns)+Rn\(H*CovTx*H'))));
