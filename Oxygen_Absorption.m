function OxygenLoss = Oxygen_Absorption(F)
%%% This is the function to estimate the power loss due to the oxygen absorption.
%%% 
%%% The oxygen absorption model is specified in Section 7.6.1 in the 3GPP TR 38.901.
%%% The power loss values are obtained from the measurements given in [1].
%%% 
%%% The oxygen loss is frequency-dependent and is measured in dB per kilometer, dB/km.
% 
% Input
%           F           : carrier frequency in Hz
% Output
%           OxygenLoss  : dB/km
%
% [1] T. S. Rappaport el at.,"Millimeter wave mobile communications for 5G cellular: It will work!",2013.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F = F*1e-9; % GHz

switch(floor(F))
    
    case 20  
        OxygenLoss = 0.1; 
    case 23 
        OxygenLoss = 0.2;
    case 25
        OxygenLoss = 0.15;
    case 30 
        OxygenLoss = 0.09;
    case 43 
        OxygenLoss = 0.15;
    case 50
        OxygenLoss = 0.3;
    case 53
        OxygenLoss = 1; 
    case 54
        OxygenLoss = 2.2;
    case 55
        OxygenLoss = 4;
    case 56
        OxygenLoss = 6.6;
    case 57
        OxygenLoss = 9.7;
    case 58
        OxygenLoss = 12.6;
    case 59
        OxygenLoss = 14.6;
    case 60
        OxygenLoss = 15;
    case 61
        OxygenLoss = 14.6;
    case 62
        OxygenLoss = 14.3;
    case 63
        OxygenLoss = 10.5;
    case 64
        OxygenLoss = 6.8;
    case 65
        OxygenLoss = 3.9;
    case 66
        OxygenLoss = 1.9;
    case 67
        OxygenLoss = 1;
    case 70
        OxygenLoss = 0.5;
    case 75
        OxygenLoss = 0.4;
    case 80
        OxygenLoss = 0.3;
    case 90
        OxygenLoss = 0.35;
    otherwise 
        OxygenLoss = 0;
        
end

