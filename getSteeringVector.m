function [Ar,At] = getSteeringVector(RayAOD,RayAOA,RayZOD,RayZOA,ElePosBS,ElePosUE,WaveLength)
% Obtain the array response vectors at the BS and UE for all propagation angles.

% RayAOD/AOA/ZOD/ZOA: MN x 1
% ElePosBS/UE: 3 x Nt/Nr

% 3 x MN
UnitVectorBS = SphericalUnitVector(RayZOD,RayAOD); 
UnitVectorUE = SphericalUnitVector(RayZOA,RayAOA);

% Steering vectors at BS and UE
At = exp(1i*2*pi*UnitVectorBS.'*ElePosBS/WaveLength); % MN x Nt
At = At.'; % Nt x MN
Ar = exp(1i*2*pi*UnitVectorUE.'*ElePosUE/WaveLength); % MN x Nr
Ar = Ar.'; % Nr x MN


