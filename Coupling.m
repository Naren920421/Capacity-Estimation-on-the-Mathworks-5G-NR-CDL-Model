function Angles_ray = Coupling(Angles_ray_temp,M,N)
% Randomly coupling depature and arrival subpath angles within one cluste
%
% Input
%           Angles_ray_temp : M x N
%           M               : number of subpaths
%           N               : number of clusters
% Output
%           Angles_ray      : M x N
%

[~,order] = sort(rand(M,N),1);
index = order+repmat([1:M:M*N],M,1)-1;
Angles_ray_temp = Angles_ray_temp(index);
Angles_ray = reshape(Angles_ray_temp,M,N);

end

