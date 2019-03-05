function UnitVector = SphericalUnitVector(theta,phi)

% theta:  P x 1
% phi : P x 1

UnitVector = [sind(theta.').*cosd(phi.');...
                sind(theta.').*sind(phi.');...
                cosd(theta.')]; % 3 x P

end

