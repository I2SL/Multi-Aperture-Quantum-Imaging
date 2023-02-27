function z = FTZernike(r,theta,n,m)
%function z = FTZernike(r,theta,radius,n,m)
% returns the evaluation of Fourier Transform of the Zernike Polynomials
%function and its linear indices
%
% INPUTS
% n - radial index          (row vector 1xK)
% m - azimuthal index       (row vector 1xK)
% r - radial argument       (col vector Dx1)
% theta - angular argument  (col vector Dx1)
% (r,theta) are coordinate pairs, (n,m) are index pairs
%
% OUTPUTS
% z - the function output   (matrix DxK)

%z = radius*FTzRadial(radius*r,n) .* FTzAngle(theta,m);    
z = FTzRadial(r,n) .* FTzAngle(theta,m);    
end