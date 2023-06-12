function mx = PhotonFlux_to_VisualMagnitude(phi)

% -----------------------------------------------------------------------
% Input:
% -----------------------------------------------------------------------
% mx -  visual magnitude
% -----------------------------------------------------------------------
% References :
% -----------------------------------------------------------------------
% https://en.wikipedia.org/wiki/Apparent_magnitude#Standard_reference_values
% https://lweb.cfa.harvard.edu/~dfabricant/huchra/ay145/mags.html
% -----------------------------------------------------------------------

% Physical Constants
h_planck = 6.626e-34;   % planck's constant [J s]
c = 3e8;                % speed of light [m / s]

% Reference Values
J0 = 3640;          % [Jansky] = [10^(-26) W / m^2 / Hz] Reference flux per unit spectral frequency for magnitude 0 star in spectral filter band
lambda = 550e-9;    % [m] spectral filter center wavelength
dlambda = 86e-9;    % [m] spectral filter bandwidth

% Photon Energy
nu = c/lambda;                      % center frequency of spectral filter [Hz]  
dnu = c/(lambda^2)*dlambda;         % frequency bandwidth of spectral filter [Hz]  
E_pho = h_planck*nu;                % energy of one photon at central wavelength [J]

% flux to visual magnitude
Fx0 = 10^(-26) * J0;                % reference flux per unit spectral frequency [W / m^2 / Hz] 
phi0 = Fx0 * dnu / (E_pho);         % reference flux [Photons / m^2 / s]
mx = -2.5*log10(phi / phi0);        % visual magnitude given photon flux [Photons / m^2 / s]

end

