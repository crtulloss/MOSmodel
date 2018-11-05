% Caleb Rees Tulloss
% Jagannaath Shiva Letchumanan
% ELEN 6302 MOS
% Project: Simplified All-Region MOSFET Model

% Parameters
classdef parameters
    properties (Constant)
        % reasonably "assumed" parameters
        NA = 5e17;                      % assumption
        tox = 10.5e-7;                  % from website
        VFB = -1.75;                    % assumption

        % extracted parameters
        u0 = 500;

        % calculated values
        phiF = constants.phit*log(parameters.NA/1e10);
        Cox = constants.epox/parameters.tox;
        gamma = constants.sqrt2qeps*sqrt(parameters.NA)/parameters.Cox;
    end
end