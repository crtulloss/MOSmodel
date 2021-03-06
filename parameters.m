% Caleb Rees Tulloss
% Jagannaath Shiva Letchumanan
% ELEN 6302 MOS
% Project: Simplified All-Region MOSFET Model

% Parameters
classdef parameters
    properties (Constant)
        % reasonably "assumed" parameters -
        % used as starting assumptions
        NA = 5e17;                      % assumption
        tox = 10.5e-7;                  % from website
        VFB = -1.75;                    % assumption
        u0 = 500;
        
        % makes vel sat  and CLM negligible for extraction of other stuff
        Ec = 1e20;
        la = 0;
        VE = 1e20;

        % calculated values
        phiF = constants.phit*log(parameters.NA/1e10);
        Cox = constants.epox/parameters.tox;
        gamma = sqrt(2*constants.q*constants.eps*parameters.NA)/parameters.Cox;
    end
end