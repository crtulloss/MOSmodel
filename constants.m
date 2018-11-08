% Caleb Rees Tulloss
% Jagannaath Shiva Letchumanan
% ELEN 6302 MOS
% Project: Simplified All-Region MOSFET Model

% Physical Constants
classdef constants
    properties (Constant)
        epox = 8.854e-14 * 3.9;
        eps = 8.854e-14 * 11.9;
        k = 1.38e-23;
        q = 1.602e-19;
        roomTemp = 300;
        phit = constants.k*constants.roomTemp/constants.q;
    end
end