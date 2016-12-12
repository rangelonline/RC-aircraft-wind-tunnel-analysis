%{
Function to calculate aerodynamic center. Returns the aerodynamic center
location normalized by the aircraft reference chord. 

REFERENCE:
"Derivation of Aerodynamic Center Location from Wind Tunnel Data," Yogesh 
Babbar.

VARIABLES:
mystruct:           a wind tunnel test data structure.
xb:                 distance from aircraft nose to balance moment center. [ft]
cref:               aircraft reference chord.
cutoff:             the index of the last wind tunnel test data point
                    before separation begins. Usually found by visual
                    inspection; where the CL vs alpha curve is no longer
                    linear. 

By John Rangel 11/29/16
%}

function S = aerodynamic_center(mystruct, xb, cref, cutoff)

% Check if mystruct has necessary fields:
    if isfield(mystruct, 'alpha') == 0
        fprintf('ERROR: passed structure does not contain required WT data field (alpha). \n')
        S=0.0;
    end
        
    if isfield(mystruct, 'Cm') == 0
        fprintf('ERROR: passed structure does not contain required WT data field (Cm). \n')
        S=0.0;
    end
    
    if isfield(mystruct, 'CL') == 0
        fprintf('ERROR: passed structure does not contain required WT data field (CL). \n')
        S=0.0;
    end        

% ignore separation data points:
len = length(mystruct.alpha);
cutoff_point = len-cutoff;
alpha = mystruct.alpha(1:cutoff_point);
Cm = mystruct.Cm(1:cutoff_point);
CL = mystruct.CL(1:cutoff_point);

% calculate derivatives:
dCm_dalpha = polyfit(alpha, Cm, 1);
dCm_dalpha = dCm_dalpha(1);
dCL_dalpha = polyfit(alpha, CL, 1);
dCL_dalpha = dCL_dalpha(1);

% calculate aerodynamic center:
xb_bar = xb/cref;
xac_bar = xb_bar - dCm_dalpha/dCL_dalpha;

% return aerodynamic center:
S = xac_bar;

end %end function


