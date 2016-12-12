%{
Function to transfer moments from one point to another.

VARIABLES:
mystruct:               A wind tunnel test data structure. 
field:                  The field within mystruct that stores the moment
                        data to be transferred. [string or char vector]
point1:                 Original point where moments are measured. 
point2:                 Desired transfer point. Assumed to be on same axis
                        as point1. 
norm:                   Normalizing length. For pitching moments the 
                        normalizing length is the aircraft reference chord 
                        (cref),and for yawing moments the normalizing 
                        length is the aircraft wingspan (b).

By John Rangel 11/29/16
%}

function S = moment_transfer(mystruct, field, point1, point2, norm)

    if strcmp(field, 'Cm') == 1
        % transfer moment coefficients to desired location:
        new_moment = mystruct.(field)+(point1-point2).*mystruct.CL./norm;

        % pass new moment:
        S = new_moment;

    elseif strcmp(field, 'Cn') == 1
        % transfer moment coefficients to desired location:
        new_moment = mystruct.(field)-(point1-point2).*mystruct.CY./norm;

        % pass new moment:
        S = new_moment;

    else
        % print error message:
        fprintf('ERROR: fieldname not recognized. Put (Cm) or (Cn) as field.\n')
        S = 1.0;
        
    end %end if
end %end function