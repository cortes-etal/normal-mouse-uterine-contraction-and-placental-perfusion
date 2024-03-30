
function [TE,Energy,IDM] = textural3D(PN)                  %output TE, Energy, and IDM
d = size(PN);                                            %calculate size of PN
TE = 0;                                                  %clear variables
Energy = 0;
IDM = 0;
for X = 1: d(1)
    for Y = 1:d(2)
        if PN(X,Y) ~= 0
            TE = TE - PN(X,Y)*log(PN(X,Y));              %sum of PN*log(PN) for each x,y          
        end %if
        Energy = Energy + PN(X,Y)^2;                     %sum of PN^2 for all x,y
        IDM = IDM + 1/(1+(X-Y)^2)*PN(X,Y);               %sum of 1/(1+X-Y)^2*PN for all x,y
    end 
end 