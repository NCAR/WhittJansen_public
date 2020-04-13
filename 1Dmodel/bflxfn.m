function b = bflxfn(D1,D2,delta)
% calculate fraction of nutrient consumption in the top box that is remineralized in the pycnocline box 
% between -(D2+D1) and -(D1)
%b=aflxfn(D1,delta)-(delta./D1.*(exp(-(D2./delta))-exp(-(D1+D2)./delta)));
b=aflxfn(D1,delta).*(1-exp(-(D2./delta)));
end

