function a = aflxfn(D1,delta)
% calculate fraction of nutrient consumption that is remineralized locally
% in the top box above D1
a=delta./D1.*(1-exp(-D1./delta));
end

