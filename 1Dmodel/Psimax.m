% This function computes the overturning transport (in units m^3/s)
% using the approximate form discused in the overleaf document
function Psi=Psimax(h,b_s,b_n)
   %b_n=0.0005; b_s=0.02+0.0005;
   f0=1e-4;
   D=h.*log(b_s./b_n);
   Psi=0.052*(b_s-b_n)./(f0.*(D.^(-2)+0.073.*h.^(-2)));
end