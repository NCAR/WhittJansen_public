% This function computes the maximum mixed layer depth via a matching
% condition, assuming an exponential profile with surface buoyancy b_s and
% e-folding depth h at low latitudes, and a minimum surface buoyancy
% b_winter in the north
function Dwinter=D_from_b(h,b_s,b_winter)
   Dwinter=h.*log(b_s./b_winter);
end