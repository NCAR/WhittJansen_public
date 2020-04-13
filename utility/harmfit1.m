function out = harmfit1(beta,t)
out=beta(1).*sin(beta(2).*t+beta(3))+beta(4);
if abs(beta(3))>=(pi-.03) || beta(1)<0 || beta(2)<0
    out=out.*1e19;
end
end

