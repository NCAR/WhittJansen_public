function out=heaviside(in)
out=zeros(size(in));
out(in>0)=1;
out(in<0)=0;
out(in==0)=0.5;
end
