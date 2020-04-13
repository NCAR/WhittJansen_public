function smth=smooth(in,n)
smth=filtfilt(ones(n,1)./n,1,in);
end