function dDdt=dDdtfn_disc(Db,Dwinter,t)
dDdt=0.5.*(Dfn(Db,Dwinter,t+1)-Dfn(Db,Dwinter,t-1))./86400;
end