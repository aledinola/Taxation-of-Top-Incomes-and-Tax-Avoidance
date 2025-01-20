function y=quantili(x,w,q)
% Authors: Cagetti and De Nardi (2006)
% computes quantiles q of x with weights w
% w is discrete prob function
% need not be sorted or normalized
% Note: I translated this function in Fortran,
% see: https://github.com/aledinola/ToolboxFortran/blob/main/mod_numerical.f90

[xs,ix]=sort(x);
ws=w(ix);
ws=ws/sum(ws);
cums=cumsum(ws);
[xs_u,ind_u] = unique(xs,'last');
cums_u = cums(ind_u);

y=zeros(length(q),1);
for i=1:length(q)
	y(i)=interp1q(cums_u,xs_u,q(i));
end
