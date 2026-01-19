function [x] = sweeplogspace(flow,fhigh)

x = round(logspace(log10(flow),log10(fhigh),8));
x(3) = x(3)-1;
x(4) = x(4)+1;
xx = nan(8,8);
xx(:,1) = x;
for kk = 1:length(x)
    xx(kk,2:end) = xx(kk,1).*[2:8];
end
xx(xx>fhigh) = 0;
xx
min(diff(sort(xx(xx>0))))
end