function[param] = create_param(name, units, val, bandwidth, npts)

param.name = name;
param.units = units;
param.val = val;
param.bandwidth = bandwidth;
param.npts = npts;
param.range = val*[1-bandwidth/2 1+bandwidth/2];
if npts > 1
    param.vals = linspace(param.range(1), param.range(2), npts);
else
    param.vals = param.val;
end
