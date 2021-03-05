function r = round2even(x)
% rounds input to nearest even integer
if mod(x,2)<1
r = fix(x);
else
r =fix(x) + 1;
end