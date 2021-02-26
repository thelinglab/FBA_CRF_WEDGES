function minSOA = calcMinSOA(onsets)
% input: a vector of onset times
% output: minimum SOA in sequence

onsets = sort(onsets); % sort in case the input is not already sorted
nOnsets = length(onsets);
soa = [];
for i = 2:nOnsets
    soa = [soa onsets(i)-onsets(i-1)];
end
minSOA = min(soa);
