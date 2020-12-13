function ratio = area_stat(record,seen_map)

% compute area of each sub surface
for ii = 1:8
    area_sub(ii) = sum(record == ii);
end
% compute area of each seen sub surface
for ii = 1:8
    area_sub_seen(ii) = sum(seen_map&(record == ii));
end
ratio = area_sub_seen./area_sub;
end