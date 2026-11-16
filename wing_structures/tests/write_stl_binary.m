function write_stl_binary(fname, T)
% WRITE_STL_BINARY  Minimal binary STL writer (test helper).
n = size(T,1);
fid = fopen(fname, 'w');
fwrite(fid, zeros(1,80), 'uint8');
fwrite(fid, n, 'uint32');
for i = 1:n
    nrm = cross(T(i,4:6)-T(i,1:3), T(i,7:9)-T(i,1:3)); nrm = nrm/norm(nrm);
    fwrite(fid, [nrm T(i,:)], 'single');
    fwrite(fid, 0, 'uint16');
end
fclose(fid);
end
