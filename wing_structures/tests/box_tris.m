function T = box_tris(lo, hi)
% BOX_TRIS  Axis-aligned box as 12 triangles with outward normals.
x = [lo(1) hi(1)]; y = [lo(2) hi(2)]; z = [lo(3) hi(3)];
P = [x(1) y(1) z(1); x(2) y(1) z(1); x(2) y(2) z(1); x(1) y(2) z(1);
     x(1) y(1) z(2); x(2) y(1) z(2); x(2) y(2) z(2); x(1) y(2) z(2)];
F = [1 3 2; 1 4 3;   5 6 7; 5 7 8;   1 2 6; 1 6 5;
     2 3 7; 2 7 6;   3 4 8; 3 8 7;   4 1 5; 4 5 8];
T = [P(F(:,1),:) P(F(:,2),:) P(F(:,3),:)];
end
