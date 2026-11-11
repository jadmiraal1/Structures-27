function xz = rect_tube_midline(b_out, h_out, t)
% RECT_TUBE_MIDLINE  Midline corner points of a rectangular tube, centred
% on the origin. b_out = width (chordwise), h_out = depth (z), t = wall.
bm = b_out - t;
hm = h_out - t;
xz = [ -bm/2 -hm/2;
        bm/2 -hm/2;
        bm/2  hm/2;
       -bm/2  hm/2 ];
end
