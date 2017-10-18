function x = PBC_pos(x,L)
%x = PBC_pos(x,L)
% gives position (x) with PBCs with length L
x = mod(x-1,L)+1;
end