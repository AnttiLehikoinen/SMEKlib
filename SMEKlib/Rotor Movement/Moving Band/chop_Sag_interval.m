function limits = chop_Sag_interval(angle_start, angle_end, msh)

if angle_end <= angle_start
    error('Baaaaaad interval!');
end

shiftTol = msh.bandData.shiftTol;
l_1 = ceil( (angle_start - shiftTol/2) / shiftTol ) * shiftTol + shiftTol/2;
l_2 = floor( (angle_end - shiftTol/2) / shiftTol ) * shiftTol + shiftTol/2;

N_ints = round( (l_2 - l_1) / shiftTol ) + 1;

limits = linspace(l_1, l_2, N_ints);
if abs(l_1-angle_start) > (1e-2 * shiftTol)
    limits = [angle_start limits];
end
if abs(angle_end-l_2) > (1e-2 * shiftTol)
    limits = [limits angle_end];
end


end