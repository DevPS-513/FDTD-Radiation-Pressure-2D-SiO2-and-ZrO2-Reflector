function reset_axis(src,eventdata,h)

global names
eval(['global ' names])

xlim([cx_1-.7*Lx_temp cx_1+.7*Lx_temp])
ylim([cy_1-.7*Ly_temp cy_1+.7*Ly_temp])
zlim([cz_1-.7*Lz_temp cz_1+.7*Lz_temp])

end

