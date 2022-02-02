clc
clear all
close all
uiopen('E:\University\project\Skin Modeling\blood vessel in skin\1-Derivative Gaussian\2-non_dispersive\3-three_layer_non_dispersive_skin_fdtd_cpml_3D_gaussian\vessel_in_three_layer_non_dispersive_skin_fdtd_cpml_gaussian_3D_3D.fig',1)
for i=0:36
    saveas(gcf,['vessel_in_three_layer_non_dispersive_skin_fdtd_cpml_gaussian_3D_3D_' num2str(37-i) '.tiff'])
    close(gcf)
end