clc
clear all
close all
uiopen('E:\University\project\Skin Modeling\blood vessel in skin\1-Derivative Gaussian\1-dispersive\vessel_in_three_layer_dispersive_skin_fdtd_cpml_gaussian_30G_3D_3D.fig',1)
for i=0:64
    saveas(gcf,['vessel_in_three_layer_non_dispersive_skin_fdtd_cpml_gaussian_3D_3D_' num2str(65-i) '.tiff'])
    close(gcf)
end