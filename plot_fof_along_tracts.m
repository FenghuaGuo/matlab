% function plot_fof_along_tracts
% addpath 'C:\Users\Fenghua\Downloads\code_Chantal_2019\code_Chantal_2019'
% h_f = findobj('Tag','MainExploreDTI');
% data = get(h_f, 'userdata');
% nl = E_DTI_n2lmax(size(data.CSD_FOD,4));
% nd = data.HARDI_GL.Dirs_render.the_dirs;
% sh_object = SH(nl,nd);

% 'D:\Matlab_files\100307\100307_tracts_t1'
% load('Tracts_th01_and.mat','TractsCSDFOD');
% data.CSD_FOD = E_DTI_read_nifti_file('100307_DWIb3000_MD_C_native_CSD_FOD.nii');
% data.HARDI_GL.Dirs_render.the_dirs = textread('Grad_dirs_0512.txt');

sh_object = SH(E_DTI_n2lmax(size(data.CSD_FOD,4)),data.HARDI_GL.Dirs_render.the_dirs);
vt = TractsCSDFOD{1,1};
SHPrecomp.init(8)
[dir_ini,val_ini] = SHPrecomp.all_peaks(vt',0,5);
angle = (180/pi)*real(acos(abs(sum(dir_ini{1,end-2}(:,1).*dir_ini{1,end-1}(:,1),1))));

fof = sh_object.amp(vt');
nd = data.HARDI_GL.Dirs_render.the_dirs;
figure;plot_odf_Chantal([fof(:,1);fof(:,1)],[nd;-nd])


% vt = TractsCSDFOD{1,2};
% fof = sh_object.amp(vt');

% [dir_ini,val_ini] = SHPrecomp.all_peaks(vt',0,5);
% angle = (180/pi)*real(acos(abs(sum(dir_ini{1,1}(:,1).*dir_ini{1,2}(:,1),1))));
% angle = (180/pi)*real(acos(abs(sum(dir_ini{1,1}(:,1).*dir_ini{1,2}(:,2),1))));
% figure;plot_odf_Chantal(fof(:,1),data.HARDI_GL.Dirs_render.the_dirs)
% figure;plot_odf_Chantal(fof(:,2),data.HARDI_GL.Dirs_render.the_dirs)


