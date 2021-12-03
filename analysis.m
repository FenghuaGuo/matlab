sh = SH(8, dirs);
csd = CSD(r_sh, 8);
dwi_s1_sh = sh.coef(sig1);
D1 = zeros(size(dwi_s1_sh), 'single');
parfor i=1:size(D1, 2)
D1(:,i) = csd.deconv(dwi_s1_sh(:,i));
end

XX = ndirez(752);
sh_obj = SH(8,XX);
D1_amp = sh_obj.amp(D1);
figure;plot_odf_Chantal(D1_amp(:,1),XX)
title('EDTI')
sh_new = SH_newBasis(8,XX);
D1_new = sh_new.amp(D1);
figure;plot_odf_Chantal(D1_new(:,1),XX)

% system([lib_ ' ' path_ 'mrconvert dev.nii.gz dev.mif'])
% system([lib_ ' ' path_ 'dwi2fod csd dev.mif EDTI_rf.txt fod_dwidev_rfnodev_t1.mif'])
% system([lib_ ' ' path_ 'dwi2fod csd dev.nii.gz EDTI_rf.txt fod_dwidev_rfnodev.mif -force '...
%      '-fslgrad no_dev.bvec no_dev.bval -shell '  num2str(1000*bv) ' -lmax 8 -threshold 0.1']);

system([lib_ ' ' path_ 'mrconvert fod_dwidev_rfnodev.mif fod_dwidev_rfnodev.nii -force'])
fod_mr = E_DTI_read_nifti_file('fod_dwidev_rfnodev.nii');
fod_mr = vec(fod_mr,mask);

fod_mr_map = sh_obj.amp(fod_mr);
figure;plot_odf_Chantal(fod_mr_map(:,1),XX)
title('Mrttrix th0.1')

system([lib_ ' ' path_ 'dwi2fod csd dev.nii.gz EDTI_rf.txt fod_dwidev_rfnodev_th0.mif -force '...
    '-fslgrad no_dev.bvec no_dev.bval -shell ' num2str(1000*bv) ' -lmax 8 -threshold 0'])
system([lib_ ' ' path_ 'mrconvert fod_dwidev_rfnodev_th0.mif fod_dwidev_rfnodev_th0.nii -force'])
fod_th = E_DTI_read_nifti_file('fod_dwidev_rfnodev_th0.nii');
fod_th = vec(fod_th,mask);

fod_th_map = sh_obj.amp(fod_th);
figure;plot_odf_Chantal(fod_th_map(:,1),XX)
title('mrtrix th0')

 SHPrecomp.init(8);
[p_csd_b1, v_csd_b1] = SHPrecomp.all_peaks( D1, 0, 4 );
[p_csd_b1_mr, v_csd_b1_mr] = SHPrecomp.all_peaks(fod_mr, 0, 4 );
[p_csd,v_csd] = SHPrecomp.all_peaks(fod_th,0,4);

x=0:0.1:10;

[pCSDu,pCSDu_2] = plotAngDat_2( p_csd_b1, fibredir, x ); % CSD, deviation, uncorrected RF
[pCSDc,pCSDc_2] = plotAngDat_2( p_csd_b1_mr, fibredir, x ); % CSD, deviation, corrected RF
[p1,v1] = plotAngDat_2(p_csd, fibredir,x);

figure;plot(x,pCSDu,x,pCSDc,x,p1)
legend({'EDTI','Mrtrix th0.1','Mrtrix th0'})
title('b1000 mrtrix')

% x = 0:0.1:10;
% [h_csd_unc,e] = plotAngDat2(p_csd_unc_1(1:3,:),fibredir,x); 
% h_csd_unm = plotAngDat2(p_csd_unm_1(1:3,:),fibredir,x);
% h_csd_cor = plotAngDat2(p_csd_cor_1(1:3,:),fibredir,x);
% h_dRL_unc = plotAngDat2(p_dRL_unc_1(1:3,:),fibredir,x); 
% h_dRL_unm = plotAngDat2(p_dRL_unm_1(1:3,:),fibredir,x);
% h_dRL_cor = plotAngDat2(p_dRL_cor_1(1:3,:),fibredir,x);
% m = e(1:end-1)+diff(e)/2;
% 
% figure; hold on;
% plot(m,h_dRL_unm); plot(m,h_dRL_unc); plot(m,h_dRL_cor);
% plot(m,h_csd_unm); plot(m,h_csd_unc); plot(m,h_csd_cor); % unc and unm the same
% legend({'dRL no deviation', 'dRL with deviation', 'dRL modified', 'CSD no deviation', 'CSD with deviation', 'CSD modified'})
% saveas(gcf, [nam '.fig']);



