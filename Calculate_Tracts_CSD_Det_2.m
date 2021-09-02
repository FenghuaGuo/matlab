function ret = Calculate_Tracts_CSD_Det_2

ret = 1<0;

h_f = findobj('Tag','MainExploreDTI');
data = get(h_f, 'userdata');

if isempty(data.CSD_FOD)
    my_msgbox('Load CSD FOD data first.','CSD FOD tracking','modal');
    return;
end

Lm = size(data.CSD_FOD,4);

if ~any(Lm == SH.lmax2n(2:2:2000))
    my_msgbox('Load *proper* CSD FOD data first.','CSD FOD tracking','modal');
    return;
end

seedPoint = CreateSeedPointList(data);
% seedPoint = CreateSeedPointList_FH_fromAtlas(data); % modified by F.guo on 2020-0110

if isempty(seedPoint)
    ret = false;
    return;
end

h_wb = my_waitbar(0,'Processing diffusion info along the tracts...');pause(0.01)

stepSize = data.tracts.params.SS;
threshold = data.HARDI.CSD.FODThresh;
maxAngle = data.tracts.params.maxAngle;
lengthRange = [data.tracts.params.minFiber data.tracts.params.maxFiber];
v2w = diag(data.DTI.VDims); v2w(4,4) = 1;


t = SHTracker_FG(v2w); % t = SHTracker_FG(v2w);
t.setData(data.CSD_FOD);
t.setParameters(stepSize, threshold, maxAngle, lengthRange); t.setProgressbar(false);

% [filename, pathname] = uigetfile({'*.nii';'*.mat'}, 'Select an ExploreDTI AFD nii-file...'); pause(0.01);
% if isa(filename,'char');
%     filename_in = fullfile(pathname, filename);
%     afd = E_DTI_read_nifti_file(filename_in);
% else
%     return;
% end
% t.setAFD(afd)

my_waitbar(1);close(h_wb);pause(0.01);

h_wb = my_waitbar(0,'Tracking...');

% [Tracts, data.tracts.TractsFOD] = t.track(seedPoint);
[Tracts, data.tracts.TractsFOD,TractsEnd,TractsCSDFOD,tracte_a,tracte_v,tract_angle,tract_dir,tract_val] = t.track(seedPoint);
% [Tracts, data.tracts.TractsFOD,data.tracts.tractCSD_FOD,data.tracts.tractDir,data.tracts.averageDir,data.tracts.stopVal,data.tracts.stopAngle] = t.track(seedPoint);

my_waitbar(1);close(h_wb);pause(0.01);

h_w = my_waitbar(0,'Processing diffusion info along the tracts...');pause(0.01)

if isempty(data.DTI.DT)

    load(data.DTI.MatfilePath,'DT');
    data.DTI.DT = DT; clear DT;

end

[data.tracts.FA, data.tracts.FE, data.tracts.Angel, data.tracts.GEO, data.tracts.Lambdas, data.tracts.MD, data.tracts.CSD_FOD] =...
    E_DTI_diff_measures_vectorized_FG_FOD(Tracts, data.DTI.VDims, data.DTI.DT, data.CSD_FOD);

num_tracts = size(Tracts,2);

data.tracts.Length = cell(1,num_tracts);

for i = 1:num_tracts
    data.tracts.Length{i} = single(size(Tracts{i},1));
end

my_waitbar(1);close(h_w);pause(0.01)

data.tracto = Tracts;
clear Tracts;

data.tracts.FList = (1:size(data.tracto,2))';
data.tracts.TractsEnd = TractsEnd;%
data.tracts.TractsCSDFOD = TractsCSDFOD;%
data.tracts.tracte_a = tracte_a;
data.tracts.tracte_v = tracte_v;
data.tracts.tract_angle = tract_angle;
data.tracts.tract_dir = tract_dir;
data.tracts.tract_val = tract_val;

set(h_f, 'userdata', []);
set(h_f, 'userdata', data);
ret = true;



% try
%
%     f_sh = zeros(size(dwi_sh));
%
%     E_DTI_open_matlabpool;
%
%     parfor i=1:size(f_sh,2)
%         f_sh(:,i) = csd.deconv(dwi_sh(:,i));
%     end
%
%     E_DTI_close_matlabpool;
%
% catch
%     f_sh = csd.deconv(dwi_sh);
% end

% clear dwi_sh;
%
% data.CSD_track_1.f_sh = f_sh;
%
% clear f_sh;
%
% set(h_f, 'userdata', []);
% set(h_f, 'userdata', data);



% load('demo_data.mat','dwi', 'grad', 'vdims');
% seedPoint = [40 41 42 43; 64 64 64 64; 30 30 30 30].*[vdims; vdims; vdims; vdims]';
% stepSize = 1; threshold = 0.1; maxAngle = 45; lengthRange = [1 1000];
% v2w = diag(vdims); v2w(4,4) = 1;
%
% bval = max(grad(:,4)); minFA = 0.8; lmax = 8;
% r_sh = SD.response(dwi, grad, bval, minFA, lmax);
% dwi = dwi(:,:,:,(grad(:,4) == bval));
% grad = grad((grad(:,4) == bval),1:3);
% sh = SH(lmax, grad);
% sd = CSD(r_sh, lmax);
% f_sh = sh.coef(dwi);
%
% t = SHTracker(v2w, sd);
% t.setData(f_sh);
% t.setParameters(stepSize, threshold, maxAngle, lengthRange); t.setProgressbar(true);
% [tract, tractVal] = t.track(seedPoint);
