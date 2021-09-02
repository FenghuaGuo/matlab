function E_DTI_Save_Tracts(pathname, filename)
pause(0.05);
pause(0.05);
h_f = findobj('Tag', 'MainExploreDTI');
data = get(h_f, 'userdata');

if isempty(data.tracts.FList)
    h = my_msgbox('There are no tracts!','Saving tracts...');
    set(h,'windowstyle','modal');
    return;
end

if nargin==0
    button =  my_questdlg('Which format for output?','Saving tract info','Matlab (*.mat)','Nifti (*.nii)','text file (*.txt)','Matlab (*.mat)');
else
    button = 'Matlab (*.mat)';
end

if isempty(button)
    return;
end


if strcmp(button,'Matlab (*.mat)')

    if nargin==0
        [filename, pathname] = uiputfile({'*.mat','Mat Files (*.mat)'},'Save tracts as ...','Tracts.mat');
    end

    exist_DTI = isa(filename,'char');

    if exist_DTI
        fname = fullfile(pathname,filename);
    else
        return;
    end

    h = my_waitbar(0,'Saving tracts...');

    VDims = data.DTI.VDims;
    FList = (1:length(data.tracts.FList'))';
    save(fname,'FList');
    clear FList;
    my_waitbar(0.1, h);
    Tracts = data.tracto(data.tracts.FList');
    save(fname,'Tracts','-append');
    clear Tracts;
    my_waitbar(0.2, h);
    TractL = data.tracts.Length(data.tracts.FList');
    save(fname,'TractL','-append');
    clear TractL;
    my_waitbar(0.3, h);
    TractFE = data.tracts.FE(data.tracts.FList');
    save(fname,'TractFE','-append');
    clear TractFE;
    my_waitbar(0.4, h);
    TractFA = data.tracts.FA(data.tracts.FList');
    save(fname,'TractFA','-append');
    clear TractFA;
    my_waitbar(0.5, h);
    TractAng = data.tracts.Angel(data.tracts.FList');
    save(fname,'TractAng','-append');
    clear TractAng;
    my_waitbar(0.6, h);
    TractGEO = data.tracts.GEO(data.tracts.FList');
    save(fname,'TractGEO','-append');
    clear TractGEO;
    my_waitbar(0.7, h);
    TractMD = data.tracts.MD(data.tracts.FList');
    save(fname,'TractMD','-append');
    clear TractMD;
    my_waitbar(0.8, h);
    TractLambdas = data.tracts.Lambdas(data.tracts.FList');
    save(fname,'TractLambdas','VDims','-append');
    clear TractLambdas;
    my_waitbar(0.9, h);
%     E_DTI_TractMask;
%     h_f = findobj('Tag', 'MainExploreDTI');
%     data = get(h_f, 'userdata');
    TractMask = data.tracts.TractMask;
    save(fname,'TractMask','-append');
    clear TractMask;

    if isfield(data.tracts,'TractMK')
        TractMK = data.tracts.TractMK(data.tracts.FList');
        save(fname,'TractMK','-append');
    end
    if isfield(data.tracts,'TractRK')
        TractRK = data.tracts.TractRK(data.tracts.FList');
        save(fname,'TractRK','-append');
    end
    if isfield(data.tracts,'TractAK')
        TractAK = data.tracts.TractAK(data.tracts.FList');
        save(fname,'TractAK','-append');
    end
    if isfield(data.tracts,'TractKA')
        TractKA = data.tracts.TractKA(data.tracts.FList');
        save(fname,'TractKA','-append');
    end
    if isfield(data.tracts,'TractsFOD')
        TractsFOD = data.tracts.TractsFOD(data.tracts.FList');
        save(fname,'TractsFOD','-append');
    end
    %
    if isfield(data.tracts,'TractsEnd')
        TractsEnd = data.tracts.TractsEnd(:,data.tracts.FList');
        save(fname,'TractsEnd','-append');
    end
    if isfield(data.tracts,'TractsCSDFOD')
        TractsCSDFOD = data.tracts.TractsCSDFOD(data.tracts.FList');
        save(fname,'TractsCSDFOD','-append');
    end

    if isfield(data.tracts,'tract_a')
        tract_a = data.tracts.tract_a(data.tracts.FList');
        save(fname,'tract_a','-append');
    end

    if isfield(data.tracts,'tract_angle')
        tract_angle = data.tracts.tract_angle(data.tracts.FList');
        save(fname,'tract_angle','-append');
    end
    if isfield(data.tracts,'tract_v')
        tract_v = data.tracts.tract_v(data.tracts.FList');
        save(fname,'tract_v','-append');
    end
    if isfield(data.tracts,'tract_dir')
        tract_dir = data.tracts.tract_dir(data.tracts.FList');
        save(fname,'tract_dir','-append');
    end
    if isfield(data.tracts,'tract_val')
        tract_val = data.tracts.tract_val(data.tracts.FList');
        save(fname,'tract_val','-append');
    end


    my_waitbar(1, h);pause(0.05);
    close(h);

elseif strcmp(button,'Nifti (*.nii)')


    [filename, pathname] = uiputfile({'*.nii','Nifti files (*.nii)'},'Save tract mask as','Tract_mask.nii');
    exist_f = isa(filename,'char');
    if exist_f
        targetfile = fullfile(pathname,filename);
    else
        return;
    end
%     E_DTI_TractMask;
%     h_f = findobj('Tag', 'MainExploreDTI');
%     data = get(h_f, 'userdata');
    TractMask = single(data.tracts.TractMask);
    VDims = data.DTI.VDims;
    E_DTI_save_nii(TractMask,VDims,targetfile)


elseif strcmp(button,'text file (*.txt)')

    pres_wd = pwd;

    dirname = uigetdir(pres_wd, 'Save *.txt files in ...');
    exist_f = isa(dirname,'char');
    if ~exist_f
        return;
    end

    [filename, pathname] = uiputfile({'*'},'Give name of txt-files (no extension)','Tracts');
    exist_f = isa(filename,'char');
    if ~exist_f
        return;
    end

    if strcmp(filename(end-3:end), '.txt')
        filename = filename(1:end-4);
    end


    A = [dirname filesep filename '_'];

    [PATHSTR,NAME,EXT] = fileparts(A);
    A = [PATHSTR filesep NAME];


    h = my_waitbar(0,'Saving tracts...');

    Tracts = data.tracto(data.tracts.FList');
    A1 = [A 'coordinates.txt'];
    fid = fopen(A1,'wt');
    for i=1:length(Tracts)
        fprintf(fid, '%12.8f %12.8f %12.8f\n', Tracts{i}');
        fprintf(fid, '\n');
    end
    fclose(fid);
    clear Tracts

    my_waitbar(0.2, h);

    TractL = data.tracts.Length(data.tracts.FList');
    A1 = [A 'lengths.txt'];
    fid = fopen(A1,'wt');
    for i=1:length(TractL)
        fprintf(fid, '%12.8f\n', TractL{i});
    end
    fclose(fid);
    clear TractL

    my_waitbar(0.4, h);

    TractFA = data.tracts.FA(data.tracts.FList');
    A1 = [A 'FA.txt'];
    fid = fopen(A1,'wt');
    for i=1:length(TractFA)
        fprintf(fid, '%12.8f\n', TractFA{i}'/sqrt(3));
        fprintf(fid, '\n');
    end
    fclose(fid);
    clear TractFA

    my_waitbar(0.6, h);

    TractAng = data.tracts.Angel(data.tracts.FList');
    A1 = [A 'angles.txt'];
    fid = fopen(A1,'wt');
    for i=1:length(TractAng)
        fprintf(fid, '%12.8f\n', TractAng{i}');
        fprintf(fid, '\n');
    end
    fclose(fid);
    clear TractAng

    my_waitbar(0.7, h);

    TractGEO = data.tracts.GEO(data.tracts.FList');
    A1 = [A 'GEO.txt'];
    fid = fopen(A1,'wt');
    for i=1:length(TractGEO)
        fprintf(fid, '%12.8f %12.8f %12.8f\n', TractGEO{i}');
        fprintf(fid, '\n');
    end
    fclose(fid);
    clear TractGEO

    my_waitbar(0.8, h);

    TractMD = data.tracts.MD(data.tracts.FList');
    A1 = [A 'MD.txt'];
    fid = fopen(A1,'wt');
    for i=1:length(TractMD)
        fprintf(fid, '%12.8f\n', TractMD{i}');
        fprintf(fid, '\n');
    end
    fclose(fid);
    clear TractMD

    my_waitbar(0.9, h);

    TractLambdas = data.tracts.Lambdas(data.tracts.FList');
    A1 = [A 'eigenvalues.txt'];
    fid = fopen(A1,'wt');
    for i=1:length(TractLambdas)
        fprintf(fid, '%12.8f %12.8f %12.8f\n', TractLambdas{i}');
        fprintf(fid, '\n');
    end
    fclose(fid);
    clear TractLambdas

    my_waitbar(1, h);pause(0.05);
    close(h);


end


set(h_f, 'userdata', []);
set(h_f, 'userdata', data);
