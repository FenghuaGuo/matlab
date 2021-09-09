% adapted from E_DTI_PlotHARDI_blobs, fguo, use TractCSDFOD,TractMask, notROINr
% Calculate_Tracts_CSD_Det_2
% TractCSDFOD = E_DTI_read_nifti_file('100307_DWIb3000_MD_C_native_CSD_FOD.nii');

% to use random input Tracts.mat, without/with Tractend points.
% former:plot_HARDI_along_tract2, and plot_tract_along_with_glyphs
% t_vox = cellfun(@(x) x/1.25,Tracts,'UniformOutput',false);

file = 'D:\Matlab_files\100307\100307_tracts_t1\Tracts_th01_and_test2.mat'; %Default,depend on seed point selection
%file = 'Tracts_ROI5_afd1tr7.mat';
file_info = who('-file',file);
if ismember('TractsEnd',file_info)
    load(file, 'TractsCSDFOD','Tracts','TractsEnd')
else
    load(file,'Tracts')
end

h_f = findobj('Tag','MainExploreDTI');
data = get(h_f, 'userdata');
sh_object = SH(E_DTI_n2lmax(size(data.CSD_FOD,4)),data.HARDI_GL.Dirs_render.the_dirs);

vertices = data.HARDI_GL.Dirs_render.the_dirs;
faces = convhulln(double(vertices));
sv = size(vertices);
sf = size(faces);

h_sc = 1;
spco = data.tracts.fancy.spco; % Specular Strength [0 1]
spexp = data.tracts.fancy.spexp; % Specular Exponent [1 50]
Amb_str = data.tracts.fancy.Amb_str; % Ambient Strength [0 1]
Diff_str = data.tracts.fancy.Diff_str; % Diffuse Strength [0 1]

the_mask = ones(174,145,145);
count = 1;
% Tract points/voxels from .mat Tracts, or data.tracto
for i = 1:size(TractsCSDFOD,2) %1*nTracst cell, voxels+2end x sh_coef

%     [x y z] = E_DTI_Get_ROI_entries(i,the_mask); %data.HARDI_GL.Use_mask.the_mask);
    [x y z] = E_DTI_Get_tract_entries2(TractsCSDFOD{1,i},the_mask, Tracts{1,i},TractsEnd(:,i));
    
    if ~isempty(x)
        L = length(x);
        Vert = [x y z].*repmat(data.DTI.VDims,[L,1]);

        VERTICES = repmat(single(0),[sv(1)*size(Vert,1) 3]);
        COL = VERTICES;
        FACES = repmat(single(0),[sf(1)*size(Vert,1) sf(2)]);

%         h_wb = my_waitbar(0,'Drawing HARDI glyphs...');
                
        for k = 1:size(Vert,1)
            FACES((k-1)*sf(1)+1:k*sf(1),:) = faces + (k-1)*sv(1);
                        
%                 sh_ini = data.CSD_FOD(x(k),y(k),z(k),:);
                sh_ini = data.CSD_FOD(round(x(k)),round(y(k)),round(z(k)),:);
                sh_ini = sh_ini(:);

            f = sh_object.amp(sh_ini);

            f(f<0)=0;

            if h_sc==1
                f = 0.5*(f/max(f(:))).*data.HARDI_GL.Scale;
            else
                f = 2*f.*data.HARDI_GL.Scale;
            end
            
            if all(f==0) || any(isnan(f))
                f(:)=sqrt(eps);
                disp('Oops...')
            end
            
            f(f<0.01)=0.001;

            VERTICES((k-1)*sv(1)+1:k*sv(1),1) = min(data.DTI.VDims)*vertices(:,1).*f+Vert(k,1);
            VERTICES((k-1)*sv(1)+1:k*sv(1),2) = min(data.DTI.VDims)*vertices(:,2).*f+Vert(k,2);
            VERTICES((k-1)*sv(1)+1:k*sv(1),3) = min(data.DTI.VDims)*vertices(:,3).*f+Vert(k,3);


            COL((k-1)*sv(1)+1:k*sv(1),:)=abs(vertices(:,[2 1 3]));

%             my_waitbar(k/size(Vert,1));

        end
        drawnow;
%         my_waitbar(1);close(h_wb);drawnow;pause(0.01);

        if all(isnan(COL(:)))
            my_msgbox('No HARDI glyphs found (due to mask).','HARDI glyphs...','modal')
            continue;
        end

%         h_wb = my_waitbar(0,'Rendering HARDI glyphs...');

        set(0, 'currentfigure', h_f);
        main_axes = findobj('Tag','Main Axes');
        set(h_f, 'currentaxes', main_axes);

        if count>1
        delete(data.hDTI_tract{i-1})
%         delete(data.hTracts_FOD{i-1})
        end
        
        data.hDTI_tract{i} = patch('faces', FACES, 'vertices', VERTICES, 'FaceVertexCData', COL,...
            'FaceColor','interp', 'FaceLighting', 'phong','linestyle','none',...
            'BackFaceLighting','unlit','SpecularStrength',spco,'FaceAlpha',data.HARDI_GL.Transp,...
            'AmbientStrength',Amb_str,'DiffuseStrength',Diff_str,'SpecularExponent',spexp,'Tag','diff_glyph_object');
%          pause(5);
%        pause;

%         data.hTracts_FOD{i} = plot_tract_along_with_glyphs(file,i,data);
        data.descrGlyph_t{i} = 'HARDI';
        drawnow;
%         my_waitbar(1);close(h_wb);drawnow;pause(0.01);
%         tr_obj = findobj('Tag','diff_glyph_object');
%         delete(tr_obj(1))
%         tr_obj = findobj('Tag','Tracts_diff');

    end
       count  = count+1;
end

        delete(data.hDTI_tract{i})
%         delete(data.hTracts_FOD{i})
        
set(h_f, 'userdata', []);
set(h_f, 'userdata', data);