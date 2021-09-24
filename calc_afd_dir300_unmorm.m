function [afd,msk,ndir,nval,angle] = calc_afd_dir300_unmorm(fodi,adir,aval,newDir)
% [afd1,afd2,dir2] = calc_afd_dir300_unmorm(fodi,adir,aval,ndir,nval);
% fodi 45x4-7, adir cell 1x4-7 (4-7 points in the mask), ndir 1x4-7 (current tracking direction), 4/7 = number of voxels val<th in this step of the tracking process 

% dirs = load('dodeca-3-272.txt');
dirs = load('dodeca-5-752.txt');

lmax = 8;
sh_object = SH(lmax,dirs);
sz = size(adir,2); 
msk = zeros(sz,1); 
thre = 0; %0.7;

fod_amp = sh_object.amp(fodi);   
afd = cell(1,sz);
afd_1 = zeros(1,sz);
afd_2 = zeros(1,sz);
nval = zeros(1,sz);
ndir = newDir;
angle = nan(4,sz);

for p = 1:sz % p = number of voxels

dir_peaks = adir{1,p}; 
fodd = fod_amp(:,p);
np = size(dir_peaks,2); % number of peaks>0.1 for each voxel  
for ii = 1:np
    angle(ii,p) = (180/pi)*real(acos(abs(sum(ndir(:,p).*dir_peaks(:,ii),1))));
end

if np > 1
    id = cell(1,np);
    amdd = cell(1,np);
    ii = nan(1,size(dirs,1));

    for i = 1:size(dirs,1) 

        dirs_i = repmat(dirs(i,:)',[1 np]);  
        inp = abs(dot(dir_peaks,dirs_i)); 

        if max(inp) > thre
        [~,ii(i)] = max(inp);
        end

    end

    for j = 1:np
        id{1,j} = find(ii==j);
        amdd{1,j} = fodd(id{1,j});
        afd{1,p}(j) = sum(amdd{1,j});
    end
    
    sp = afd{1,p};
    [afd_1(p),ind] = max(sp); 
    sp(ind) = [];
    afd_2(p) = max(sp);  
                          
    if afd_1(p)> 3*afd_2(p)
        msk(p) = 1;
        ndir(:,p) = dir_peaks(:,ind); 
        nval(:,p) = aval{1,p}(:,ind);
    end

elseif np == 1
        afd{1,p} = sum(fodd);
        msk(p) = 1;
        ndir(:,p) = dir_peaks(:,1);
        nval(:,p) = aval{1,p}(:,1);
else
    break
end


end
        
     
