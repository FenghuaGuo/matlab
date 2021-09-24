 function [tract, tractVal,tracte,tractFOD,tracte_a,tracte_v,tract_angle,tract_dir,tract_val] = onepoint_trackDir(point)
 
 % [ix,iy,iz] = ind2sub(size(fod(:,:,:,1)),idk); %[71 71 5]
 % setData
 
h_f = findobj('Tag','MainExploreDTI');
data = get(h_f, 'userdata');

stepSize = data.tracts.params.SS;
threshold = data.HARDI.CSD.FODThresh;
maxAngle = data.tracts.params.maxAngle;
lengthRange = [data.tracts.params.minFiber data.tracts.params.maxFiber];
v2w = diag(data.DTI.VDims); v2w(4,4) = 1;

t = SHTracker_FG(v2w); 
t.setData(data.CSD_FOD);
t.setParameters(stepSize, threshold, maxAngle, lengthRange); 

% tracking from the point, dir and -dir

fod = data.CSD_FOD; 
VDims = data.DTI.VDims;
point = single(point);
fod_itp = interpolate_neighbour(fod,point,VDims);                    
[point, dir, val, idk] = getInitDir(point,fod_itp,threshold);

point = point(:,1);
dir = dir(:,1);

fod_itp = interpolate_neighbour(fod,point,VDims); 
[tract1, tractVal1,tractCSD_FOD1,tracte1,tracte_a1,tracte_v1,tract_angle1,tract_dir1,tract_val1] = trackOneDir(point, dir,fod,lengthRange,stepSize,threshold,VDims,maxAngle); 
[tract2, tractVal2,tractCSD_FOD2,tracte2,tracte_a2,tracte_v2,tract_angle2,tract_dir2,tract_val2] = trackOneDir(point, -dir, fod,lengthRange,stepSize,threshold,VDims,maxAngle); 


    % join tracts from both directions
    tractCSD_FOD2 = cellfun(@flipud,tractCSD_FOD2,'UniformOutput',false);
    tractFOD = cell(size(tractCSD_FOD2));

    tract2 = cellfun(@flipud,tract2,'UniformOutput',false);
    tractVal2 = cellfun(@flipud,tractVal2,'UniformOutput',false);
    tract = cell(size(tract2));
    tractVal = cell(size(tractVal2));

    tract_angle2 = cellfun(@flipud,tract_angle2,'UniformOutput',false);
    tract_angle = cell(size(tract_angle2));
    tract_val2 = cellfun(@flipud,tract_val2,'UniformOutput',false);
    tract_val = cell(size(tract_val2));
    tract_dir2 = cellfun(@flipud,tract_dir2,'UniformOutput',false);
    tract_dir = cell(size(tract_dir2));

    for j = 1:size(tract2,2)
        if ~isempty(tract2{j})
            tract{j} = [tract2{j}; point(:,j)'];
            tractVal{j} = [tractVal2{j}; val(:,j)'];
            tractFOD{j} = [tractCSD_FOD2{j}; fod_itp(:,j)'];
            tract_val{j} = [tract_val2{j};val(:,j)'];
            tract_dir{j} = [tract_dir2{j};dir(:,j)'];
            tract_angle{j} = [tract_angle2{j};0]; 

            if ~isempty(tract1{j})
                tract{j} = [tract{j}; tract1{j}];
                tractVal{j} = [tractVal{j}; tractVal1{j}];
                tractFOD{j} = [tractFOD{j}; tractCSD_FOD1{j}];
                tract_val{j} = [tract_val{j}; tract_val1{j}];
                tract_dir{j} = [tract_dir{j}; tract_dir1{j}];
                tract_angle{j} = [tract_angle{j}; tract_angle1{j}];
            end
        else
            if ~isempty(tract1{j})
                tract{j} = [point(:,j)'; tract1{j}];
                tractVal{j} = [val(:,j)'; tractVal1{j}];
                tractFOD{j} = [fod_itp(:,j)'; tractCSD_FOD1{j}];
                tract_val{j} = [val(:,j)'; tract_val1{j}];
                tract_dir{j} = [dir(:,j)'; tract_dir1{j}];
                tract_angle{j} = [0; tract_angle1{j}];
            end
        end
    end

     fode1 =  interpolate_neighbour(fod,single(tracte1),VDims); 
     fode2 = interpolate_neighbour(fod,single(tracte2),VDims);
     for j = 1:size(tract2,2)
         tractFOD{j} = [fode2(:,j)';tractFOD{j}; fode1(:,j)'];
     end

    % enforce length limitations
    mask = cellfun('size',tract,1)>1;
    tract = tract(mask);
    tractVal = tractVal(mask);
    tractFOD = tractFOD(mask);
    tract_val = tract_val(mask);
    tract_dir = tract_dir(mask);
    tract_angle = tract_angle(mask);

    tracte = [tracte2;tracte1];
    tracte = tracte(:,mask);

    tracte_a = [tracte_a2;tracte_a1];
    tracte_a = tracte_a(:,mask);
    tracte_v = [tracte_v2;tracte_v1];
    tracte_v = tracte_v(:,mask);

end

    

function [tract, tractVal,tractCSD_FOD,tracte,tracte_a,tracte_v,tract_angle,tract_dir,tract_val] = trackOneDir(point, dir,fod, lengthRange,stepSize,threshold,VDims,maxAngle)
            tract = cell(1,size(point,2));
            tractVal = cell(1,size(point,2));
            flist = 1:size(point,2);
            flist_end = 1:size(point,2);
            
            tracte_v = nan(1,size(point,2));
            tracte_a = nan(1,size(point,2));
            tracte = nan(3,size(point,2));
            tractCSD_FOD = cell(1,size(point,2)); 
            tract_angle = cell(1,size(point,2));
            tract_dir = cell(1,size(point,2));
            tract_val = cell(1,size(point,2));
                       
            for it = 1:(lengthRange(2)/stepSize)

                point = point + stepSize .* dir; 
                fun = interpolate_neighbour(fod,point,VDims);
                [newDir, val, angle] = getDir(fun,dir);
                flipsign = sign(sum(dir.*newDir,1));               
                newDir = flipsign([1 1 1],:).*newDir;
                mask = val > 0.5*threshold; 
                
                % include neighbour voxels
                [mask,ndir,nval,newangle,idk] = include_nbtrack3 (point,mask,dir,fod,threshold,stepSize,VDims,val);
                newDir(:,idk) = ndir(:,idk);
                angle(idk) = newangle(idk); 
                val(idk) = nval(idk); 
                

                 
                % continue
                temp = ~mask;
                tracte(:,flist_end(temp)) = point(:,temp);
                tracte_a(:,flist_end(temp)) = angle(:,temp);
                tracte_v(:,flist_end(temp)) = val(:,temp);
                flist_end = flist_end(mask);
                point = point(:,mask);           
                flist = flist(mask);
                dir = dir(:,mask);
                newDir = newDir(:,mask);
                angle = angle(mask); 
                val = val(:,mask);

                % mask out large angles, a2              
                mask = angle < maxAngle;                     
                temp = ~mask;
                tracte(:,flist_end(temp)) = point(:,temp);
                tracte_a(:,flist_end(temp)) = angle(:,temp);
                tracte_v(:,flist_end(temp)) = val(:,temp);
                flist_end = flist_end(mask);
                point = point(:,mask);
                dir = dir(:,mask);
                newDir = newDir(:,mask);
                flist = flist(mask);
                val = val(:,mask);
                angle = angle(:,mask);

                flipsign = sign(sum(dir.*newDir,1));
                dir = flipsign([1 1 1],:).*newDir;
                
                % stop if we are out of points
                if isempty(point)
                    break
                end

                % add points to the tracts
                 for i=1:length(flist)
                    tract{flist(i)}(it,:) = point(:,i);
                    tractVal{flist(i)}(it,:) = val(:,i);
                    tractCSD_FOD{flist(i)}(it,:) = fun(:,i);
                    tract_angle{flist(i)}(it,:) = angle(:,i);
                    tract_dir{flist(i)}(it,:) = dir(:,i);
                    tract_val{flist(i)}(it,:) = val(:,i);
            
                end

            end

end
        

        function [point, dir, val, idx] = getInitDir(point,fun,threshold)

            [dir, val] = SHPrecomp.all_peaks(fun, threshold, 4);
            c = cellfun('size', dir, 2);
            idx = [];
            for i = 1:size(c,2)
               idx = [idx ones(1,c(i))*i];
            end
            point = point(:,idx);
            
            dir = cellfun(@single,dir,'UniformOutput',false);
            val = cellfun(@single,val,'UniformOutput',false);
            
            dir = cell2mat(dir);
            val = cell2mat(val);
        end
        
        function [dir, val, angle] = getDir(fun,prevDir)

            [dir, val] = SHPrecomp.peaks(fun, prevDir);
            angle = (180/pi)*real(acos(abs(sum(prevDir.*dir,1))));
        end  