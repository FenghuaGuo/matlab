 classdef Tracker_FG_8 < handle
    %
    % Abstract Tractography class
    %
    % see also DTITracker, FODTracker
    %
    % Copyright Ben Jeurissen (ben.jeurissen@ua.ac.be)
    %
    properties (Access = public) %protected)
        f;
        v2w;
        stepSize;
        threshold;
        maxAngle;
        lengthRange;
        fun;
        iterations = 1;
        pb = false;
        afd;
    end
    
    methods (Access = public)
        function this = Tracker_FG_8(v2w)
            this.v2w = single(v2w);
        end 
        
        function this = setData(this, f)
            this.f = single(f);
        end
        
       function this = setAFD(this, afd)
            this.afd = single(afd);
        end
        
        function this = setParameters(this, stepSize, threshold, maxAngle, lengthRange) 
            if stepSize > 0
                this.stepSize = single(stepSize);
            else
                error('stepSize must be > 0');
            end
            if threshold >= 0
                this.threshold = single(threshold);
            else
                error('threshold must be >= 0');
            end
            if maxAngle > 0 && maxAngle <= 90
                this.maxAngle = single(maxAngle);
            else
                error('maxAngle must be in [0 90]');
            end
            if lengthRange(1) >= stepSize
                if (lengthRange(1) <= lengthRange(2))
                    this.lengthRange = single(lengthRange);
                else
                    error('lengthRange(2) must be >= lengthRange(1)');
                end
            else
                error('lengthRange(1) must be >= stepSize');
            end
        end
        
        function this = setProgressbar(this, val)
            this.pb = val;
        end
        
        function [tract, tractVal,tracte,tractFOD,tracte_a,tracte_v,tract_angle,tract_dir,tract_val] = track(this, point) %tractDir,averageDir,stopVal,stopAngle
            point = single(point);
            if (isempty(this.stepSize) || isempty(this.threshold) || isempty(this.maxAngle) || isempty (this.lengthRange))
                error('set tracking parameters first, using setParameters');
            end
            if isempty(this.f)
                error('set input data first, using setData');
            end
%             VDims = [this.v2w(1,1) this.v2w(2,2) this.v2w(3,3)];
%             file_sz = size(this.f(:,:,:,1));
%             mk = from_point_2_mask(point,file_sz,VDims');
%            E_DTI_write_nifti_file(mk,VDims,'point_mk_seed1.nii')

            % interpolate
            this.interpolate(point);
            
            % mask out NaN values
            mask_ful = ~isnan(this.fun(1,:));
            point = point(:,mask_ful);
            this.fun = this.fun(:,mask_ful);
            
            % process function
            this.process();
           
            % determine initial track direction(s)
            [point, dir, val, idk] = this.getInitDir(point);
%             [ix,iy,iz] = ind2sub(size(this.f(:,:,:,1)),idk); %[71 71 5]

% % %             [udk,ik] = unique(idk);
% % %             [apoint,ip,nb,nbv] = include_nb(old_point,udk,point);                  
            
%   v1          ibb = find(nb);
%   v1          [ia,ib] = ismember(nbv(:,(ibb-1)*2+1,:)',point','rows');
%   v1          db = dir(:,ib);   

% % %             new_nv = [];
% % %             db = []; ab = [];
% % %             for ii = 1:size(ip,2)
% % %                ind = find(nb(:,ii));
% % %                new_nv = [new_nv repmat(apoint(:,ii),[1 length(ind)])];
% % %  %  v1              [ia,ib] = ismember(nbv(:,(ind-1)*2+1,:)',point','rows');
% % %             [ia,ib] = ismember(nbv(:,(ind-1)*2+1,ii)',point','rows');
% % %             db = [db dir(:,ib)];
% % %             ab = [ab val(:,ib)];
% % %             end
            
%             % add part too extract CSD_FOD, this.fun histogram
%             nfac = this.fun(1,:);
%             voxel_th = 0.1*nfac;
                
            % mask out all small peaks
            mask = val > this.threshold;  
            this.interpolate(point);%
            fod = this.fun(:,mask);%
            
%             % add, mask &voxel_th              
%             mask_th = val > voxel_th;
%             mask = mask & mask_th;
                
            point = point(:,mask);
            dir = dir(:,mask);
% % %             point = [point new_nv];
% % %             dir = [dir db];
% % %             val = [val ab];
            
            % repeat data for probabilistic tractography
            point = repmat(point,[1 this.iterations]);
            dir = repmat(dir,[1 this.iterations]);
            val = repmat(val,[1 this.iterations]);
            
            % Track in both directions
            if this.pb; progressbar('start', [size(point,2) 1], 'Tracking in first direction'); end
            [tract1, tractVal1,tractCSD_FOD1,tracte1,tracte_a1,tracte_v1,tract_angle1,tract_dir1,tract_val1] = this.trackOneDir(point, dir); %add voxel_th, tractDir1,averageDir1,
            if this.pb; progressbar('ready'); end
            
            if this.pb; progressbar('start', [size(point,2) 1], 'Tracking in second direction'); end
            [tract2, tractVal2,tractCSD_FOD2,tracte2,tracte_a2,tracte_v2,tract_angle2,tract_dir2,tract_val2] = this.trackOneDir(point, -dir); %add voxel_th, tractDir2,averageDir2,
            if this.pb; progressbar('ready'); end
            
            % join tracts from both directions
%             tracte2 = cellfun(@flipud,tracte2,'UniformOutput',false);
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
           
%             tracte = cell(size(tracte2));
            for j = 1:size(tract2,2)
                if ~isempty(tract2{j})
                    tract{j} = [tract2{j}; point(:,j)'];
                    tractVal{j} = [tractVal2{j}; val(:,j)'];
                    tractFOD{j} = [tractCSD_FOD2{j}; fod(:,j)'];
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
                        tractFOD{j} = [fod(:,j)'; tractCSD_FOD1{j}];
                        tract_val{j} = [val(:,j)'; tract_val1{j}];
                        tract_dir{j} = [dir(:,j)'; tract_dir1{j}];
                        tract_angle{j} = [0; tract_angle1{j}];
                    end
                end
            end
            
             this.interpolate(tracte1);
             fode1 = this.fun;
             this.interpolate(tracte2);
             fode2 = this.fun;
             for j = 1:size(tract2,2)
                 tractFOD{j} = [fode2(:,j)';tractFOD{j}; fode1(:,j)'];
             end
            
            % enforce length limitations
            %             maska = cellfun('size',tract,1)*this.stepSize >= this.lengthRange(1);
            %             maskb = cellfun('size',tract,1)*this.stepSize <= this.lengthRange(2);
            %             mask = maska & maskb;
            %             tract = tract(mask);
            %             tractVal = tractVal(mask);
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
    end
    
    methods (Access = private)
        function [tract, tractVal,tractCSD_FOD,tracte,tracte_a,tracte_v,tract_angle,tract_dir,tract_val] = trackOneDir(this, point, dir) %tractDir,averageDir,
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
            
            for it = 1:(this.lengthRange(2)/this.stepSize)
                if this.pb; progressbar(size(point,2),int2str(size(point,2))); end
                % advance streamline

                point = point + this.stepSize .* dir; 

                % interpolate
                this.interpolate(point); %this.f->this.fun, point
              
                % mask out NaN values
                mask = ~isnan(this.fun(1,:)); 
                temp = ~mask;
                tracte(:,flist_end(temp)) = point(:,temp); %
                flist_end = flist_end(mask);%

                
                point = point(:,mask);
                dir = dir(:,mask);
                this.fun = this.fun(:,mask);
                flist = flist(mask);

                % process function,[adir,aval] = SHPrecomp.all_peaks(this.fun,0.1,4);
                % getDir, [point, dir, val, idx] = getInitDir(this, point)
                this.process(); 
                 % get new direction
                [newDir, val, angle] = this.getDir(dir);
               
                % make sure we don't move back in the streamline
                flipsign = sign(sum(dir.*newDir,1));               
                % update dir
                newDir = flipsign([1 1 1],:).*newDir;

                % mask out small peaks
                mask = val > this.threshold;
% ======================================================
                fod = this.f;
                th = this.threshold;
                stepsize = this.stepSize; 
                VDims = [this.v2w(1,1) this.v2w(2,2) this.v2w(3,3)];
                [mask,ndir,nval,newangle,idk] = include_nbtrack3 (point,mask,dir,fod,th,stepsize,VDims);
                newDir(:,idk) = ndir(:,idk);
                angle(idk) = newangle(idk); 
                val(idk) = nval(idk); 
                
% --------------------------------------------------------                 
                omask = mask;
                fun_ = this.fun;
                [mask,a_dir,id,a_val,ir_angle] = include_high_afd_branches2(omask,newDir,fun_,val);
                [s2_dir,s2_val] = SHPrecomp.all_peaks(this.fun,0.1,4);
                % SHPrecomp.all_peaks
                newDir(:,id) = a_dir(:,id);
                val(id) = a_val(id);
                 
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

                % mask out large angles               
                mask = angle < this.maxAngle;                     
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

                % make sure we don't move back in the streamline
                flipsign = sign(sum(dir.*newDir,1));
                
                % update dir
                dir = flipsign([1 1 1],:).*newDir;
                
                % stop if we are out of points
                if isempty(point)
                    break
                end

                % add points to the tracts
                 for i=1:length(flist)
                    tract{flist(i)}(it,:) = point(:,i);
                    tractVal{flist(i)}(it,:) = val(:,i);
%                     tractDir{flist(i)}(it,:) = theta_dir(:,i);
%                     averageDir(:,i) = nanmean(tractDir{flist(i)}); % this can be used earlier as an additional angular threshold
                    tractCSD_FOD{flist(i)}(it,:) = this.fun(:,i);
                    tract_angle{flist(i)}(it,:) = angle(:,i);
                    tract_dir{flist(i)}(it,:) = dir(:,i);
                    tract_val{flist(i)}(it,:) = val(:,i);
            
                end

            end
            % what to add between every step and whole trackOneDir function
        end
        
        function interpolate(this, point)
            point(4,:) = 1; voxel = this.v2w\point; voxel = voxel(1:3,:);
            this.fun = mex_interp(this.f, voxel);
        end
        %
        function point_n6 = Get_neighbour (this,point)
                point_n6 = cell(1,6);
                temp = point;
                temp(1,:) = point(1,:)-1;
                point_n6{1,1} = temp;
                
                temp = point;
                temp(1,:) = point(1,:)+1;
                point_n6{1,2} = temp;
                
                temp = point;
                temp(2,:) = point(2,:)-1;
                point_n6{1,3} = temp;
                
                temp = point;
                temp(2,:) = point(2,:)+1;
                point_n6{1,4} = temp;
                
                temp = point;
                temp(3,:) = point(3,:)-1;
                point_n6{1,5} = temp;
                
                temp = point;
                temp(3,:) = point(3,:)+1;
                point_n6{1,6} = temp; clear temp;
        end
    end
    methods (Access = public, Abstract = true) %protected
        process(this);       
        [point, dir, val,idk] = getInitDir(this, point);
        [dir, val, angle] = getDir(this, prevDir);
    end
end