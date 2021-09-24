function [mask,ndir,nval,angle,inb] = include_nbtrack3 (point,mask,dir,fod,th,stepsize,VDims,val)    

                temp = ~mask;

                npoint =  point + stepsize.* dir; 
                opoint = point - stepsize.* dir;
                ffod = interpolate_neighbour(fod,point,VDims);
%                 [dir_,val_] = SHPrecomp.all_peaks(ffod,0,4);
 
                fodi = interpolate_neighbour(fod,npoint,VDims);
                [ndir, nval] = SHPrecomp.peaks(fodi,dir);
                
                fodo = interpolate_neighbour(fod,opoint,VDims);
                [odir,oval] = SHPrecomp.peaks(fodo,dir);
                angle = (180/pi)*real(acos(abs(sum(odir.*ndir,1))));
                               
                maskn = nval>th; 
                masko = oval>th;
                msk = maskn&masko;         
                mask(msk) = true;  
                inb = find(temp&mask);
                
                % afd informed tracking 
                [mask,a_dir,a_val,id] = include_high_afd_branches2(mask,dir,ffod,val);
                ndir(:,id) = a_dir(:,id);
                nval(id) = a_val(id);

%                 idk = found(temp);
%                 mask(idk) = msk;
                
                
                % make sure we don't move back in the streamline
                flipsign = sign(sum(dir.*ndir,1));               
                % update dir
                ndir = flipsign([1 1 1],:).*ndir;
end               
%              point_n6 = Get_neighbour(point);
%                nb_mat = cell2mat(point_n6);
%                nb_mat = reshape(nb_mat,[3 size(point,2) 6]);
%                nb_mat = permute(nb_mat,[1 3 2]);
%                 nb_mat = zeros(3,6,size(op_ex,2));
%                 for ii =1:6
%                     nb_mat(:,ii,:) = point_n6{1,ii};
%                 end
%                 [ndir, nval] = SHPrecomp.peaks(fod, dir);


%                 for ii = 1:size(point,2)
%                     if temp(:,ii)
%                     pb6 = squeeze(nb_mat(:,:,ii)); %3*6
                    
%                     for jj = 1:6
%                         fodi(:,jj) = fod(pb6(1,jj),pb6(2,jj),pb6(3,jj),:);
%                     end
%                     nd = repmat(dir(:,ii),[1 6]);
%                     [ndir, nval] = SHPrecomp.peaks(fodi, nd); %6x
                    
                    
%                     if sum(nval> th)>0 % at least one nb val is true
%                         maskn(:,ii) = true;
%                         dir_nb = [dir_nb ndir];
%                         val_nb = [val_nb nval];
%                         angle_nb = [angle_nb angle];
%                         
%                         ms = nval> th & angle <45; % 1-6
%                         ndir = ndir(:,ms);
%                         nval = nval(:,ms);
%                         angle = angle(:,ms);
%                         if mode == 1
%                             [~,iv] = max(nval); % ms/1-6
%                             newd(:,ii) = ndir(:,iv);%
%                             newangle(:,ii) = angle(:,iv);
%                         else
%                             [~,iv] = min(angle);
%                             newd(:,ii) = ndir(:,iv);
%                             newangle(:,ii) = angle(:,iv);
%                             
%                         end
                        
                        
%                         ap = [ap point(:,ii)];
%                         ip = [ip ii];
%                         nbv = [nbv pb6];
%                     end
%                     end
%                 end
