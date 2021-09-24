function [mask,dir_new,val_new,id] = include_high_afd_branches2(mask,dir,fod,val)

        dir_new = dir; 
        val_new = val;
        temp = ~mask;

        dir_ = dir(:,temp);
        fodi = fod(:,temp);
        id = find(temp);
        if isempty(id)
             return
        else                 
            [adir,aval] = SHPrecomp.all_peaks(fodi,0,4);
            [afd, mask_afd,dir2,val2,angle] = calc_afd_dir300_unmorm(fodi,adir,aval,dir_);
            mask(id) = mask_afd;
            dir_new(:,id) = dir2;
            val_new(id) = val2;
        end
end                             
%                 for ii = 1:sz
%                       
%                     if afd1(ii)> 3*afd2(ii)
%                         mask_afd(ii) = 1;
%                         ndir(ii) = dir2(ii); 
%                     end
%                 end
                
%                 % make sure we don't move back in the streamline
%                 flipsign = sign(sum(dir.*ndir,1));               
%                 % update dir
%                 ndir = flipsign([1 1 1],:).*ndir;
%                 ndir(ii) = dir{1,ii};
%                 val(ii) = val{1,ii};
%                 angle = (180/pi)*real(acos(abs(sum(odir.*ndir,1))));
                               

                


