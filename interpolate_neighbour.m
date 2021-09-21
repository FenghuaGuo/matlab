        function fun = interpolate_neighbour(fod,point,VDims)
            v2w = diag(VDims); v2w(4,4) = 1;
            point(4,:) = 1; voxel = v2w\point; voxel = voxel(1:3,:);
            fun = mex_interp(fod, voxel);

        end