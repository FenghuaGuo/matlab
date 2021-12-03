% code in wholebraintracking_exe to build .TractMask
TractMask = repmat(0,MDims);

% for i = 1:length(FList)
%     IND = unique(sub2ind(MDims,...
%         round(double(Tracts{i}(:,1))./VDims(1)),...
%         round(double(Tracts{i}(:,2))./VDims(2)),...
%         round(double(Tracts{i}(:,3))./VDims(3))));
%     TractMask(IND) = TractMask(IND) + 1;
% end

for i = 1:length(FList)
    IND = unique(sub2ind(MDims,...
        round(double(TractsEnd(4,i))./VDims(1)),...
        round(double(TractsEnd(5,i))./VDims(2)),...
        round(double(TractsEnd(6,i))./VDims(3))));
    TractMask(IND) = TractMask(IND) + 1;
end

TractMask_b = repmat(0,MDims);
for i = 1:length(FList)
    IND = unique(sub2ind(MDims,...
        round(double(TractsEnd(1,i))./VDims(1)),...
        round(double(TractsEnd(2,i))./VDims(2)),...
        round(double(TractsEnd(3,i))./VDims(3))));
    TractMask_b(IND) = TractMask_b(IND) + 1;
end

E_DTI_write_nifti_file(TractMask,VDims,'TractMask_tractstart.nii')
E_DTI_write_nifti_file(TractMask_b,VDims,'TractMask_tractend.nii')

