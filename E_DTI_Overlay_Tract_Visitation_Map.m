function E_DTI_Overlay_Tract_Visitation_Map(id)

h_f = findobj('Tag','MainExploreDTI');
data = get(h_f, 'userdata');

range = [0 100];


cmap = data.Color_map;

L = size(cmap,1);

% try
VM = data.tracts.TractMask;
VM = VM+1;
% VM(VM==0)=1;
VM = log(single(VM));


% catch
%     disp('No tracts available')
%     return;
% end

s_VM = size(VM);

if id==2
    mip_x = max((VM),[],1);
    VM = repmat(mip_x,[s_VM(1) 1 1]);
elseif id==3
    mip_y = max((VM),[],2);
    VM = repmat(mip_y,[1 s_VM(2) 1]);
elseif id==4
    mip_z = max((VM),[],3);
    VM = repmat(mip_z,[1 1 s_VM(3)]);
end

mask_VM = VM>0;

Vals = double(VM(mask_VM));

y = prctile(Vals,range);

M = y(2);
m = y(1);
Vals(Vals<m)=m;
Vals(Vals>M)=M;


if M==m;
    M = m+1;
end

Index = round((L-1)*((Vals-m)/(M-m))) + 1;


data.CV.X(mask_VM) = cmap(Index,1);
data.CV.Y(mask_VM) = cmap(Index,2);
data.CV.Z(mask_VM) = cmap(Index,3);


xField = findobj('Tag', 'XField');
yField = findobj('Tag', 'YField');
zField = findobj('Tag', 'ZField');
ix = get(xField, 'String');
iy = get(yField, 'String');
iz = get(zField, 'String');

ix = round(str2double(ix));
iy = round(str2double(iy));
iz = round(str2double(iz));


map_x = double(cat(3, squeeze(data.CV.X(ix,:,:)), squeeze(data.CV.Y(ix,:,:)), squeeze(data.CV.Z(ix,:,:))));
map_y = double(cat(3, squeeze(data.CV.X(:,iy,:)), squeeze(data.CV.Y(:,iy,:)), squeeze(data.CV.Z(:,iy,:))));
map_z = double(cat(3, squeeze(data.CV.X(:,:,iz)), squeeze(data.CV.Y(:,:,iz)), squeeze(data.CV.Z(:,:,iz))));

set(data.surfaceX,'CData',double(map_x))
set(data.surfaceY,'CData',double(map_y))
set(data.surfaceZ,'CData',double(map_z))


set(h_f, 'userdata', []);
set(h_f, 'userdata', data);