function [t_FA, t_FE, t_Ang, t_GEO, t_L, t_MD, tCSDFOD] = E_DTI_diff_measures_vectorized_FG_FOD(Tracts, VDims, DT, CSDFOD)

L = length(Tracts);

Lts = zeros(L,1);

for i=1:L
    Lts(i)=size(Tracts{i},1);
end

CLts = [0 ;cumsum(Lts)];

Tmat = repmat(single(0),[sum(Lts) 3]);

for i=1:L
    Tmat(CLts(i)+1:CLts(i+1),:) = Tracts{i};
end

Tmat = Tmat';

DT = E_DTI_DT_cell2mat(DT);

DT = Interpolate(single(DT), single(Tmat), single(VDims));

%%% CSDFOD

% set an SH

tCSDFOD = Interpolate(single(CSDFOD),single(Tmat), single(VDims)); %

%%% end CSDFOD

DT = DT([1 2 3 5 6 9],:);

[V, D, PD] = E_DTI_eig_Lex_H(DT);

clear DT PD;

V = [squeeze(V(1,1,:))'; squeeze(V(2,1,:))'; squeeze(V(3,1,:))'];

dummy = E_DTI_FrAn(D(1,:),D(2,:),D(3,:))*sqrt(3);

dummy(dummy>sqrt(3)) = sqrt(3);
dummy(dummy<0) = 0;

t_FA = cell(1,L);

for i=1:L
    t_FA{i} = dummy(1,CLts(i)+1:CLts(i+1))'; 
end


t_L = cell(1,L);
for i=1:L
    t_L{i} = D(:,CLts(i)+1:CLts(i+1))'; 
end

dummy = sum(D,1)/3;

t_MD = cell(1,L);

for i=1:L
    t_MD{i} = dummy(1,CLts(i)+1:CLts(i+1))'; 
end

[CL, CP, CS] = E_DTI_Westin_measures(D');

clear D;

t_GEO = cell(1,L);

for i=1:L
    t_GEO{i} = [CL(CLts(i)+1:CLts(i+1)) CP(CLts(i)+1:CLts(i+1)) CS(CLts(i)+1:CLts(i+1))] ; 
end



t_FE = E_DTI_Get_tract_dir(Tracts);
t_Ang = cell(1,L);

for i=1:L
    
    T = Tracts{i}(2:end,:)-Tracts{i}(1:end-1,:);
    T = T./repmat(sqrt(sum(T.*T,2)),[1 3]);
    
    t_FE{i} = abs(t_FE{i});
    
    t_Ang{i} = t_FA{i};
    
    t_Ang{i}(2:end-1)=(180/pi)*acos(abs(sum(T(1:end-1,:).*T(2:end,:),2)));
    t_Ang{i}(1) = (180/pi)*acos(abs(sum(T(1,:).*V(:,CLts(i)+1)',2)));
    t_Ang{i}(end) = (180/pi)*acos(abs(sum(T(end,:).*V(:,CLts(i+1))',2)));
    t_Ang{i} = real(t_Ang{i});
    
end
