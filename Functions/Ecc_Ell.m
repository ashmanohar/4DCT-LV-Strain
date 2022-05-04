function Mesh = Ecc_Ell(Mesh,info)

% Copyright (C) 2021-2022 Gabrielle Colvert (gcolvert@eng.ucsd.edu)

% Create Graph
MeshG = Mesh(info.reference).rotated_verts.*info.desired_res;
MeshG(isnan(MeshG(:,3)),:) = -1000;

faces = Mesh(info.reference).faces;

G = CreateGraph(MeshG,faces);

% Find patches
clear KeepPts

Mesh1 = Mesh(info.reference).rotated_verts.*info.desired_res;

id1 = find(~isnan(Mesh1(:,3)));

b = size(id1,1);

N = round(0.15*b);

pt = unique(randsample(id1,N));

% Find closest points within radius r of each of NCalc points
r = 7.5;
[Idx,~] = rangesearch(Mesh1,Mesh1(pt,:),r);

% Keep only points along surface (i.e. avoid going across paps)
d = G.distances;

thresh = 300;

for i = 1:length(pt)
    s = i;
    t = Idx{i}(Idx{i}~=s);
    dist = d(s,t);
    KeepPts{i} = t(dist <= thresh); 
end

for i = 1:length(KeepPts)
    if length(KeepPts{i}) >= 25
        id_dum(i) = i;
    end
end
id_crit = id_dum(id_dum~=0);

pt = pt(id_crit);

NCalc = length(pt);

MeshSelect = Mesh1(pt,:);

%%% Calculating strain tensor
for j = info.timeframes
    for i = 1:NCalc
        
        XX = [Mesh1(Idx{id_crit(i)}(1),:); Mesh1(KeepPts{id_crit(i)},:)]';
        YY = [Mesh(j).rotated_verts(Idx{id_crit(i)}(1),:).*info.desired_res;...
            Mesh(j).rotated_verts(KeepPts{id_crit(i)},:).*info.desired_res]';
        
        C = [0 0 0 0 0 0 0 0 0 0 0 0 -1];
        c = [0];
        
        A0 = eye(3);
        b0 = zeros(3,1);
        
        sg0 = 1;

        [~,Aopt{j,i},~] = findAffine(XX,YY,A0,b0,sg0,C,c);

    end
    
    disp(['Done with timeframe ',num2str(j)])
    
end

for j = info.timeframes
    for i = 1:NCalc
            E{j,i} = 0.5*(Aopt{j,i}'*Aopt{j,i}-eye(3,3));
    end
end

%%% Calculating cartesian elements of strain tensor
[~,NCalc] = size(E);

for i = 1:NCalc
    Y = [MeshSelect(i,:);Mesh1(KeepPts{i},:)];
    if length(unique(Y(:,1))) <=2 || length(unique(Y(:,2))) <=2 || length(unique(Y(:,3))) <=2
        test(i) = NaN;
    else
        test(i) = 1;
    end
end

%%% Rotating strain tensors of each patch into cardiac coordinates
r = 15;
clear KeepPts2

MeshSelect = Mesh1(pt,:);

[Idx,~] = rangesearch(Mesh(info.reference).rotated_verts.*info.desired_res,MeshSelect,r);

d = G.distances;
thresh = thresh*2;

for i = 1:length(MeshSelect)
    s = i;
    t = Idx{i}(Idx{i}~=s);
    dist = d(s,t);
    KeepPts2{i} = t(dist <= thresh);
end

for i = 1:NCalc
    if ~isempty(KeepPts2{i}) && length(KeepPts2{i}) >=50
        Y = [MeshSelect(i,:);Mesh1(KeepPts2{i},:)];
        [coeff,~,~] = pca(Y);
        PNorm(:,i) = coeff(:,3);
    else
        PNorm(:,i) = [NaN NaN NaN];
    end
    
end

clear R

for i = 1:NCalc
    % normal vector to plane

    t = PNorm(:,i)';
    
    % Calculate Circ vector, assume Zc = 0, c dot t = 0, and abs(c) = 1
    xt = t(1); yt = t(2);
    yc = sqrt(1/(1+(yt/xt)^2));
    xc = -(yc*yt)/xt;
    zc = 0;
    c = [xc,yc,zc];
    
    % Calculate Long vector as cross prod bw t and c
    l = cross(t,c)/norm(cross(t,c));
    
    % Rotate all points into new coordinate system
    R{i}  = [c;l;t];
end

% rotate strain tensor
for j = info.timeframes
    for i = 1:NCalc
        E_rot{j,i} = R{i}*E{j,i}*R{i}';
    end
end

for j = info.timeframes
    for i = 1:NCalc
        Ecc(j,i) = E_rot{j,i}(1,1);
        Ell(j,i) = E_rot{j,i}(2,2);
    end
    Ecc_new(j,:) = Ecc(j,:).*test;
    Ell_new(j,:) = Ell(j,:).*test;
end

%%% Interpolating strains from patch centers to every vertex
for j = info.timeframes
    cent = nanmean(Mesh(j).rotated_verts(pt,:).*info.desired_res);
    th = atan2d((Mesh(j).rotated_verts(pt,2).*info.desired_res - cent(2)),...
        (Mesh(j).rotated_verts(pt,1).*info.desired_res - cent(1)));
    th(th < 0) = th(th < 0) + 360;
    z = Mesh(j).rotated_verts(pt,3).*info.desired_res;
    
    F_Ecc{j} = scatteredInterpolant(z(~isnan(Ecc_new(j,:))),th(~isnan(Ecc_new(j,:))),Ecc_new(j,~isnan(Ecc_new(j,:)))');
    F_Ell{j} = scatteredInterpolant(z(~isnan(Ell_new(j,:))),th(~isnan(Ell_new(j,:))),Ell_new(j,~isnan(Ell_new(j,:)))');
end

for j = info.timeframes
    cent = nanmean(Mesh(j).rotated_verts.*info.desired_res);
    All_th = atan2d((Mesh(j).rotated_verts(:,2).*info.desired_res - cent(2)),...
        (Mesh(j).rotated_verts(:,1).*info.desired_res - cent(1)));
    All_th(All_th < 0) = All_th(All_th < 0)  + 360;
    All_z = Mesh(j).rotated_verts(:,3).*info.desired_res;
    
    EccI(j,:) = F_Ecc{j}(All_z,All_th);
    EllI(j,:) = F_Ell{j}(All_z,All_th);
end

%%% Smoothing the strain map
r2 = 7.5;
for j = info.timeframes
    [Idx2,~] = rangesearch(Mesh(j).rotated_verts.*info.desired_res,Mesh(j).rotated_verts.*info.desired_res,r2);
    for i = 1:length(EllI)
        if isempty(Idx2{i})
            mu_Ell(j,i) = NaN;
            mu_Ecc(j,i) = NaN;
        else
            pts2 = Mesh(j).rotated_verts(Idx2{i},:).*info.desired_res;
            dist = sqrt(sum((pts2(1,:) - pts2).^2,2));
            psize = r2;
            mu_Ell(j,i) = sum(exp(-dist/psize^2).*EllI(j,Idx2{i})')/sum(exp(-dist/psize^2));
            mu_Ecc(j,i) = sum(exp(-dist/psize^2).*EccI(j,Idx2{i})')/sum(exp(-dist/psize^2));
        end
    end
    
    Mesh(j).Ecc = mu_Ecc(j,:);
    Mesh(j).Ell = mu_Ell(j,:);
    
end
