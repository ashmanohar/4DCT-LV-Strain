function Mesh = RSCT(Mesh,info)

faces = Mesh(info.template).faces;
template = Mesh(info.reference).CPD;

%Calculating area of ED (template) faces
ab = template(faces(:,2),:) - template(faces(:,1),:);
ac = template(faces(:,3),:) - template(faces(:,1),:);
c = cross(ab,ac,2);
ED_Ar = 0.5*sqrt(sum(c.^2,2));
clear ab ac c

%Identifying faces belonging to planes
nan_ind = isnan(ED_Ar);

%Calculating Patches for RSct
for j = 1:size(faces,1)

    verts = faces(j,:);
    for k = 1:length(verts)
        [r,~] = find(faces==verts(k));
        ind{k} = r';
    end
    indxs{j} = [ind{1} ind{2} ind{3}];
    indxs{j} = unique(indxs{j},'stable');
    clear ind verts r
    
    % Removing NaN values in calculating the mean areas of patches
    temp = ED_Ar(indxs{j}); temp = temp(~isnan(temp));
    ED_Avg_Ar(j,1) = mean(temp);
    clear temp

end

% NaN-ing values of faces belonging to planes 
ED_Avg_Ar(nan_ind) = NaN;

for j = 1:length(info.timeframes)
    
    % Calculating mean areas of registered meshes
    temp_verts = Mesh(info.timeframes(j)).CPD;
    
    ab = temp_verts(faces(:,2),:) - temp_verts(faces(:,1),:);
    ac = temp_verts(faces(:,3),:) - temp_verts(faces(:,1),:);
    c = cross(ab,ac,2);
    Ar = 0.5*sqrt(sum(c.^2,2));
    
    for k = 1:size(faces,1)
        % Removing NaN values in calculating the mean areas of patches
        temp = Ar(indxs{k}); temp = temp(~isnan(temp));
        Avg_Ar(k,1) = mean(temp); clear temp
    end
    
    % NaN-ing values of faces belonging to planes 
    Avg_Ar(nan_ind) = NaN;
    
    % RSct
    Mesh(info.timeframes(j)).RSct = sqrt(Avg_Ar./ED_Avg_Ar) - 1;
    
    %RSct values at vertices
    for j4 = 1:length(template)
        [row,~] = find(Mesh(info.template).faces == j4);

        Mesh(info.timeframes(j)).RSct_vertex(j4) = nanmean(Mesh(info.timeframes(j)).RSct(row));

    end
    
    clear ab ac c Ar Avg_Ar

end    
    