function [G] = CreateGraph(Y,faces)

Npts = size(Y,1);
A = [0, 0, 0];

for i = 1:Npts
    pt = i;
    idx = find(i == faces(:,1) | i == faces(:,2) | i == faces(:,3));
    for j = 1:length(idx)
        face = faces(idx(j),:);
        oth = face(face ~= i);
        d1 = sqrt(((Y(i,1)-Y(oth(1),1))^2)+((Y(i,2)-Y(oth(1),2))^2)+((Y(i,3)-Y(oth(1),3))^2));
        d2 = sqrt(((Y(i,1)-Y(oth(2),1))^2)+((Y(i,2)-Y(oth(2),2))^2)+((Y(i,3)-Y(oth(2),3))^2));
        A = [A; ...
            i, oth(1), d1; ...
            i, oth(2), d2];
    end
end

Mat = A(2:end,:);

for i = 1:size(Mat,1)
    pair = Mat(i,1:2);
    idx1 = find(Mat(:,1) == pair(1) & Mat(:,2) == pair(2));
    idx2 = find(Mat(:,1) == pair(2) & Mat(:,2) == pair(1));
    idx = [idx1;idx2];
    idx = idx(idx~=i);
    Mat(idx,:) = NaN;
end

MatR = Mat(~isnan(Mat(:,1)),:);

G = graph(MatR(:,1),MatR(:,2),MatR(:,3));

end