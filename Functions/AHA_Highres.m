function [Mesh, info] = AHA_Highres(Mesh,info)

seg = [];

% Order of AHA segments from 0 degrees defined at inferior wall to 90
% degrees at anterior,but keeping the anterior wall as segment 1...
for j = 1:info.aha_highres(2)
    seg = [seg,((info.aha_highres(1)/2+1) + (j-1)*info.aha_highres(1)):j*info.aha_highres(1),...
        (1 + (j-1)*info.aha_highres(1)):(info.aha_highres(1)/2 + (j-1)*info.aha_highres(1))];
end

data = Mesh(info.template).Polar_Data;

r = mod(numel(data),info.aha_highres(2));
q = (numel(data) - r)/info.aha_highres(2);
dummy = ones(1,info.aha_highres(2)).*q;

counter = 0;
while r > 0
    dummy(end-counter) = q+1;
    counter = counter + 1;
    r = r - 1;
end

clear q r counter

%Determining list of slice numbers from base to apex
list = mat2cell([1:numel(data)]',dummy); clear dummy

%Transposing list
list = cellfun(@transpose,list,'UniformOutput',false);

count = 1; % Counter for segments
%%% Section defining the AHA segment in the template mesh and saving indices of those points
% Basal-1, apex-info.aha_highres(2)
for j1 = 1:info.aha_highres(2)
    
    angles = []; indices = [];

    % Extracting all angles in the respective section slices and their corresponding strains
    for j2 = list{j1}
        
        angles = [angles, data{j2}(2,:)];
        indices = [indices, data{j2}(3,:)];

    end

    % rotating aha highres segments for easy extraction of values
    % example: if we have 18 segments in theta, rotating the plot such that
    % segment 10 starts perfectly at 6 o'clock, and segment 1 at 12 o'clock
    angles = angles + pi/info.aha_highres(1);
    angles(angles>=2*pi) = angles(angles>=2*pi) - 2*pi;
    
    c = 1; dummy = 0;
    while dummy < 2*pi   
        
        % Finding indices of points in the reference mesh corresponding to the particular AHA segment
        aha{seg(count)} = indices(angles >= dummy & angles < c*(2*pi/info.aha_highres(1)));
        dummy = c*(2*pi/info.aha_highres(1));
        c = c + 1;
        count = count + 1;
    end

end
    
for j = info.timeframes
    
    % Finding the mean rsct of all points on the registered meshes within the defined AHA segments
    for j1 = 1:dot(info.aha_highres(1),info.aha_highres(2))
        
        % If there are less than 2 vertices in a segment, especially due to LVOT - then calling that segment a NaN
        if nnz(~isnan(aha{j1})) < 2
            rsct(j1) = NaN;
            if info.endo_strains
                ecc(j1) = NaN;
                ell(j1) = NaN;
            end    
        else    
            rsct(j1) = nanmean(Mesh(j).RSct_vertex(aha{j1}(~isnan(aha{j1}))));
            if info.endo_strains
                ecc(j1) = nanmean(Mesh(j).Ecc(aha{j1}(~isnan(aha{j1}))));
                ell(j1) = nanmean(Mesh(j).Ell(aha{j1}(~isnan(aha{j1}))));
            end                
        end    
    end    
    
    Mesh(j).AHA_Highres = rsct;
    if info.endo_strains
        Mesh(j).AHA_Highres_Ecc = ecc;
        Mesh(j).AHA_Highres_Ell = ell;
    end  
    
    clear rsct ecc ell
end