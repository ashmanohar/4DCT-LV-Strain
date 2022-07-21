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
            rsct(j1) = mean(Mesh(j).RSct_vertex(aha{j1}(~isnan(aha{j1}))),'omitnan');
            if info.endo_strains
                ecc(j1) = mean(Mesh(j).Ecc(aha{j1}(~isnan(aha{j1}))),'omitnan');
                ell(j1) = mean(Mesh(j).Ell(aha{j1}(~isnan(aha{j1}))),'omitnan');
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

% Extracting all the hi-res RSct values and reshaping them into segments x
% time
hires_rsct = [Mesh.AHA_Highres];
hires_rsct = reshape(hires_rsct,[length(Mesh(1).AHA_Highres) length(Mesh)]);

% Plotting
figure('pos',[0 0 2560 1800]);
tile_plot = tiledlayout(info.aha_highres(2),info.aha_highres(1));

for j = 1:size(hires_rsct,1)
        
        nexttile
        plot(info.percent_rr,hires_rsct(j,:),'LineWidth',3); hold on;    
        ax = gca; ax.FontSize = 12; ax.FontWeight= 'bold';
        yticks([-1:0.1:1]); xticks(0:25:100);
        ylim([info.RSct_limits]);xlim([info.percent_rr(1) info.percent_rr(end)])
        grid on; grid minor
        
end

tile_plot.TileSpacing = 'compact';
tile_plot.Padding = 'normal';

annotation('textarrow',[0.02 0.02],[0.185 0.185],'string','Apex', ...
  'HeadStyle','none','LineStyle', 'none','FontSize',30,'FontWeight','bold','TextRotation',90);

annotation('textarrow',[0.02 0.02],[0.525 0.525],'string','Mid', ...
  'HeadStyle','none','LineStyle', 'none','FontSize',30,'FontWeight','bold','TextRotation',90);

annotation('textarrow',[0.02 0.02],[0.88 0.88],'string','Base', ...
  'HeadStyle','none','LineStyle', 'none','FontSize',30,'FontWeight','bold','TextRotation',90);

annotation('textarrow',[0.13 0.13],[0.95 0.95],'string','Anterior', ...
  'HeadStyle','none','LineStyle', 'none','FontSize',25,'FontWeight','bold');

annotation('textarrow',[0.297 0.297],[0.95 0.95],'string','Anteroseptal', ...
  'HeadStyle','none','LineStyle', 'none','FontSize',25,'FontWeight','bold');

annotation('textarrow',[0.448 0.448],[0.95 0.95],'string','Inferoseptal', ...
  'HeadStyle','none','LineStyle', 'none','FontSize',25,'FontWeight','bold');

annotation('textarrow',[0.595 0.595],[0.95 0.95],'string','Inferior', ...
  'HeadStyle','none','LineStyle', 'none','FontSize',25,'FontWeight','bold');

annotation('textarrow',[0.761 0.761],[0.95 0.95],'string','Inferolateral', ...
  'HeadStyle','none','LineStyle', 'none','FontSize',25,'FontWeight','bold');

annotation('textarrow',[0.92 0.92],[0.95 0.95],'string','Anterolateral', ...
  'HeadStyle','none','LineStyle', 'none','FontSize',25,'FontWeight','bold');

sgtitle([info.patient,' RS_{CT}'],'FontSize',35,'FontWeight','bold');

savefig([info.save_path,info.patient,'_AHA_HighRes_RSct.fig'])
close;

% Plotting Ecc and Ell
if info.endo_strains
    strain_name = {'Ecc','Ell'};
    
    for j = 1:2
        figure('pos',[0 0 2560 1800]);
        tile_plot = tiledlayout(info.aha_highres(2),info.aha_highres(1));
        
        eval(['hires_strain = [Mesh.AHA_Highres_',strain_name{j},'];'])
        hires_strain = reshape(hires_strain,[length(Mesh(1).AHA_Highres) length(Mesh)]);
        
        for j3 = 1:size(hires_strain,1)

            nexttile
            plot(info.percent_rr,hires_strain(j3,:),'LineWidth',3);
            ax = gca; ax.FontSize = 12; ax.FontWeight= 'bold';
            yticks([-1:0.1:1]); xticks(0:25:100);
            ylim([info.RSct_limits]);xlim([info.percent_rr(1) info.percent_rr(end)])
            grid on; grid minor

        end
        
        tile_plot.TileSpacing = 'compact';
        tile_plot.Padding = 'normal';

        annotation('textarrow',[0.02 0.02],[0.185 0.185],'string','Apex', ...
          'HeadStyle','none','LineStyle', 'none','FontSize',30,'FontWeight','bold','TextRotation',90);

        annotation('textarrow',[0.02 0.02],[0.525 0.525],'string','Mid', ...
          'HeadStyle','none','LineStyle', 'none','FontSize',30,'FontWeight','bold','TextRotation',90);

        annotation('textarrow',[0.02 0.02],[0.88 0.88],'string','Base', ...
          'HeadStyle','none','LineStyle', 'none','FontSize',30,'FontWeight','bold','TextRotation',90);

        annotation('textarrow',[0.13 0.13],[0.95 0.95],'string','Anterior', ...
          'HeadStyle','none','LineStyle', 'none','FontSize',25,'FontWeight','bold');

        annotation('textarrow',[0.297 0.297],[0.95 0.95],'string','Anteroseptal', ...
          'HeadStyle','none','LineStyle', 'none','FontSize',25,'FontWeight','bold');

        annotation('textarrow',[0.448 0.448],[0.95 0.95],'string','Inferoseptal', ...
          'HeadStyle','none','LineStyle', 'none','FontSize',25,'FontWeight','bold');

        annotation('textarrow',[0.595 0.595],[0.95 0.95],'string','Inferior', ...
          'HeadStyle','none','LineStyle', 'none','FontSize',25,'FontWeight','bold');

        annotation('textarrow',[0.761 0.761],[0.95 0.95],'string','Inferolateral', ...
          'HeadStyle','none','LineStyle', 'none','FontSize',25,'FontWeight','bold');

        annotation('textarrow',[0.92 0.92],[0.95 0.95],'string','Anterolateral', ...
          'HeadStyle','none','LineStyle', 'none','FontSize',25,'FontWeight','bold');

        sgtitle([info.patient,' ',strain_name{j}],'FontSize',35,'FontWeight','bold');

        savefig([info.save_path,info.patient,'_AHA_HighRes_',strain_name{j},'.fig'])
        close;
        
        clear hires_strain tile_plot
    end
end
    