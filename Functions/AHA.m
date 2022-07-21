function Mesh = AHA(Mesh,info)

seg = [4 5 6 1 2 3 10 11 12 7 8 9 15 16 13 14]; % Order of AHA segments from 0 degrees

name = {'Basal Anterior','Basal Anteroseptal','Basal Inferoseptal','Basal Inferior','Basal Inferolateral','Basal Anterolateral',...
    'Mid Anterior','Mid Anteroseptal','Mid Inferoseptal','Mid Inferior','Mid Inferolateral','Mid Anterolateral',...
    'Apical Anterior','Apical Septal','Apical Inferior','Apical Lateral'};

count = 1; % Counter for segments

data = Mesh(info.reference).Polar_Data;

%Defining basal, mid, and apical chunk lengths
chunks = round((numel(data) - info.lvot_limit(info.reference) + 1)/3);

%Defining basal, mid, and apical slices
list = {info.lvot_limit(info.reference):info.lvot_limit(info.reference) + chunks-1,
    info.lvot_limit(info.reference) + chunks:info.lvot_limit(info.reference) + 2*chunks - 1,
    info.lvot_limit(info.reference) + 2*chunks:numel(data)};

%%% Section defining the AHA segment in the template mesh and saving
%%% indices of those points
% Basal-1, mid-2, apex-3
for j1 = 1:3

    angles = []; indices = [];

    % Extracting all angles in the respective section slices and their corresponding strains
    for j2 = list{j1}

        angles = [angles, data{j2}(2,:)];
        indices = [indices, data{j2}(3,:)];

    end

    if j1 == 1 || j1 ==2

        % rotating aha 16 segment plot by 30 for easy extraction of values, making segments 4 and 10 start at 6 o'clock instead of 7 o'clock
        angles = angles + pi/6;
        angles(angles>=2*pi) = angles(angles>=2*pi) - 2*pi;

        c = 1; dummy = 0;
        while dummy < 2*pi
            % Finding indices of points in the reference mesh corresponding to the particular AHA segment
            aha{seg(count)} = indices(angles >= dummy & angles < c*(pi/3));
            dummy = c*(pi/3);
            c = c + 1;
            count = count + 1;
        end

    else

        % rotating aha 16 segment plot by 45 for easy extraction of values, making segment 15 start at 6 o'clock instead of 7:30 o'clock
        angles = angles + pi/4;
        angles(angles>=2*pi) = angles(angles>=2*pi) - 2*pi;

        c = 1; dummy = 0;
        while dummy < 2*pi
            % Finding indices of points in the reference mesh corresponding to the particular AHA segment
            aha{seg(count)} = indices(angles >= dummy & angles < c*(pi/2));
            dummy = c*(pi/2);
            c = c + 1;
            count = count + 1;
        end

    end

end
    
for j = info.timeframes
    
    % Finding the mean rsct of all points on the registered meshes within the defined AHA segments
    for j1 = 1:16
        rsct(j1) = mean(Mesh(j).RSct_vertex(aha{j1}(~isnan(aha{j1}))),'omitnan');
        
        if info.endo_strains
            ecc(j1) = mean(Mesh(j).Ecc(aha{j1}(~isnan(aha{j1}))),'omitnan');
            ell(j1) = mean(Mesh(j).Ell(aha{j1}(~isnan(aha{j1}))),'omitnan');
        end     
    end    
    
    Mesh(j).AHA = rsct;
    
    if info.endo_strains
        Mesh(j).AHA_Ecc = ecc;
        Mesh(j).AHA_Ell = ell;
    end    
    
    clear chunks list rsct ecc ell
end

clear aha

aha_rsct = [Mesh.AHA];
aha_rsct = reshape(aha_rsct,[length(Mesh(1).AHA) length(Mesh)]);

figure('pos',[10 10 2400 1800])
% Plotting loop
for j3 = 1:16
      
    subplot(3,6,j3)
    plot(info.percent_rr,aha_rsct(j3,:),'LineWidth',3);
    ax = gca; ax.FontSize = 20; ax.FontWeight= 'bold';
    yticks([-1:0.1:1]); xticks(0:20:100)
    ylim([info.RSct_limits]); xlim([info.percent_rr(1) info.percent_rr(end)])
    ylabel('RS_{CT}','FontSize',20); xlabel('%R-R Phase','FontSize',20)
    title([num2str(j3),': ',name{j3}],'FontSize',32)
        
end
savefig([info.save_path,info.patient,'_AHA_RSct.fig'])
close;

% Plotting Ecc and Ell
if info.endo_strains
    strain_name = {'Ecc','Ell'};
    
    for j = 1:2
        figure('pos',[10 10 2400 1800])
        
        eval(['aha_strain = [Mesh.AHA_',strain_name{j},'];'])
        aha_strain = reshape(aha_strain,[length(Mesh(1).AHA) length(Mesh)]);
        
        for j3 = 1:16

            subplot(3,6,j3)
            plot(info.percent_rr,aha_strain(j3,:),'LineWidth',3);
            ax = gca; ax.FontSize = 20; ax.FontWeight= 'bold';
            yticks([-1:0.1:1]); xticks(0:20:100)
            ylim([info.RSct_limits]); xlim([info.percent_rr(1) info.percent_rr(end)])
            ylabel(strain_name{j},'FontSize',20); xlabel('%R-R Phase','FontSize',20)
            title([num2str(j3),': ',name{j3}],'FontSize',32)

        end
        savefig([info.save_path,info.patient,'_AHA_',strain_name{j},'.fig'])
        close;
        
        clear aha_strain
    end
end
end