function Mesh = Bullseye_Plots(Mesh,info)

load('RSct_Colormap.mat')

rows = ceil(length(info.timeframes)/info.polar_NoOfCols);
c = 1;
figure('pos',[10 10 2400 1800])

for j = info.timeframes
    
    data = Mesh(j).Polar_Data;
    
    % Interpolating raw data to achieve desired bullseye plot resolution in the azimuthal direction
    for k = 1:numel(data)
                
        if k < info.lvot_limit(j)
            temp = interp1(data{k}(2,:),data{k}(1,:),linspace(min(data{k}(2,:)),max(data{k}(2,:)),info.polar_res(1)));
        else
           ind = ~isnan(data{k}(1,:));
           temp = interp1(data{k}(2,ind),data{k}(1,ind),linspace(min(data{k}(2,ind)),max(data{k}(2,ind)),info.polar_res(1)));
        end   
        
        bull_data(k,:) = temp;
        clear temp
    end
    
    %Interpolating bull_data to achieve desired resolution in slice direction
    temp = interp1(1:numel(data),bull_data,linspace(1,numel(data),info.polar_res(2)));
    temp = flip(temp);
       
    % Putting last value at first and first value at last to ensure a smooth bullseye plot
    temp = [temp(:,end), temp];
    temp = [temp, temp(:,2)];
    temp = temp';
    
    Mesh(j).Bullseye = temp;
    
    subplot(rows,info.polar_NoOfCols,c)
    bullseye(temp); colormap(cmap); caxis(info.RSct_limits);
    h = colorbar; set(h,'position',[0.93 0.05 0.015 0.9]); h.FontSize = 25; h.FontWeight = 'bold'; h.Ticks = -1:0.05:1;
    set(get(h,'Title'),'string','RS_{CT}');
    title([num2str(info.percent_rr(c)),' %'],'FontSize',30','FontWeight','bold')
    
    c = c + 1;
    clear bull_data temp
end
savefig([info.save_path,info.patient,'_Bullseye_RSct.fig'])
close;

if info.endo_strains
        
    c = 1;
    h1 = figure('pos',[10 10 2400 1800]);
    h2 = figure('pos',[10 10 2400 1800]);

    for j = info.timeframes

        data = Mesh(j).Polar_Data;

        % Interpolating raw data to achieve desired bullseye plot resolution in the azimuthal direction
        for k = 1:numel(data)

            if k < info.lvot_limit(j)
                temp_ecc = interp1(data{k}(2,:),data{k}(4,:),linspace(min(data{k}(2,:)),max(data{k}(2,:)),info.polar_res(1)));
                temp_ell = interp1(data{k}(2,:),data{k}(5,:),linspace(min(data{k}(2,:)),max(data{k}(2,:)),info.polar_res(1)));
            else
               ind = ~isnan(data{k}(4,:));
               temp_ecc = interp1(data{k}(2,ind),data{k}(4,ind),linspace(min(data{k}(2,ind)),max(data{k}(2,ind)),info.polar_res(1)));
               temp_ell = interp1(data{k}(2,ind),data{k}(5,ind),linspace(min(data{k}(2,ind)),max(data{k}(2,ind)),info.polar_res(1)));
            end   

            bull_data_ecc(k,:) = temp_ecc;
            bull_data_ell(k,:) = temp_ell;
            clear temp_ecc temp_ell
        end

        %Interpolating bull_data to achieve desired resolution in slice direction
        temp_ecc = interp1(1:numel(data),bull_data_ecc,linspace(1,numel(data),info.polar_res(2)));
        temp_ecc = flip(temp_ecc);
        
        temp_ell = interp1(1:numel(data),bull_data_ell,linspace(1,numel(data),info.polar_res(2)));
        temp_ell = flip(temp_ell);

        % Putting last value at first and first value at last to ensure a smooth bullseye plot
        temp_ecc = [temp_ecc(:,end), temp_ecc];
        temp_ecc = [temp_ecc, temp_ecc(:,2)];
        temp_ecc = temp_ecc';
        
        temp_ell = [temp_ell(:,end), temp_ell];
        temp_ell = [temp_ell, temp_ell(:,2)];
        temp_ell = temp_ell';

        Mesh(j).Bullseye_Ecc = temp_ecc;
        Mesh(j).Bullseye_Ell = temp_ell;

        set(0,'CurrentFigure',h1)
        subplot(rows,info.polar_NoOfCols,c)
        bullseye(temp_ecc); colormap(cmap); caxis(info.RSct_limits);
        h = colorbar; set(h,'position',[0.93 0.05 0.015 0.9]); h.FontSize = 25; h.FontWeight = 'bold'; h.Ticks = -1:0.05:1;
        set(get(h,'Title'),'string','E_{cc}');
        title([num2str(info.percent_rr(c)),' %'],'FontSize',30','FontWeight','bold')
        
        set(0,'CurrentFigure',h2)
        subplot(rows,info.polar_NoOfCols,c)
        bullseye(temp_ell); colormap(cmap); caxis(info.RSct_limits);
        h = colorbar; set(h,'position',[0.93 0.05 0.015 0.9]); h.FontSize = 25; h.FontWeight = 'bold'; h.Ticks = -1:0.05:1;
        set(get(h,'Title'),'string','E_{ll}');
        title([num2str(info.percent_rr(c)),' %'],'FontSize',30','FontWeight','bold')

        c = c + 1;
        clear bull_data_ecc bull_data_ell temp_ecc temp_ell
    end

    savefig(h1,[info.save_path,info.patient,'_Bullseye_Ecc.fig'])
    savefig(h2,[info.save_path,info.patient,'_Bullseye_Ell.fig'])
    
    close all;
end

end