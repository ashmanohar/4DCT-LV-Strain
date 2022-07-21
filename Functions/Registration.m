function Mesh = Registration(Mesh,info,opts)

% Defining a time vector without the chosen template phase for registration
time_vector = info.timeframes; time_vector(info.timeframes == info.template) = [];

%Initializing template CPD to template mesh
Mesh(info.template).CPD = Mesh(info.template).crop_verts;
Mesh(info.template).CPD_noSmooth = Mesh(info.template).crop_verts;

for j = time_vector %Time frame loop
    
    [T,C] = cpd_register(Mesh(j).vertices(Mesh(j).indxs,:),Mesh(info.template).vertices(Mesh(info.template).indxs,:),opts);
    
    % Initializing the registered meshes
    Mesh(j).CPD = zeros(size(Mesh(info.template).vertices));
    % NaN-ing vertices belonging to planes and initializing the rest to the output of CPD
    Mesh(j).CPD(~Mesh(info.template).indxs,:) = NaN; Mesh(j).CPD(Mesh(info.template).indxs,:) = T.Y;
    
    % Saving a copy of the non-smoothed CPD vertices
    Mesh(j).CPD_noSmooth = Mesh(j).CPD;
    
    dummy = find(Mesh(j).indxs);
    
    %Initializiing the correspondence matrix
    Mesh(j).Correspondence = zeros(size(Mesh(info.template).vertices,1),1);
    % NaN-ing vertices belonging to planes and initializing the rest to the output correspondence of CPD
    Mesh(j).Correspondence(~Mesh(info.template).indxs) = NaN; Mesh(j).Correspondence(Mesh(info.template).indxs) = dummy(C);
    
    clear dummy
    
    disp(['Done registering time frame ',num2str(j)])
    
end


if info.smooth_verts
    
    th_sort = linspace(0,2*pi,length(info.timeframes)+1);
    nHARMO = 3; %no. of modes

    % Tpeak period of my signal
    tI = linspace(th_sort(1),th_sort(1)+2*pi,length(info.timeframes)+1)';
    
    verts_vect = 1:length(Mesh(info.template).CPD);
    verts_vect = verts_vect(Mesh(info.template).indxs);

    for j1 = verts_vect

        for j = 1:length(info.timeframes)

            X(1,j) = Mesh(info.timeframes(j)).CPD(j1,1);
            X(2,j) = Mesh(info.timeframes(j)).CPD(j1,2);
            X(3,j) = Mesh(info.timeframes(j)).CPD(j1,3);

        end
        
        % Enforcing periodicity
        X = [X X(:,1)];
        
        % Average value
        Avg_X(:,1) = trapz(th_sort,X(1,:))./(2*pi).*ones(length(info.timeframes)+1,1);
        Avg_X(:,2) = trapz(th_sort,X(2,:))./(2*pi).*ones(length(info.timeframes)+1,1);
        Avg_X(:,3) = trapz(th_sort,X(3,:))./(2*pi).*ones(length(info.timeframes)+1,1);

        for i = 1:nHARMO
            Si = sin(i*th_sort);
            Ci = cos(i*th_sort);

            Avg_X(:,1) = Avg_X(:,1) + trapz(th_sort,X(1,:).*Si)./trapz(th_sort,Si.^2).*...
                sin(i*tI);
            Avg_X(:,1) = Avg_X(:,1) + trapz(th_sort,X(1,:).*Ci)./trapz(th_sort,Ci.^2).*...
               cos(i*tI);

            Avg_X(:,2) = Avg_X(:,2) + trapz(th_sort,X(2,:).*Si)./trapz(th_sort,Si.^2).*...
                sin(i*tI);
            Avg_X(:,2) = Avg_X(:,2) + trapz(th_sort,X(2,:).*Ci)./trapz(th_sort,Ci.^2).*...
               cos(i*tI);

            Avg_X(:,3) = Avg_X(:,3) + trapz(th_sort,X(3,:).*Si)./trapz(th_sort,Si.^2).*...
                sin(i*tI);
            Avg_X(:,3) = Avg_X(:,3) + trapz(th_sort,X(3,:).*Ci)./trapz(th_sort,Ci.^2).*...
               cos(i*tI);
        end

        for j = 1:length(info.timeframes)

            Mesh(info.timeframes(j)).CPD(j1,1) = Avg_X(j,1);
            Mesh(info.timeframes(j)).CPD(j1,2) = Avg_X(j,2);
            Mesh(info.timeframes(j)).CPD(j1,3) = Avg_X(j,3);

        end

        clear Avg_X X

    end
end    