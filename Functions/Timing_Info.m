function info = Timing_Info(info)

dcm_flds = dir([info.dcm_path,'*img*']);

for j = 1:length(dcm_flds)
    
    dcm_files = dir([dcm_flds(j).folder,'/',dcm_flds(j).name,'/*IM*']);
    
    temp = dicominfo([dcm_files(1).folder,'/',dcm_files(1).name]);
    
    info.percent_rr(j) = temp.NominalPercentageOfCardiacPhase;
    info.gantry_angle(j) = temp.Private_0045_1050;
    info.time_ms(j) = temp.Private_0045_103f;
    
    clear temp
    
end
