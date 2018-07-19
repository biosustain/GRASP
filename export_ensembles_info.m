clear

strain = 'HMP1489';
replicate_list = [0];
time_point_list = [0];

label = '_unconstrained';
input_folder = strcat('./output', label, '/');
output_folder = strcat('./output', label, '/dat_files/');

n_models = 10000;

for time_i = time_point_list
    disp('time_i');
    disp(time_i);
    
    for rep_i = replicate_list
        disp('rep_i');
        disp(rep_i);
    
        model_id = strcat(strain, '_r', int2str(rep_i), '_t', int2str(time_i));
        
        clearvars ensemble
        load(strcat(input_folder, 'ensembleSMC_rejection_', model_id,'.mat'));

        write(cell2table(ensemble.rxns), strcat(output_folder, 'ensembleSMC_rejection_', model_id, '_rxns.dat'));    
        write(cell2table({ensemble.fluxRef; ensemble.fluxRefStd}), strcat(output_folder, 'ensembleSMC_rejection_', model_id, '_fluxRef.dat'));    
              
        fileID = fopen(strcat(output_folder, 'ensembleSMC_rejection_', model_id, '_gibbsRanges.dat'),'w');
        fprintf(fileID,'%6s,%6s\n','dG_min', 'dG_max');
        fprintf(fileID,'%6.5f,%6.5f\n', ensemble.gibbsRanges);
        fclose(fileID);
        parpool(2);
        parfor model_i = 1:n_models
            
            writetable(struct2table(ensemble.populations.models(model_i).rxnParams), strcat(output_folder, 'ensembleSMC_rejection_',model_id, '_m', num2str(model_i), '.dat'))
            %writetable(struct2table(ensemble.populations.models(model_i).rxnParams), strcat(output_folder, 'ensembleSMC_rejection_HMP1489_r0_t0_MA_m', num2str(model_i), '.dat'))
        end     
    end
end
