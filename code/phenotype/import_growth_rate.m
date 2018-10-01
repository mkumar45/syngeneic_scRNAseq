function growth_rates = import_growth_rate( data_path )

    growth_tbl = readtable(data_path, 'Sheet','combined',....
                           'ReadRowNames',false,'ReadVariableNames',true);
                       
    unique_models = unique( growth_tbl.Model );
    num_models = length(unique_models);
    rates = cell(num_models,1);
%     fig1 = figure();
%     fig2 = figure();
    
    for n = 1:num_models

        model_name = unique_models{n};
        model_idx = strcmp( growth_tbl.Model, model_name );
        
        model_tbl = growth_tbl( model_idx, : );
            
        treat_idx = strcmp( model_tbl.expGroup,  'Saline' );
        treat_tbl = model_tbl( treat_idx, : );
        % ddd
        unique_animal_id = unique( treat_tbl.animalId );
        num_animals = length(unique_animal_id);

        rates{n} = zeros(num_animals,1); 
        for a = 1:num_animals

            animal_id = unique_animal_id(a);
            animal_idx = treat_tbl.animalId == animal_id;

            time = treat_tbl{ animal_idx, 'fromInjection' };
            volume = treat_tbl{ animal_idx, 'volume' };
            log_volume = log10(1+volume); 


            linear_coeff = polyfit( time, log_volume, 1);
            rates{n}(a) = linear_coeff(1);
%             figure(fig2)
%             subplot(2,num_models,n);
%             plot(time,volume); hold on;
%             box on; axis square;
%             xlabel('Time (days)'); ylabel('Volumne (mm^3)');
%             title(model_name)

%             subplot(2,num_models, num_models+n)
%             plot(time,log_volume,'k','LineWidth',1.5); hold on;
%             box on; axis square;
%             xlabel('Time (days)'); ylabel('log_{10} Volumne (mm^3)');
%             ylim([0,4])
%             t = 1:max(time);
%             v = linear_coeff(1)*t + linear_coeff(2);
            %plot(t,v,'k.')

        end 
        

%         figure(fig1);
%         scatter( n.*ones(num_animals,1), rates{n} ,500,'k', 'filled' )
%         hold on    
%         plot( [n-0.2,n+.2], [median( rates{n}), median( rates{n})], 'k', 'LineWidth', 2 )
%         xlim([0 num_models+1])
%         set(gca,'XTick',1:num_models)
%         set(gca,'XTickLabels',unique_models)
%         box on; axis square;
%         ylabel('Growth rates')
%         set(gca,'FontSize',16)

    end

    growth_rates = table( cellfun(@mean,rates), 'VariableNames', {'growth_rate'},...
                      'RowNames', unique_models);
    
end