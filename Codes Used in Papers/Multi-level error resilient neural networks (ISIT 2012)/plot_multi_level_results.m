%==========================================================================
%********************************READ ME***********************************
%==========================================================================


%--------------------------------INPUTS------------------------------------
% N: The number of pattern nodes in each LOCAL network. It can be a matrix.
% K: The dimension of the subspace where LOCAL sub-patterns come from. It can be a matrix.
% -------------------------------------------------------------------------


%--------------------------------OUTPUTS-----------------------------------
% none
% -------------------------------------------------------------------------


%--------------------------FUNCTION DESCRIPTION----------------------------
% This piece of code reads the results of the recall process of the
% multi-level neural constraint enforcing network from the specified file 
% and plots the bit and block error rates.
% It plots Pattern Error Rate (PER) and Bit Error Rate (BER) as a function
% of the initial number of erroneous symbols as well as the gain we get by
% adding a second correction level (in version 3 and higher).
%--------------------------------------------------------------------------

%==========================================================================
%==========================================================================

function plot_multi_level_results(N,K,version)

%============================INITIALIZATION================================
figure_legend = [];                         % Stores the legend for the PER and BER figures
figure_legend_gain = [];                    % Stores the legend for the gain figure
%==========================================================================


%==============PROCESS THE INPUT FILES AND DISPLAY THE RESULTS=============
for l = 1:length(N)

    %---------Read Error Performance Results From the Results File---------
    fid = fopen(['neural_multi_level_results_v',num2str(version),'_N_',num2str(N(l)),'_K_',num2str(K(l)),'.txt'], 'r');             % This is for single-level versions
    if (version <= 2)                   
        results = fscanf(fid, '%s \t %d \t %s \t %f \t %s \t %f',[43,inf]);
    else                                                    
        results = fscanf(fid, '%s \t %d \t %s \t %f \t %s \t %f \t %s \t %f \t \t %s \t %f \t %s \t %f \t %s \t %d',[155,inf]);     % This is for multi-level versions
    end
    fclose(fid);
    %----------------------------------------------------------------------

    %--------------Read Number of Learn Phases From the File---------------
    % fid = fopen(['learn_no_N_',num2str(N(l)),'_K_',num2str(K(l)),'.txt'], 'r');
    % learn_results = fscanf(fid, '%s \t %d \t %s \t %d',[22,inf]);
    % fclose(fid);
    %----------------------------------------------------------------------

    %-----------------------Pre-Process the Results------------------------        
    if (version <=2)                                                % For single-level networks
        unprocessed_error_bits = results(9,:);
        first_level_unprocessed_pattern_error_rate = results(28,:);
        first_level_unprocessed_bit_error_rate = results(43,:);
    elseif (version ==3)                                            % For single-level networks
        unprocessed_success_count = results(155,:);
        unprocessed_error_bits = results(9,:);
        first_level_unprocessed_pattern_error_rate = results(36,:);
        first_level_unprocessed_bit_error_rate = results(59,:);
        second_level_unprocessed_pattern_error_rate = results(113,:);
        second_level_unprocessed_bit_error_rate = results(141,:);
    end
    %     unprocessed_learn_no_err = learn_results(9,:);
    %     unprocessed_learn_no = learn_results(22,:);   
    %----------------------------------------------------------------------

    %-------------------------Other Initialization-------------------------
    error_bits = [0];
    error_bits_count = [1];
    first_level_pattern_error_rate = [0];
    first_level_bit_error_rate = [0];
    second_level_pattern_error_rate = [0];
    second_level_bit_error_rate = [0];
    %learn_no_err = [0];
    %learn_no = [0];
    %----------------------------------------------------------------------

    %======================================================================

    

    %=====================PROCESS THE RESULTS==============================
    for i = 1:length(unprocessed_error_bits)
    
        %-Check If This Number of Initial Errors Has Been Considered Before-
        processed_flag = 0;
        for j = 1:length(error_bits)
            if (error_bits(j) == unprocessed_error_bits(i))                     % If this number of initial errors has been encountered before...
                processed_flag = 1;                                             %...then trigeer the appropriate flag.
                break;
            end
        end
        %------------------------------------------------------------------
    
        %-------------------If Not Processed Earlier-----------------------
        if (processed_flag == 0)
            error_bits = [error_bits,unprocessed_error_bits(i)];                            % Add the new number of initial errors to the appropriate list.
            error_bits_count = [error_bits_count,1];                                        % Set the number of times this error number has been visited to 1.
            first_level_pattern_error_rate = [first_level_pattern_error_rate,first_level_unprocessed_pattern_error_rate(i)];    % Add the pattern error rate to the appropriate list.
            first_level_bit_error_rate = [first_level_bit_error_rate,first_level_unprocessed_bit_error_rate(i)];                % Add the bit error rate to the appropriate list.
            if (version >=3)
                second_level_pattern_error_rate = [second_level_pattern_error_rate,second_level_unprocessed_pattern_error_rate(i)];    % Add the pattern error rate for the second level to the appropriate list.
                second_level_bit_error_rate = [second_level_bit_error_rate,second_level_unprocessed_bit_error_rate(i)];                % Add the bit error rate for the second level to the appropriate list.
            end
            % learn_no = [learn_no,unprocessed_learn_no(i)];
            % learn_no_err = [learn_no_err,unprocessed_learn_no_err(i)];
        %------------------------------------------------------------------
    
        %---------------------If Encountered Earlier-----------------------
        else
            error_bits_count(j) = error_bits_count(j)+1;                                    % Increase the number of times this error number has been visited to by 1.
            first_level_pattern_error_rate(j)=first_level_pattern_error_rate(j) + first_level_unprocessed_pattern_error_rate(i);% Add the pattern error rate to the right place at the appropriate list.
            first_level_bit_error_rate(j) = first_level_bit_error_rate(j)+first_level_unprocessed_bit_error_rate(i);            % Add the bit error rate to the right place at the appropriate list.
            if (version >=3)
                second_level_pattern_error_rate(j)=second_level_pattern_error_rate(j) + second_level_unprocessed_pattern_error_rate(i);% Add the pattern error rate for the second level to the right place at the appropriate list.
                second_level_bit_error_rate(j) = second_level_bit_error_rate(j)+second_level_unprocessed_bit_error_rate(i);            % Add the bit error rate for the second level to the right place at the appropriate list.
            end
            % learn_no(j) = learn_no(j)+unprocessed_learn_no(i);
            %  learn_no_err(j) = learn_no_err(j)+unprocessed_learn_no_err(i);
        end
        %------------------------------------------------------------------
    end

    %----------------------Normalize the Error Rates-----------------------
    first_level_pattern_error_rate = first_level_pattern_error_rate./error_bits_count;                  % Normalize the pattern error rate for the first level
    first_level_bit_error_rate = first_level_bit_error_rate./error_bits_count;                          % Normalize the bit error rate for the first level
    second_level_pattern_error_rate = second_level_pattern_error_rate./error_bits_count;                % Normalize the pattern error rate for the second level
    second_level_bit_error_rate = second_level_bit_error_rate./error_bits_count;                        % Normalize the bit error rate for the second level
    
    error_bits = sort(error_bits);                                                                      % Sort the initial number of erroneous nodes for plotting purposes                                                                  
    first_level_bit_error_rate = sort(first_level_bit_error_rate);                                      % Sort the BER values for the first level for plotting purposes
    second_level_bit_error_rate = sort(second_level_bit_error_rate);                                    % Sort the BER values for the second level for plotting purposes
    
    first_level_pattern_error_rate = sort(first_level_pattern_error_rate);                              % Sort the PER values for the first level for plotting purposes                              
    second_level_pattern_error_rate = sort(second_level_pattern_error_rate);                            % Sort the PER values for the second level for plotting purposes      
    %----------------------------------------------------------------------

    %======================================================================

    %========================DISPLAY THE RESULTS===========================

    %---------------------Display Bit Error Rate---------------------------
    
    %-----------Select the approprate figure-------------
    if (l ==1)
        h_bit = figure;
    else
        figure(h_bit);
    end
    %----------------------------------------------------
    if (version <=2)                                        % For single-level networks
        plot(error_bits,first_level_bit_error_rate,'b--s','LineWidth',3,'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','y','Color',[tanh(l/length(N)),l/length(N),1-(l/length(N))]);
    else                                                    % For double-level networks
        plot(error_bits,first_level_bit_error_rate,'b--s','LineWidth',3,'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','y','Color',[tanh(l/length(N)),l/length(N),1-(l/length(N))]);   % Plot the BER of the 1st level
        hold on
        plot(error_bits,second_level_bit_error_rate,'b-s','LineWidth',3,'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','y','Color',[tanh(l/length(N)),l/length(N),1-(l/length(N))]);   % Plot the BER of the 2nd level
    end
    title('Bit error rate vs. initial number of erroneous nodes','fontsize',16);
    xlabel('Initial number of erroneous nodes','fontsize',16);
    ylabel('Final bit error rate','fontsize',16);
    set(gca,'FontSize',16);
    hold on
    %----------------------------------------------------------------------

    %--------------------Display Bit Error Rate----------------------------
    
    %-----------Select the approprate figure-------------
    if (l ==1)
        h_pat = figure;
    else
        figure(h_pat);
    end
    %----------------------------------------------------
    
    if (version <=2)
        plot(error_bits,first_level_pattern_error_rate,'r--s','LineWidth',3,'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','y','Color',[tanh(l/length(N)),l/length(N),1-(l/length(N))]);
    else
        plot(error_bits,first_level_pattern_error_rate,'r--s','LineWidth',3,'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','y','Color',[tanh(l/length(N)),l/length(N),1-(l/length(N))]);   % Plot the PER of the 1st level   
        hold on
        plot(error_bits,second_level_pattern_error_rate,'r-s','LineWidth',3,'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','y','Color',[tanh(l/length(N)),l/length(N),1-(l/length(N))]);   % Plot the PER of the 2nd level
    end
    title('Pattern error rate vs. initial number of erroneous nodes','fontsize',16);
    xlabel('Initial number of erroneous nodes','fontsize',16);
    ylabel('Final pattern error rate','fontsize',16);
    set(gca,'FontSize',16);
    hold on
    %----------------------------------------------------------------------


    %----------------Display the Gain in Pattern Error Rate----------------
    if (version > 2)                                        % For multi-level networks
        %-----------Select the approprate figure-------------
        if (l ==1)
            h_gain=figure;
        else
            figure(h_gain);
        end
        %----------------------------------------------------
        
        
        %---------------Calculate the Gain-------------------
        BER_Gain = first_level_bit_error_rate./second_level_bit_error_rate;
        PER_Gain = first_level_pattern_error_rate./second_level_pattern_error_rate;
        for i = 1:length(PER_Gain)
            if (isnan(PER_Gain(i)))
                PER_Gain(i) = 0;
            end
             if (isnan(BER_Gain(i)))
                BER_Gain(i) = 0;
            end
        end
        %----------------------------------------------------
        
        plot(error_bits,(PER_Gain),'g--s','LineWidth',3,'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','y','Color',[tanh(l/length(N)),l/length(N),1-(l/length(N))]);                           % Plot the gain for PER
        hold on 
        plot(error_bits,(BER_Gain),'g-s','LineWidth',3,'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','y','Color',[tanh(l/length(N)),l/length(N),1-(l/length(N))]);                            % Plot the gain for BER
        title('The gain in bit and pattern error rates due to addition of a second layer','fontsize',16);
        xlabel('Initial number of erroneous nodes','fontsize',16);
        ylabel('Gain','fontsize',16);
        set(gca,'FontSize',16);
        hold on    
        %----------------------------------------------------------------------
    end
    
    %-------------------------Constructu the Legend------------------------
    if (version <=2)
        temp_leg = ['N=',num2str(N(l))];
        figure_legend = [figure_legend;temp_leg];
    else
        temp_leg = ['N=',num2str(N(l)),' - Level 1'];
        figure_legend = [figure_legend;temp_leg];
        temp_leg = ['N=',num2str(N(l)),' - Level 2'];
        figure_legend = [figure_legend;temp_leg];
        
        temp_leg_gain = ['N=',num2str(N(l)),' - PER Gain'];
        figure_legend_gain = [figure_legend_gain;temp_leg_gain];
        temp_leg_gain = ['N=',num2str(N(l)),' - BER Gain'];
        figure_legend_gain = [figure_legend_gain;temp_leg_gain];
    end
    
    %----------------------------------------------------------------------
    
    %======================================================================
end

%======================ADD THE LEGEND TO THE FIGURES=======================
figure(h_bit)
legend(figure_legend)
figure(h_pat)
legend(figure_legend)
figure(h_gain)
legend(figure_legend_gain)
%==========================================================================
