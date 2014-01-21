% for target_class = 1:no_of_classes
    clc
    target_class = 1;
    no_of_classes = 5;
    alpha0 = 0.95;
    beta0= 0.75;
    theta0 = 0.008;
    N = 576;
    test_flag = 0;
    if (test_flag)
        db_file_in = (['/scratch/amir/Databases/CIFAR_10/Gray_Scale_Class_',num2str(target_class),'/Test_Set/Sparse Filtered/Layer_2/Pooled_3_3/Pooled_Sparse_Filtered_CIFAR_10_Gray_Class_',num2str(target_class),'_Test_Layer_2.mat']);
    else
        db_file_in = (['/scratch/amir/Databases/CIFAR_10/Gray_Scale_Class_',num2str(target_class),'/Train_Set/Sparse Filtered/Layer_2/Pooled_3_3/Pooled_Sparse_Filtered_CIFAR_10_Gray_Class_',num2str(target_class),'_Train_Layer_2.mat']);        
    end
    load(db_file_in);
    CIFAR_DB_orig=final_database_vectorized;
    
    for j = 1:no_of_classes
        %------------------Read the Already Learned Constraints------------
        fid = fopen(['/scratch/amir/Databases/CIFAR_10/Gray_Scale_Class_',num2str(j),'/Train_Set/Sparse Filtered/Layer_2/Pooled_3_3/Learn_Results/Zero_One/N_',num2str(N),'/Weigh_matrix_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'.txt'], 'r');        
        if (fid > -1)                
            W = fscanf(fid, '%f',[N,inf]);                    
            W = W';                    
            fclose(fid);                                                                
            [m,~] = size(W);               
        else    
            error('Invalid input matrix');       
        end
        %------------------------------------------------------------------
        eval(['W_class',num2str(j),'=W;']);
          
%         db_file_in = (['/scratch/amir/Databases/CIFAR_10/Gray_Scale_Class_',num2str(j),'/Train_Set/Sparse Filtered/Layer_2/Pooled_3_3/Pooled_Sparse_Filtered_CIFAR_10_Gray_Class_',num2str(j),'_Train_Layer_2.mat']);            
%         
%         load(db_file_in)
        
        
%         CIFAR_DB=final_database_vectorized;
%         a = sqrt(sum(CIFAR_DB'.*CIFAR_DB'));    
%         dataset_zero_thr = .05;
%         [dataset_size,pattern_length] = size(CIFAR_DB);
% %     CIFAR_DB = CIFAR_DB./(a'*ones(1,pattern_length));
% %     CIFAR_DB = CIFAR_DB.*(abs(CIFAR_DB)>.03);
% %     
%         CIFAR_DB = CIFAR_DB .*(abs(CIFAR_DB)>dataset_zero_thr);
%         CIFAR_DB  = ones(dataset_size,pattern_length).*((CIFAR_DB>dataset_zero_thr)-(CIFAR_DB<-dataset_zero_thr));
%         eval(['CIFAR_DB_class',num2str(j),'=CIFAR_DB;']);
        
    end
    
%     load(['/scratch/amir/Databases/CIFAR_10/Gray_Scale_Class_',num2str(target_class),'/Train_Set/CIFAR_10_train_gray_scale_class_',num2str(target_class),'.mat']);
    
    
%     eval(['CIFAR_DB=CIFAR_10_Gray_DB_class_',num2str(target_class),';']);
    
        
%     load(['/scratch/amir/Databases/CIFAR_10/Gray_Scale_Class_',num2str(target_class),'/Test_Set/CIFAR_10_test_gray_scale_class_',num2str(target_class),'.mat']);
%     
%     eval(['CIFAR_DB=CIFAR_10_Gray_DB_test_class_',num2str(target_class),';']);
 
    a = sqrt(sum(CIFAR_DB_orig'.*CIFAR_DB_orig'));    
    dataset_zero_thr = .05;
    [dataset_size,pattern_length] = size(CIFAR_DB_orig);
%     CIFAR_DB = CIFAR_DB./(a'*ones(1,pattern_length));
%     CIFAR_DB = CIFAR_DB.*(abs(CIFAR_DB)>.03);
%     CIFAR_DB_orig
    CIFAR_DB_orig = CIFAR_DB_orig .*(abs(CIFAR_DB_orig)>dataset_zero_thr);
    CIFAR_DB_orig  = ones(dataset_size,pattern_length).*((CIFAR_DB_orig>dataset_zero_thr)-(CIFAR_DB_orig<-dataset_zero_thr));
    
    
%     a = sqrt(sum(CIFAR_DB'.*CIFAR_DB'));    
%     CIFAR_DB = CIFAR_DB./(a'*ones(1,pattern_length));
    success_count = 0;
%     display(' Mean Projection          Norm(proj)          Sum(proj==0)        mean(Sum(proj<.005))           Length')            
    display('-----------------      --------------      ---------------     ------------------------        ----------')
    for i = 1:dataset_size
        mu = 1+floor(rand*dataset_size);
        pattern = CIFAR_DB_orig(mu,:);
        
        
        max_sum = 0;
        min_sum = 10000;
        for j = 1:no_of_classes
            
                        
            
            eval(['W_tot = W_class',num2str(j),';']);
            
            [m,~] = size(W_tot);
%             if (mean(abs(pattern*W_tot')<.02)>max_sum)
%             if (sum(((pattern*W_tot').*abs(pattern*W_tot')<.02))/length((pattern*W_tot'))<min_sum)
%             proj = zeros(1,m);
            mat_norm = sqrt(sum(W_tot'.*W_tot'));
            proj = abs(W_tot*pattern')./(norm(pattern)*mat_norm');
%             for iki = 1:m
%                 proj(iki) = abs(W_tot(iki,:)*pattern')/norm(W_tot(iki,:))/norm(pattern);
% 
%             end
            
            
%             display(['   ',num2str(mean(proj)),'                ',num2str(norm(proj)),'                 ',num2str(sum(proj==0)),'                   ',num2str(mean(proj<.001)),'                    ',num2str(length(proj))]);
            if (mean(proj)<min_sum)
%             eval(['CIFAR_DB = CIFAR_DB_class',num2str(j),';']);
%             dist = sum(abs(sign((CIFAR_DB-ones(5000,1)*pattern)')));
%             if (mean(dist)<min_sum)
%             if (mean(sum(proj<.001))>max_sum)
%             if (mean(abs(pattern*W_tot'))<min_sum)
%             if (mean(abs(pattern*W_tot').*(abs(pattern*W_tot')>1))<max_sum)
                  pattern_class = j;
%                 max_sum = mean(abs(pattern*W_tot')<.02);
%                   max_sum = mean(abs(pattern*W_tot').*(abs(pattern*W_tot')>1));
%                 min_sum = sum(((pattern*W_tot').*abs(pattern*W_tot')<.02))/length((pattern*W_tot'));
%                 min_sum = mean(abs(pattern*W_tot'));
                    min_sum = mean(proj);
%                     min_sum = mean(dist);
%                     max_sum = mean(sum(proj<.001));
            end            
        end
%         display(' ')
        if (pattern_class==target_class)
            success_count = success_count+1;
        end
        
        %-------------------------Display Progress-------------------------
        if (mod(i,50)==0)
            display(' ')
            display(['Success rate so far = ',num2str( success_count/i)]);
            
        end
        %------------------------------------------------------------------
    end
    success_rate = success_count/i
    