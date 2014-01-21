load('/scratch/amir/Databases/CIFAR_10/Gray_Scale_Class_1/Train_Set/Sparse Filtered/Layer_2/Pooled_3_3/Pooled_Sparse_Filtered_CIFAR_10_Gray_Class_1_Train_Layer_2.mat')
[COEFF,SCORE,latent] = princomp(final_database_vectorized);
 SCORE2=SCORE(:,1:3);
 COEFF2=COEFF(:,1:3);
 temp = final_database_vectorized-ones(5000,1)*mean(final_database_vectorized);