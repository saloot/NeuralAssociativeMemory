function gabor_feat_ext_exec(L_horiz,L_vert,cluster_index_horiz,cluster_index_vert)


%=============================INITIALIZATION===============================

%-----------------------------Load Dataset---------------------------------
db_file = ['/home/amir/Hesam/Research/Biocoding Project/Databases/STL-10/Projected/Class_1/L_horiz_',num2str(L_horiz),'_L_vert_',num2str(L_vert),'/dataset_chopped_',num2str(cluster_index_horiz),'_',num2str(cluster_index_vert),'.mat'];
db_name = 'dataset_chopped';
load(db_file);
eval(['dataset_learn=',db_name,';']);
[dataset_size,pattern_length] = size(dataset_learn);
%--------------------------------------------------------------------------

%------------------------Create Directory If Necessary---------------------
if (~exist(['/scratch/amir/Databases/STL_10/L_horiz_',num2str(L_horiz),'_L_vert_',num2str(L_vert),'/Patch_',num2str(cluster_index_horiz),'_',num2str(cluster_index_vert)],'dir'))
    mkdir(['/scratch/amir/Databases/STL_10/L_horiz_',num2str(L_horiz),'_L_vert_',num2str(L_vert),'/Patch_',num2str(cluster_index_horiz),'_',num2str(cluster_index_vert)]);
end
%--------------------------------------------------------------------------


%-----------------------Gabor Filter Parameters----------------------------
lambda  = 8;
theta   = 0;
psi     = [pi/2];
gamma   = 0.5;
bw      = 1;
N       = 32;
sigma = 1;
%--------------------------------------------------------------------------

%-------------------------Other Initializations----------------------------
addpath(genpath('/home1/amir/cluster/Common_Library'));                 % Include the library of common functions in the search path
width = sqrt(pattern_length);
%--------------------------------------------------------------------------

%==========================================================================

aa = [];
fid = fopen(['/scratch/amir/Databases/STL_10/L_horiz_',num2str(L_horiz),'_L_vert_',num2str(L_vert),'/Patch_',num2str(cluster_index_horiz),'_',num2str(cluster_index_vert),'/Features_vector.txt'], 'w');
fclose(fid);
for n = 1:N
    fid = fopen(['/scratch/amir/Databases/STL_10/L_horiz_',num2str(L_horiz),'_L_vert_',num2str(L_vert),'/Patch_',num2str(cluster_index_horiz),'_',num2str(cluster_index_vert),'/Gabor_out_',num2str(n),'.txt'], 'w');                           
   fclose(fid); 
end
%================================MAIN LOOP=================================
for i = 1:dataset_size
    %----------------------Read Image From File----------------------------
    img_in = vector2matrix(dataset_learn(i,:),width);
    %----------------------------------------------------------------------
        
    feature_vector = zeros(1,N);
    for n = 1:N
        
        %-----------------------Create Filter------------------------------
        gb = gabor_fn2(sigma,theta,lambda,psi,gamma);
        %------------------------------------------------------------------
        
        %----------------------Calculate Filter Output---------------------
        a = imfilter(img_in, gb, 'symmetric');        
        a_vector = matrix2vector(a);
        feature_vector(n) = norm(abs(tanh(a_vector)));
        %------------------------------------------------------------------
                        
        %----------------Update Filter Orientations------------------------
        theta = theta + pi/N;
        %------------------------------------------------------------------
        
        %----------------------Store the Results---------------------------
        fid = fopen(['/scratch/amir/Databases/STL_10/L_horiz_',num2str(L_horiz),'_L_vert_',num2str(L_vert),'/Patch_',num2str(cluster_index_horiz),'_',num2str(cluster_index_vert),'/Gabor_out_',num2str(n),'.txt'], 'a');                           
        fprintf(fid, '%f \t',a_vector);
        fprintf(fid, '\n');
        fclose(fid);                  
        %------------------------------------------------------------------
        
    end
       
    %--------------------------Store the Results---------------------------        
    fid = fopen(['/scratch/amir/Databases/STL_10/L_horiz_',num2str(L_horiz),'_L_vert_',num2str(L_vert),'/Patch_',num2str(cluster_index_horiz),'_',num2str(cluster_index_vert),'/Features_vector.txt'], 'a');                                   
    fprintf(fid, '%f \t',feature_vector);        
    fprintf(fid, '\n');                
    fclose(fid);                          
    aa = [aa;feature_vector];
    %----------------------------------------------------------------------
    
    %------------------------Display Progress------------------------------
    if (mod(i,50) == 0)
        display(i);
        size(aa)
        rank(aa)
    end
    %----------------------------------------------------------------------
    
end   
%==========================================================================
111;