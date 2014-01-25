function image2gabor(db_file,db_name)


%=============================INITIALIZATION===============================

%-----------------------------Load Dataset---------------------------------
load(db_file);
eval(['dataset_learn=',db_name,';']);
[dataset_size,pattern_length] = size(dataset_learn);
%--------------------------------------------------------------------------

%------------------------Create Directory If Necessary---------------------
if (~exist(['/scratch/amir/Databases/STL_10/Gabor'],'dir'))
    mkdir(['/scratch/amir/Databases/STL_10/Gabor']);
end
%--------------------------------------------------------------------------


%-----------------------Gabor Filter Parameters----------------------------
lambda  = 8;
theta   = 0;
psi     = [pi/2];
gamma   = 0.5;
bw      = 1;
N       = 16;
sigma = 1;
%--------------------------------------------------------------------------

%-------------------------Other Initializations----------------------------
addpath(genpath('/home1/amir/cluster/Common_Library'));                 % Include the library of common functions in the search path
width = sqrt(pattern_length);
image_out = zeros(1,pattern_length);
%--------------------------------------------------------------------------

%==========================================================================

aa = [];
fid = fopen(['/scratch/amir/Databases/STL_10/Gabor/Filtered_image.txt'], 'w');
fclose(fid);
for n = 1:N
    fid = fopen(['/scratch/amir/Databases/STL_10/Gabor/Gabor_out_',num2str(n),'.txt'], 'w');                           
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
        image_out = image_out + a_vector.^2;
        %------------------------------------------------------------------
                        
        %----------------Update Filter Orientations------------------------
        theta = theta + pi/N;
        %------------------------------------------------------------------
        
        %----------------------Store the Results---------------------------
        fid = fopen(['/scratch/amir/Databases/STL_10/Gabor/Gabor_out_',num2str(n),'.txt'], 'a');                           
        fprintf(fid, '%f \t',a_vector);
        fprintf(fid, '\n');
        fclose(fid);                  
        %------------------------------------------------------------------
        
    end
       
    %--------------------------Store the Results---------------------------        
    image_out = image_out.^(.5);
    fid = fopen(['/scratch/amir/Databases/STL_10/Gabor/Filtered_image.txt'], 'a');                                   
    fprintf(fid, '%f \t',image_out);        
    fprintf(fid, '\n');                
    fclose(fid);                          
    aa = [aa;image_out];
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