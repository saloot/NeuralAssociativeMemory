%-------------------------Add Necessary Libraries--------------------------
addpath(genpath('/home1/amir/cluster/Common_Library'));                          % Include the library of common functions in the search path
%--------------------------------------------------------------------------


load('/home/amir/Hesam/Research/Biocoding Project/Databases/STL-10/test.mat');
Z = reshape(X,8000,96,96,3);
temp = zeros(96,96,3);
STL_10_Gray_DB = [];
class = 1;
for i = 1:8000
    if (y(i) == class)
        temp(:,:,:) = Z(i,:,:,:);                            % Read each colored image from the dataset
        temp_gray = rgb2gray(temp/255);                     % Transform it to gray scale image
        STL_10_Gray_DB = [STL_10_Gray_DB;matrix2vector(temp_gray)];   % Save the result as rows of a new dataset
    end
end
save(['/home/amir/Hesam/Research/Biocoding Project/Databases/STL-10/STL_10_test_gray_scale_class_',num2str(class),'.mat'],'STL_10_Gray_DB');