load('/scratch/amir/Spatially_Coupled/Database/CIFAR/data_batch_1')
[s,ind] = sort(labels);
class = 1;
CIFAR_DB = [];
temp = zeros(length(data(1,:))/3,3);
for i = 1:length(ind)
    if (s(i) == class)
        temp(:,1) = data(ind(i),1:1024);
        temp(:,2) = data(ind(i),1025:2048);
        temp(:,3) = data(ind(i),2049:3072);
        temp_gray = rgb2gray(temp/255);        
        CIFAR_DB = [CIFAR_DB;temp_gray(:,1)'];
    end
end
save(['/scratch/amir/Spatially_Coupled/Database/CIFAR/db_class_',num2str(class)],'CIFAR_DB');