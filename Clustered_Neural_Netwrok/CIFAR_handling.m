load('/scratch/amir/Databases/CIFAR_10/test_batch.mat')
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

a = sqrt(sum(CIFAR_DB'.*CIFAR_DB'));
CIFAR_DB = CIFAR_DB./(a'*ones(1,1024));
CIFAR_DB = CIFAR_DB.*(abs(CIFAR_DB)>.03);
% eval(['CIFAR_10_Gray_DB_test_class_',num2str(class),'=CIFAR_DB;']);
% save(['/scratch/amir/Databases/CIFAR_10/Gray_Scale_Class_',num2str(class),'/Test_Set/CIFAR_10_test_gray_scale_class_',num2str(class)],['CIFAR_10_Gray_DB_test_class_',num2str(class)]);