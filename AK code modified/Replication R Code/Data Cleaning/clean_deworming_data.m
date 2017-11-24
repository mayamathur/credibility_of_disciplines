data=csvread('../../Applications/deworming/DewormingTable1.csv',0,0,[0,0,21,1]);
X=data(:,1);
sigma=data(:,2);
n=length(X);

output=[X sigma];

%Read in list of names
names=readtable('../../Applications/deworming/DewormingLabels.csv');
%Sort list of names
[output, idx]=sortrows(output,1);
sorted_names=names(idx,:);

%create cluster identifiers
count=1;
cluster_ID=zeros(n,1);
for m=1:n
    if cluster_ID(m)==0
        for m2=1:n
            if isequal(sorted_names(m2,1),sorted_names(m,1))==1
            cluster_ID(m2,:)=count;
            end
        end
       count=count+1;
    end
end
output=[output cluster_ID];

%Output data and names
csvwrite('../../Applications/deworming/cleaned_deworming_data.csv',output)
writetable(sorted_names(:,:),'../../Applications/deworming/sorted_names.csv')

