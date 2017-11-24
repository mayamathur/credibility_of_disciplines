data=csvread('../../Applications/MinimumWagev2/BWWP-MinWageData.csv',1,0,[1,0,1008,6]);
X=data(:,6);
sigma=data(:,7);
ID=data(:,1);
published=data(:,3);
year=data(:,2);
n=length(X);

output=[X sigma ID published year];

%Read in list of names
names=readtable('../../Applications/MinimumWagev2/BWWP-MinWageText.csv');
%Sort list of names
[output, idx]=sortrows(output,1);
sorted_names=names(idx,:);

sorted_names(output(:,2)==0,:)=[];
output(output(:,2)==0,:)=[];

%Output data and names
csvwrite('../../Applications/MinimumWagev2/cleaned_minwage_data.csv',output)
writetable(sorted_names(:,:),'../../Applications/MinimumWagev2/sorted_text.csv')

