data=csvread('../../Applications/PsychExperiments/SelectedPsychColumns.csv',1,0,[1,0,167,14]);

%Read in list of names
names=readtable('../../Applications/PsychExperiments/PsychLabels.csv');

%Drop rows for which do not have replication study
names(data(:,5)==0,:)=[];
data(data(:,5)==0,:)=[];
%Drop two observations which report zero p-value in replication
%data(data(:,6)==0,:)=[];

Obs1=data(:,1);
Tails1=data(:,2);
Pval1=data(:,3);
Pval1_alt=data(:,4);

Obs2=data(:,5);
Tails2=data(:,6);
Pval2=data(:,7);
Pval2_alt=data(:,8);

signflip=data(:,9);

Estimate1=data(:,9);
Estimate2=data(:,10);

Original_researcher_approval=data(:,11);

dof1O=data(:,12);
dof1R=data(:,13);
dof2O=data(:,14);
dof2R=data(:,15);

%Impute z-statistics from p-values.  When imputation produces infinite Z
%statistic, use alternative p_value
Z1=norminv(1-Pval1/2);
Z1(Tails1==1)=norminv(1-Pval1(Tails1==1));
Z1(Z1==inf)=norminv(1-Pval1_alt(Z1==inf)/2);
Z1(Tails1==1&Z1==inf)=norminv(1-Pval1_alt(Tails1==1&Z1==inf));

Z2=norminv(1-Pval2/2);
Z2(Tails2==1)=norminv(1-Pval2(Tails2==1));
Z2(Z2==inf)=norminv(1-Pval1_alt(Z2==inf)/2);
Z2(Tails2==1&Z2==inf)=norminv(1-Pval2_alt(Tails2==1&Z2==inf));
Z2=Z2.*(1-2*(Estimate2<0));

%Use Fisher Transformation to obtain approximately normal standardized
%estimates
Fisher=@(x) 0.5*log((1+x)./(1-x));
Estimate1_final=Fisher(Estimate1);
Estimate2_final=Fisher(Estimate2);

SE1=abs(Estimate1_final)./abs(Z1);
SE2=abs(Estimate2_final)./abs(Z2);

output=[Estimate1_final SE1 Estimate2_final SE2 Z1 Original_researcher_approval==1 dof1O dof1R dof2O dof2R];
%Sort list of names
[output, idx]=sortrows(output,5);
sorted_names=names(idx,:);
%drop observations for which do not have point estimate
sorted_names(output(:,3)==0,:)=[];
output(output(:,3)==0,:)=[];

%Drop column Z statistics (which included only for sorting purposes)
output(:,5)=[];

%Output data and names
csvwrite('../../Applications/PsychExperiments/cleaned_psych_data.csv',output)
writetable(sorted_names(:,1),'../../Applications/PsychExperiments/sorted_names.csv')
writetable(sorted_names(:,3),'../../Applications/PsychExperiments/sorted_1stAuthor.csv')