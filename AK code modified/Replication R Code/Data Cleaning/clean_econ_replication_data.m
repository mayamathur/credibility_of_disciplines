data=csvread('../../Applications/EconExperiments/ReplicationEconComplete.csv',0,0,[0,0,17,14]);
n=18;

%create variables
Estimate1=data(:,2);
Estimate1_Standard=data(:,3);
SE1=data(:,4);
Obs1=data(:,5);
Z1=data(:,6);
Pval1=data(:,7);
Pval1_Analysis=data(:,8);

Estimate2=data(:,9);
Estimate2_Standard=data(:,10);
SE2=data(:,11);
Obs2=data(:,12);
Z2=data(:,13);
Pval2=data(:,14);
Pval2_Analysis=data(:,15);

%Use Fisher Transformation to obtain approximately normal standardized
%estimates
Fisher=@(x) 0.5*log((1+x)./(1-x));
Estimate1_final=Fisher(Estimate1_Standard);
Estimate2_final=Fisher(Estimate2_Standard);

Pval1_final=Pval1_Analysis;
%Where Pval1_Analysis missig (censored in original dataset), fill in with
%pval collected from replication report.
Pval1_final(Pval1_final==0)=Pval1(Pval1_final==0);
%Impute Pval1 from Z statistics wherever still missing
Pval1_impute=2*(1-normcdf(Z1));
Pval1_final(Pval1_final==0)=Pval1_impute(Pval1_final==0);

Pval2_final=Pval2_Analysis;
%Where Pval1_Analysis missig (censored in original dataset), fill in with
%pval collected from replication report.
Pval2_final(Pval2_final==0)=Pval1(Pval2_final==0);
%Impute Pval1 from Z statistics wherever still missing
Pval2_impute=2*(2-normcdf(Z1));
Pval2_final(Pval2_final==0)=Pval1_impute(Pval2_final==0);

%Impute z-statistics from p-values and standardized effect sizes
Z1_final=norminv(1-Pval1_final/2);
Z2_final=sign(Estimate2_final).*norminv(1-Pval2_final/2);

%Fill in infinite Z statistics from raw Z statistics
Z1_final(Z1_final==inf)=Z1(Z1_final==inf);
Z2_final(Z2_final==inf)=Z2(Z2_final==inf);

SE1=abs(Estimate1_final)./abs(Z1_final);
SE2=abs(Estimate2_final)./abs(Z2_final);

output=[Estimate1_final SE1 Estimate2_final SE2 Z1_final];
%Read in list of names
names=readtable('../../Applications/EconExperiments/EconLabels.csv');
%Sort list of names
[output, idx]=sortrows(output,5);
sorted_names=names(idx,:);

%Drop column Z statistics (which included only for sorting purposes)
output(:,end)=[];

%Output data and names
csvwrite('../../Applications/EconExperiments/cleaned_econ_data.csv',output)
writetable(sorted_names(:,1),'../../Applications/EconExperiments/sorted_names.csv')


