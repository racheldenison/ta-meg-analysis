% rm_anova_tests.m

%% 3-way anova example
% generate random data for the example
alpha_power = randn(24,8); % subjects x conditions

% Create a table storing the respones
varNames = {'Y1','Y2','Y3','Y4','Y5','Y6','Y7','Y8'}; % condition names
t = array2table(alpha_power,'VariableNames',varNames);

% Create a table reflecting the within subject factors 'TestCond', 'Attention', and 'TMS' and their levels
factorNames = {'TestCond','Attention','TMS'};
within = table({'M';'M';'M';'M';'N';'N';'N';'N'},{'A';'A';'V';'V';'A';'A';'V';'V'},{'T';'S';'T';'S';'T';'S';'T';'S'},'VariableNames',factorNames);

% fit the repeated measures model
rm = fitrm(t,'Y1-Y8~1','WithinDesign',within);

% run my repeated measures anova hererm
[ranovatbl] = ranova(rm, 'WithinModel','TestCond*Attention*TMS');


%% 2-way anova example
% generate random data for the example
alpha_power = randn(24,4);

% Create a table storing the respones
varNames = {'Y1','Y2','Y3','Y4'};
t = array2table(alpha_power,'VariableNames',varNames);

% Create a table reflecting the within subject factors 'TestCond', 'Attention', and 'TMS' and their levels
factorNames = {'TestCond','Attention'};
within = table({'M';'M';'N';'N'},{'A';'V';'A';'V'},'VariableNames',factorNames);

% fit the repeated measures model
rm = fitrm(t,'Y1-Y4~1','WithinDesign',within);

% run my repeated measures anova here
[ranovatbl] = ranova(rm, 'WithinModel','TestCond*Attention');

