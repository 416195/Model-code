clc
clear
sp_markov = xlsread('train.xlsx',1,'B2:V572');      %% Read the data in the range B2£ºJ62 from sheet1 of the Excel file named 'train'.
B= xlsread('train.xlsx',2,'B2:UZ572');     %%Read the data in the range B2-BJ62 from sheet2 of the Excel file named 'train' as the spatial weight matrix.

% Classification using the quartile method
Q1 = prctile(sp_markov(:),25);
Q2 = prctile(sp_markov(:),50);
Q3 = prctile(sp_markov(:),75);
tt=1;% Number of spanning years
p=tra_mark(sp_markov,tt,Q1,Q2,Q3);
pp=space_mark(sp_markov,B,tt,Q1,Q2,Q3);
results=[p;pp];  %The first four rows of the results are the traditional Markov results. Rows 5-8 represent the Low level; rows 9-12 represent the Medium level; rows 13-16 represent the Medium high level; and rows 17-20 represent the High level.
