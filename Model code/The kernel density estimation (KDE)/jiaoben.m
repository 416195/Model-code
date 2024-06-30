clc
clear
cd('D:\train'); % File location.
X = xlsread('train.xlsx'); % File name.

begin_year = 2000; % Begin year.
xx = 'CGPE'; 
yy = 'Year'; 
zz = 'Kernel estimates';
title('')
[Max_Value,Min_Value,T,X] = N_density(begin_year,xx,yy,zz,X)