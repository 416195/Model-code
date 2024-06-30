function Exresults = SBM_ex(T,info,results); 

fields = fieldnames(info);
nf = length(fields);

if nf > 0
    for i=1:nf
        if strcmp(fields{i},'model') model = info.model;
        elseif strcmp(fields{i},'rts') rts = info.rts;
        elseif strcmp(fields{i},'Super') Super = info.Super;
        elseif strcmp(fields{i},'Glob') Glob = info.Glob;
        elseif strcmp(fields{i},'Decomp') Decomp = info.Decomp;
        end
    end
end

if T==1;
    Glob=1;
elseif T<=0;
    error('SBM: wrong input number of T');
end

if rts==0;
    if Glob==0;
        title_E1={'id','year','Et_C(t)'};
        title_E2={'id','year','Et+1_C(t+1)','Et_C(t)','Et+1_C(t)','Et_C(t+1)','MI_c','EC_c','TC_c'}; 
        title_FGNZ={'id','year','MI_c','PEC','SEC','TC_c'};
        title_RD={'id','year','MI_c','PEC','SCH','PTC'};
        title_Zofio={'id','year','MI_c','PEC','SEC','PTC','STC'};    
    elseif Glob==1;
        title_E1={'id','year','EG_C(t)'};
        title_E2={'id','year','EG_C(t+1)','EG_C(t)','GMI_c'};
        title_PL={'id','year','MI_c','EC_c','BPC_c'};
    end
elseif rts==1;
    if Glob==0;
        title_E1={'id','year','Et_V(t)'};
        title_E2={'id','year','Et+1_V(t+1)','Et_V(t)','Et+1_V(t)','Et_V(t+1)','MI_v','EC_v','TC_v'};
    elseif Glob==1;
        title_E1={'id','year','EG_V(t)'};
        title_E2={'id','year','EG_V(t+1)','EG_V(t)','GMI_v'};
        title_PL={'id','year','MI_v','EC_v','BPC_v'};
    end
end

file_name1=[results.method '_eff'];
file_name2=[results.method '_Index'];
file_name3=[results.method '_Decomp'];
Exresults.method=results.method;

%% Efficiency
if rts==0;
    if Glob==0;
        E1=num2cell(results.Ec_tt);
    elseif Glob==1
        E1=num2cell(results.Ec_Gt);
    end
elseif rts==1;
    if Glob==0;
        E1=num2cell(results.Ev_tt);
    elseif Glob==1;
        E1=num2cell(results.Ev_Gt);
    end
end
E1_ex=[title_E1;E1];
xlswrite(file_name1,E1_ex,'eff');
Exresults.Efficiency=E1_ex;
%% lambda

Lambda=num2cell(results.lambda);
for i=1:(size(results.lambda,2)-2);
    title_lambda{1,i}=strcat('lambda',num2str(i))
end
title_lambda2=['id','year',title_lambda]
lambda_ex=[title_lambda2;Lambda];
try
xlswrite(file_name1,lambda_ex,'Lambda');
end
try
Exresults.Lambda=lambda_ex;
end
%% Slacks

for i=1:(size(results.s_x,2)-2);
    title_sx{1,i}=strcat('s_x',num2str(i))
end
for i=1:(size(results.s_y,2)-2);
    title_sy{1,i}=strcat('s_y',num2str(i))
end
if model==0;
    title_sxyz=['id','year',title_sx,title_sy];
    s_xyz=num2cell([results.s_x,results.s_y(:,3:end)]);
elseif model==1;
    for i=1:(size(results.s_z,2)-2);
        title_sz{1,i}=strcat('s_b',num2str(i))
    end
    title_sxyz=['id','year',title_sx,title_sy,title_sz];
    s_xyz=num2cell([results.s_x,results.s_y(:,3:end),results.s_z(:,3:end)]);
end
s_xyz_ex=[title_sxyz;s_xyz];
xlswrite(file_name1,s_xyz_ex,'Slacks');
Exresults.Slacks=s_xyz_ex;
%% Malmquist Index

if T>1;
    if rts==0;
        if Glob==0;
            E2=num2cell([results.Ec_t1t1,results.Ec_tt(1:size(results.Ec_t1t1,1),3),results.Ec_t1t(:,3),results.Ec_tt1(:,3),results.MI_c(:,3),results.EC_c(:,3),results.TC_c(:,3)]);
        elseif Glob==1;
            E2=num2cell([results.Ec_Gt1,results.Ec_Gt(1:size(results.Ec_Gt1,1),3),results.GMI_c(:,3)]);
        end
    elseif rts==1;
        if Glob==0;
            E2=num2cell([results.Ev_t1t1,results.Ev_tt(1:size(results.Ev_t1t1,1),3),results.Ev_t1t(:,3),results.Ev_tt1(:,3),results.MI_v(:,3),results.EC_v(:,3),results.TC_v(:,3)]);
        elseif Glob==1;
            E2=num2cell([results.Ev_Gt1,results.Ev_Gt(1:size(results.Ev_Gt1,1),3),results.GMI_v(:,3)]);
        end
    end
    E2_ex=[title_E2;E2];
    xlswrite(file_name2,E2_ex,'Index');
    Exresults.Productivity_Index=E2_ex;
end

%% Index Decomposition
if T>1;
    if Decomp==1
        if Glob==0;
            if rts==0;
                E_FGNZ=num2cell([results.MI_c,results.FGNZ_1994.PEC(:,3),results.FGNZ_1994.SEC(:,3),results.FGNZ_1994.TC_c(:,3)]);
                E_RD=num2cell([results.MI_c,results.RD_1997.PEC(:,3),results.RD_1997.SCH(:,3),results.RD_1997.PTC(:,3)]);
                E_Zofio=num2cell([results.MI_c,results.Zofio_2007.PEC(:,3),results.Zofio_2007.SEC(:,3),results.Zofio_2007.PTC(:,3),results.Zofio_2007.STC(:,3)]);
                FGNZ_ex=[title_FGNZ;E_FGNZ];
                RD_ex=[title_RD;E_RD];
                Zofio_ex=[title_Zofio;E_Zofio];
                xlswrite(file_name3,FGNZ_ex,'FGNZ(1994)');
                xlswrite(file_name3,RD_ex,'RD(1997)');
                xlswrite(file_name3,Zofio_ex,'Zofio(2007)');
                Decomposition.FGNZ_1994=FGNZ_ex;
                Decomposition.RD_1997=RD_ex;
                Decomposition.Zofio_2007=Zofio_ex;
                Exresults.Decomposition=Decomposition;
            end
        elseif Glob==1;
            if rts==0;
                E_PL=num2cell([results.GMI_c,results.PL_2005.EC_c(:,3),results.PL_2005.BPC_c(:,3)]);
            elseif rts==1;
                E_PL=num2cell([results.GMI_v,results.PL_2005.EC_v(:,3),results.PL_2005.BPC_v(:,3)]);
            end
            PL_ex=[title_PL;E_PL];
            xlswrite(file_name3,PL_ex,'PL(2005)');
            Decomposition.PL_2005=PL_ex;
            Exresults.Decomposition=Decomposition;
        end
    end
end
end


 
                
                
                
            
            
        

                    

                
                  
                





    

