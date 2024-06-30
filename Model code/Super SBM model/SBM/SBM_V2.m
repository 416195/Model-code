function results = SBM_V2(data,T,nx,ny,nz,info); 
time1 = 0; 
time2 = 0;
time3 = 0;
time4 = 0;

timet = clock; % start the clock for overall timing

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

if model==0 & rts==0 & Glob==0 & Super==0
    results.method='SBM_CRS'
elseif model==0 & rts==1 & Glob==0 & Super==0
    results.method='SBM_VRS'
elseif model==0 & rts==0 & Glob==0 & Super==1
    results.method='Super_SBM_CRS'
elseif model==0 & rts==1 & Glob==0 & Super==1
    results.method='Super_SBM_VRS'
elseif model==1 & rts==0 & Glob==0 & Super==0
    results.method='Un_SBM_CRS'
elseif model==1 & rts==1 & Glob==0 & Super==0
    results.method='Un_SBM_VRS'
elseif model==1 & rts==0 & Glob==0 & Super==1
    results.method='Un_Super_SBM_CRS'
elseif model==1 & rts==1 & Glob==0 & Super==1
    results.method='Un_Super_SBM_VRS'    
elseif model==0 & rts==0 & Glob==1 & Super==0
    results.method='Global_SBM_CRS'
elseif model==0 & rts==1 & Glob==1 & Super==0
    results.method='Global_SBM_VRS'
elseif model==0 & rts==0 & Glob==1 & Super==1
    results.method='Global_Super_SBM_CRS'
elseif model==0 & rts==1 & Glob==1 & Super==1
    results.method='Global_Super_SBM_VRS'
elseif model==1 & rts==0 & Glob==1 & Super==0
    results.method='Global_Un_SBM_CRS'
elseif model==1 & rts==1 & Glob==1 & Super==0
    results.method='Global_Un_SBM_VRS'
elseif model==1 & rts==0 & Glob==1 & Super==1
    results.method='Global_Un_Super_SBM_CRS'
elseif model==1 & rts==1 & Glob==1 & Super==1
    results.method='Global_Un_Super_SBM_VRS'
else
    error('SBM: wrong input number of info');
end    

if Glob==0;
    if model==0;
        s_z=[];
        [lambda,s_x,s_y,E_tt,E_t1t1]=SBM_DEtt(data,T,nx,ny,rts,Super);
        E_t1t = SBM_DEt1t(data,T,nx,ny,rts);
        E_tt1 = SBM_DEtt1(data,T,nx,ny,rts);
    elseif model==1;
        [lambda,s_x,s_y,s_z,E_tt,E_t1t1]=SBM_UNtt(data,T,nx,ny,nz,rts,Super);
        E_t1t = SBM_UNt1t(data,T,nx,ny,nz,rts);
        E_tt1 = SBM_UNtt1(data,T,nx,ny,nz,rts);
    end    
    EC=(E_tt(1:size(E_t1t1,1),3).^(-1)).*E_t1t1(:,3);
    TC=((EC.^(-1)).*((E_t1t(:,3).^(-1)).*E_tt1(:,3))).^0.5;
    MI=EC.*TC;
    EC=[E_t1t1(:,1:2) EC];
    TC=[E_t1t1(:,1:2) TC];
    MI=[E_t1t1(:,1:2) MI];
    if rts==0;
        results.Ec_tt=E_tt;
        results.Ec_t1t1=E_t1t1;
        results.Ec_t1t=E_t1t;
        results.Ec_tt1=E_tt1;
        results.EC_c=EC;
        results.TC_c=TC;
        results.MI_c=MI;
    elseif rts==1
        results.Ev_tt=E_tt;
        results.Ev_t1t1=E_t1t1;
        results.Ev_t1t=E_t1t;
        results.Ev_tt1=E_tt1;
        results.EC_v=EC;
        results.TC_v=TC;
        results.MI_v=MI;
    end
elseif Glob==1;
    if model==0;
        s_z=[];
        [lambda,s_x,s_y,E_Gt,E_Gt1]=SBM_DEGt(data,T,nx,ny,rts,Super);
    elseif model==1;
        [lambda,s_x,s_y,s_z,E_Gt,E_Gt1]=SBM_UNGt(data,T,nx,ny,nz,rts,Super);
    end
    GMI=(E_Gt(1:size(E_Gt1,1),3).^(-1)).*E_Gt1(:,3);
    GMI=[E_Gt1(:,1:2) GMI];   
    if rts==0;
        results.Ec_Gt=E_Gt;
        results.Ec_Gt1=E_Gt1;
        results.GMI_c=GMI;
    elseif rts==1;
        results.Ev_Gt=E_Gt;
        results.Ev_Gt1=E_Gt1;
        results.GMI_v=GMI;
    end
else
    error('SBM: wrong input number of info.model');
end
results.lambda=lambda;
results.s_x=s_x;
results.s_y=s_y;
results.s_z=s_z;

%% Index_Decomposition

if Decomp==0;
    disp('No Index Decomposition!!!')
elseif Decomp==1;
    if Glob==0;
        if rts==0;
            if model==0;
                [lambda,s_x,s_y,Ev_tt,Ev_t1t1]=SBM_DEtt(data,T,nx,ny,1,Super);
                Ev_t1t = SBM_DEt1t(data,T,nx,ny,1);
                Ev_tt1 = SBM_DEtt1(data,T,nx,ny,1);
            elseif model==1;
                [lambda,s_x,s_y,s_z,Ev_tt,Ev_t1t1]=SBM_UNtt(data,T,nx,ny,nz,1,Super);
                Ev_t1t = SBM_UNt1t(data,T,nx,ny,nz,1);
                Ev_tt1 = SBM_UNtt1(data,T,nx,ny,nz,1);
            end
            PEC=(Ev_tt(1:size(Ev_t1t1,1),3).^(-1)).*Ev_t1t1(:,3);
            PTC=((PEC.^(-1)).*((Ev_t1t(:,3).^(-1)).*Ev_tt1(:,3))).^0.5;
            SEC=(((Ev_tt(1:size(Ev_t1t1,1),3).^(-1)).*E_tt(1:size(Ev_t1t1,1),3)).^(-1)).*((Ev_t1t1(:,3).^(-1)).*E_t1t1(:,3));
            SCH=(((((Ev_t1t(:,3).^(-1)).*E_t1t(:,3)).^(-1)).*((Ev_tt1(:,3).^(-1)).*E_tt1(:,3))).*SEC).^0.5;
            STC=((SEC.^(-1)).*((((Ev_t1t(:,3).^(-1)).*E_t1t(:,3)).^(-1)).*((Ev_tt1(:,3).^(-1)).*E_tt1(:,3)))).^0.5;
            PEC=[Ev_t1t1(:,1:2) PEC];
            PTC=[Ev_t1t1(:,1:2) PTC];
            SEC=[Ev_t1t1(:,1:2) SEC];
            SCH=[Ev_t1t1(:,1:2) SCH];
            STC=[Ev_t1t1(:,1:2) STC];
            FGNZ_1994.PEC=PEC;
            FGNZ_1994.SEC=SEC;
            FGNZ_1994.TC_c=results.TC_c;
            RD_1997.PEC=PEC;
            RD_1997.SCH=SCH;
            RD_1997.PTC=PTC;
            Zofio_2007.PEC=PEC;
            Zofio_2007.SEC=SEC;
            Zofio_2007.PTC=PTC;
            Zofio_2007.STC=STC;
            results.FGNZ_1994=FGNZ_1994;
            results.RD_1997=RD_1997;
            results.Zofio_2007=Zofio_2007;
        elseif rts==1;
            disp('No Index Decomposition!!!')
        end
    elseif Glob==1;
        if model==0;
            [lambda,s_x,s_y,E_tt,E_t1t1]=SBM_DEtt(data,T,nx,ny,rts,Super);
        elseif model==1;
            [lambda,s_x,s_y,s_z,E_tt,E_t1t1]=SBM_UNtt(data,T,nx,ny,nz,rts,Super);
        end
        EC=(E_tt(1:size(E_t1t1,1),3).^(-1)).*E_t1t1(:,3);
        BPC=(EC.^(-1)).*GMI(:,3);
        EC=[E_t1t1(:,1:2) EC];
        BPC=[E_t1t1(:,1:2) BPC];
        if rts==0;
            PL_2005.EC_c=EC;
            PL_2005.BPC_c=BPC;
        elseif rts==1;
            PL_2005.EC_v=EC;
            PL_2005.BPC_v=BPC;      
        end    
            results.PL_2005=PL_2005;
    end
else
    error('SBM: wrong input number of info.model');
end
        
end