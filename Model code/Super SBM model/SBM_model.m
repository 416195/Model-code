%% Parameter settings
T=21; % Number of years (Panel data T > 1,otherwise T = 1) ;
nx=6; % The number of input variables ;
ny=3; % The number of desirable output variables ;
nb=2; % The number of undesirable output variables ;
info.model=1; % 0 SBM model; 1 Un_SBM model ;
info.rts=1; % 0 Constant returns to scale (CRS); 1 Variable returns to scale (VRS) ;
info.Super=1; % 0 Excluding super-efficiency; 1 Including super-efficiency ;
info.Glob=0; % 0 DEA ; 1 Global DEA ; 
info.Decomp=1; % 0 No need for index decomposition; 1 Index decomposition needed ;
results=SBM_V2(data,T,nx,ny,nb,info); % Computational model results ;
results_Ex=SBM_ex(T,info,results); % Export computational results .

