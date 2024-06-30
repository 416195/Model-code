function eff = SBM_UNtt1(data,T,xs,ys,yb,rts);

data=sortrows(data,[2 1]);
data2=data(:,3:end);
for i=1:size(data2,1);
    for j=1:size(data2,2);
        Ddata(i,j)=1/data2(i,j);
    end
end

Ddata(isnan(Ddata))=0;
Ddata(Ddata==inf)=0;
X=data2(:,1:xs)';
Y=data2(:,1+xs:xs+ys)';
Z=data2(:,1+xs+ys:xs+ys+yb)';
DX=Ddata(:,1:xs);
DY=Ddata(:,1+xs:xs+ys);
DZ=Ddata(:,1+xs+ys:xs+ys+yb);
NT=size(X',1);m=size(X,1);s1=size(Y,1);s2=size(Z,1);

n=NT/T;

LB=[zeros(n+m+s1+s2,1);1e-10];UB=[];
w=[];fval=[];
A1=zeros(1,n+m+s1+s2+1);b1=0;
options=optimset('TolFun',1e-6,'MaxIter',1e+3);
if rts==0;
    for t=1:(T-1);
        for i=1:n;
            f1=[zeros(1,n) -(1/m)*DX(i+t*n,:) zeros(1,s1) zeros(1,s2) 1];
            Aeq1=[X(:,n*(t-1)+1:n*t) eye(m) zeros(m,s1) zeros(m,s2) -X(:,i+t*n);
                Y(:,n*(t-1)+1:n*t) zeros(s1,m) -eye(s1) zeros(s1,s2) -Y(:,i+t*n);
                Z(:,n*(t-1)+1:n*t) zeros(s2,m) zeros(s2,s1) eye(s2) -Z(:,i+t*n);
                zeros(1,n) zeros(1,m) (1/(s1+s2))*DY(i+t*n,:) (1/(s1+s2))*DZ(i+t*n,:) 1];
            beq1=[zeros(m,1);
                zeros(s1,1);
                zeros(s2,1);
                1];
            f2=[zeros(1,n) (1/m)*DX(i+t*n,:) zeros(1,s1) zeros(1,s2) 1];
            Aeq2=[zeros(1,n) zeros(1,m) -(1/(s1+s2))*DY(i+t*n,:) -(1/(s1+s2))*DZ(i+t*n,:) 1];
            beq2=1;
            A2=[X(:,n*(t-1)+1:n*t) -eye(m) zeros(m,s1) zeros(m,s2) -X(:,i+t*n);
                -Y(:,n*(t-1)+1:n*t) zeros(s1,m) -eye(s1) zeros(s1,s2) Y(:,i+t*n);
                Z(:,n*(t-1)+1:n*t) zeros(s2,m) zeros(s2,s1) -eye(s2) -Z(:,i+t*n)];
            b2=[zeros(m,1);
                zeros(s1,1);
                zeros(s2,1)];
            try
                [w(:,i+(t-1)*n),fval(:,i+(t-1)*n),exitflag(:,i+(t-1)*n)]=linprog(f1,A1,b1,Aeq1,beq1,LB,UB,[],options);
            catch
                try
                    [w(:,i+(t-1)*n),fval(:,i+(t-1)*n),exitflag(:,i+(t-1)*n)]=linprog(f2,A2,b2,Aeq2,beq2,LB,UB,[],options);
                catch
                    w(:,i+(t-1)*n)=nan(n+m+s1+s2+1,1);
                    fval(:,i+(t-1)*n)=nan;
                end
            end
        end
    end
elseif rts==1;
    for t=1:(T-1);
        for i=1:n;
            f1=[zeros(1,n) -(1/m)*DX(i+t*n,:) zeros(1,s1) zeros(1,s2) 1];
            Aeq1=[X(:,n*(t-1)+1:n*t) eye(m) zeros(m,s1) zeros(m,s2) -X(:,i+t*n);
                Y(:,n*(t-1)+1:n*t) zeros(s1,m) -eye(s1) zeros(s1,s2) -Y(:,i+t*n);
                Z(:,n*(t-1)+1:n*t) zeros(s2,m) zeros(s2,s1) eye(s2) -Z(:,i+t*n);
                zeros(1,n) zeros(1,m) (1/(s1+s2))*DY(i+t*n,:) (1/(s1+s2))*DZ(i+t*n,:) 1;
                ones(1,n) zeros(1,m+s1+s2) -1];
            beq1=[zeros(m,1);
                zeros(s1,1);
                zeros(s2,1);
                1;0];
            f2=[zeros(1,n) (1/m)*DX(i+t*n,:) zeros(1,s1) zeros(1,s2) 1];
            Aeq2=[zeros(1,n) zeros(1,m) -(1/(s1+s2))*DY(i+t*n,:) -(1/(s1+s2))*DZ(i+t*n,:) 1;
                ones(1,n) zeros(1,m+s1+s2) -1];
            beq2=[1;0];
            A2=[X(:,n*(t-1)+1:n*t) -eye(m) zeros(m,s1) zeros(m,s2) -X(:,i+t*n);
                -Y(:,n*(t-1)+1:n*t) zeros(s1,m) -eye(s1) zeros(s1,s2) Y(:,i+t*n);
                Z(:,n*(t-1)+1:n*t) zeros(s2,m) zeros(s2,s1) -eye(s2) -Z(:,i+t*n)];
            b2=[zeros(m,1);
                zeros(s1,1);
                zeros(s2,1)];
            try
                [w(:,i+(t-1)*n),fval(:,i+(t-1)*n),exitflag(:,i+(t-1)*n)]=linprog(f1,A1,b1,Aeq1,beq1,LB,UB,[],options);
            catch
                try
                    [w(:,i+(t-1)*n),fval(:,i+(t-1)*n),exitflag(:,i+(t-1)*n)]=linprog(f2,A2,b2,Aeq2,beq2,LB,UB,[],options);
                catch
                    w(:,i+(t-1)*n)=nan(n+m+s1+s2+1,1);
                    fval(:,i+(t-1)*n)=nan;
                end
            end
        end
    end
end
   
eff=[data(n+1:end,1:2) fval'];
end




