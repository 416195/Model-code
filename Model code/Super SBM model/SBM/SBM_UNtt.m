function [lambda,s_x,s_y,s_b,E_tt,E_t1t1]=SBM_UNtt(data,T,xs,ys,yb,rts,Super);

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
A1=zeros(1,n+m+s1+s2+1);b1=0;
LB1=[zeros(n+m+s1+s2,1);1e-10];UB=[];
LB2=[zeros(n+m+s1+s2-1,1);1e-10];
w=[];fval=[];

options=optimset('TolFun',1e-6,'MaxIter',1e+3);
if rts==0;
    for t=1:T;
        for i=1:n;
            f1=[zeros(1,n) -(1/m)*DX(i+(t-1)*n,:) zeros(1,s1) zeros(1,s2) 1];
            Aeq1=[X(:,n*(t-1)+1:n*t) eye(m) zeros(m,s1) zeros(m,s2) -X(:,i+(t-1)*n);
                Y(:,n*(t-1)+1:n*t) zeros(s1,m) -eye(s1) zeros(s1,s2) -Y(:,i+(t-1)*n);
                Z(:,n*(t-1)+1:n*t) zeros(s2,m) zeros(s2,s1) eye(s2) -Z(:,i+(t-1)*n);
                zeros(1,n) zeros(1,m) (1/(s1+s2))*DY(i+(t-1)*n,:) (1/(s1+s2))*DZ(i+(t-1)*n,:) 1];
            beq1=[zeros(m,1);
                zeros(s1,1);
                zeros(s2,1);
                1];
            try
                [w(:,i+(t-1)*n),fval1(:,i+(t-1)*n),exitflag1(:,i+(t-1)*n)]=linprog(f1,A1,b1,Aeq1,beq1,LB1,UB,[],options);
            catch
                w(:,i+(t-1)*n)=nan(n+m+s1+s2+1,1);
                fval1(:,i+(t-1)*n)=nan;
            end
            if Super==0;
                fval(:,i+(t-1)*n)=fval1(:,i+(t-1)*n);
            elseif Super==1;
                f2=[zeros(1,n-1) (1/m)*DX(i+(t-1)*n,:) zeros(1,s1) zeros(1,s2) 1];
                Aeq2=[zeros(1,n-1) zeros(1,m) -(1/(s1+s2))*DY(i+(t-1)*n,:) -(1/(s1+s2))*DZ(i+(t-1)*n,:) 1];
                beq2=1;
                X2=X(:,n*(t-1)+1:n*t);Y2=Y(:,n*(t-1)+1:n*t);Z2=Z(:,n*(t-1)+1:n*t);
                X2(:,i)=[];Y2(:,i)=[];Z2(:,i)=[];
                A2=[X2 -eye(m) zeros(m,s1) zeros(m,s2) -X(:,i+(t-1)*n);
                    -Y2 zeros(s1,m) -eye(s1) zeros(s1,s2) Y(:,i+(t-1)*n);
                    Z2 zeros(s2,m) zeros(s2,s1) -eye(s2) -Z(:,i+(t-1)*n)];
                b2=[zeros(m,1);
                    zeros(s1,1);
                    zeros(s2,1)];
                try
                    [w2(:,i+(t-1)*n),fval2(:,i+(t-1)*n),exitflag2(:,i+(t-1)*n)]=linprog(f2,A2,b2,Aeq2,beq2,LB2,UB,[],options);
                catch
                    w2(:,i+(t-1)*n)=nan(n+m+s1+s2,1);
                    fval2(:,i+(t-1)*n)=nan;
                end
                fval(:,i+(t-1)*n)=fval1(:,i+(t-1)*n)*fval2(:,i+(t-1)*n);
            end
        end
    end
elseif rts==1;
    for t=1:T;
        for i=1:n;
            f1=[zeros(1,n) -(1/m)*DX(i+(t-1)*n,:) zeros(1,s1) zeros(1,s2) 1];
            Aeq1=[X(:,n*(t-1)+1:n*t) eye(m) zeros(m,s1) zeros(m,s2) -X(:,i+(t-1)*n);
                Y(:,n*(t-1)+1:n*t) zeros(s1,m) -eye(s1) zeros(s1,s2) -Y(:,i+(t-1)*n);
                Z(:,n*(t-1)+1:n*t) zeros(s2,m) zeros(s2,s1) eye(s2) -Z(:,i+(t-1)*n);
                zeros(1,n) zeros(1,m) (1/(s1+s2))*DY(i+(t-1)*n,:) (1/(s1+s2))*DZ(i+(t-1)*n,:) 1;
                ones(1,n) zeros(1,m+s1+s2) -1];
            beq1=[zeros(m,1);
                zeros(s1,1);
                zeros(s2,1);
                1;0];
            try
                [w(:,i+(t-1)*n),fval1(:,i+(t-1)*n),exitflag1(:,i+(t-1)*n)]=linprog(f1,A1,b1,Aeq1,beq1,LB1,UB,[],options);
            catch
                w(:,i+(t-1)*n)=nan(n+m+s1+s2+1,1);
                fval1(:,i+(t-1)*n)=nan;
            end
            if Super==0;
                fval(:,i+(t-1)*n)=fval1(:,i+(t-1)*n);
            elseif Super==1;
                f2=[zeros(1,n-1) (1/m)*DX(i+(t-1)*n,:) zeros(1,s1) zeros(1,s2) 1];
                Aeq2=[zeros(1,n-1) zeros(1,m) -(1/(s1+s2))*DY(i+(t-1)*n,:) -(1/(s1+s2))*DZ(i+(t-1)*n,:) 1;
                    ones(1,n-1) zeros(1,m+s1+s2) -1];
                beq2=[1;0];
                X2=X(:,n*(t-1)+1:n*t);Y2=Y(:,n*(t-1)+1:n*t);Z2=Z(:,n*(t-1)+1:n*t);
                X2(:,i)=[];Y2(:,i)=[];Z2(:,i)=[];
                A2=[X2 -eye(m) zeros(m,s1) zeros(m,s2) -X(:,i+(t-1)*n);
                    -Y2 zeros(s1,m) -eye(s1) zeros(s1,s2) Y(:,i+(t-1)*n);
                    Z2 zeros(s2,m) zeros(s2,s1) -eye(s2) -Z(:,i+(t-1)*n)];
                b2=[zeros(m,1);
                    zeros(s1,1);
                    zeros(s2,1)];
                try
                    [w2(:,i+(t-1)*n),fval2(:,i+(t-1)*n),exitflag2(:,i+(t-1)*n)]=linprog(f2,A2,b2,Aeq2,beq2,LB2,UB,[],options);
                catch
                    w2(:,i+(t-1)*n)=nan(n+m+s1+s2,1);
                    fval2(:,i+(t-1)*n)=nan;
                end
                fval(:,i+(t-1)*n)=fval1(:,i+(t-1)*n)*fval2(:,i+(t-1)*n);
            end
        end
    end
end
   

t2=w(n+m+s1+s2+1,:)';
Lambda=w(1:n,:)';
S_x=w(n+1:n+m,:)';
S_y=w(n+m+1:n+m+s1,:)';
S_b=w(n+m+s1+1:n+m+s1+s2,:)';
lambda=[data(:,1:2) t2.^(-1).*Lambda];
s_x=[data(:,1:2) t2.^(-1).*S_x];
s_y=[data(:,1:2) t2.^(-1).*S_y];
s_b=[data(:,1:2) t2.^(-1).*S_b];
E_tt=[data(:,1:2) fval'];
E_t1t1=[data((n+1):end,1:2) fval(:,(n+1):end)'];
end




