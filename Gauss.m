function X=Gauss(A,b)

[n,m] = size(A);

if(n~=m)
    disp('matrix should be square');
    return;
end

B=[A b];

X=zeros(size(b));
for i=1:n-1
    [p,q] = max(B(i:n,i));
    if(p==0)
        disp('sol not defined');
        return;
    end
    %a=B(i,i); 
    %X(i,1)=p;
    l=q+i-1;
    if(l~=i)
        B([l,i],:) = B([i,l],:);
    end
    for j=i+1:n
        a=B(j,i)/p;
        B(j,i:n+1)=B(j,i:n+1)-a*B(i,i:n+1);
    end
    
end

X(n,1)=B(n,n+1)/B(n,n);
for i=n-1:-1:1
    sum=0;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
    for j=i+1:n
        
        sum=sum+B(i,j)*X(j,1);
    end
    X(i,1)=(B(i,n+1)-sum)/B(i,i);
end
