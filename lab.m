

inp = 'input\n>'
x = input(inp);

dis = abs(DFT(x));
plot(dis);

s = getLocalMax(dis);
[u,sss] = size(s);
s = 100*s;

k = 3 + sss;


s
MNK(x,k,s)

function r =  getCoef(x, t,k ,sv)
    if( x == 1)
        r = t*t*t;
    end
    if( x == 2)
        r = t*t;
    end
    if( x== 3)
        r = t;
    end
    
    if(x>3 && x<= k )
        r= sin( 2* pi *sv(x-3)*t);
    end
    
    if(x>k)
        r = 1;
    end
end


function s = getLocalMax(x)
    [un,n] = size(x);
    
    s = [];
    for i = 2:(n/2-1)
        if (x(i)>x(i-1)) && (x(i)>x(i+1))
            s = [s,(i-1)/(n-1)];
        end   
    end
end

function sl = MNK(val,k,sv)
    [unused,N] = size(val);
    T = 5;
    dt = 0.01;
    B = zeros(1,k+1);
    
    A = zeros(k+1,k+1);
    
    for i=1:k+1
        for j=1:k+1
            for t=1:N
                A(i,j) = A(i,j) + getCoef(i,(t-1)*dt,k,sv)*getCoef(j,(t-1)*dt,k,sv);
            end
        end
        
        for t=1:N
            B(i) = B(i) + getCoef(i,(t-1)*dt,k,sv)*val(t);
        
        end
    end
        
    sl = linsolve(transpose(A),transpose(B));
          
end
    
    
function c = DFT(x)

	[unused1,N] = size(x);
    c = zeros(1,N);
    for k=1:N
        for j = 1:N
            c(k)= c(k) + x(j)*exp(-1i*2*pi*(k-1)*(j-1)/N);
        end
        c(k)= c(k) / N;
    end 
end

function x = revFT(c)
	[unused1,N] = size(c);
    x = zeros(1,N);
    for k=1:N
        for j = 1:N
            x(k)= x(k) + c(j)*exp(1i*2*pi*(k-1)*(j-1)/N);
        end
    end 
end

