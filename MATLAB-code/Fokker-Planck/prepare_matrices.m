clearvars

N = 100;
M = 100;
L = N*M;

dx = 1/N;

%The periodic coordinate functions
fi=@(i) mod(i-1,N)+1;
fj=@(j) mod(j-1,M)+1;
fl=@(l) mod(l-1,M*N)+1;

lij=@(i,j) fi(i)+(fj(j)-1)*N;

jl=@(l) ceil(fl(l)/N);
il=@(l) fl(l) - (jl(l)-1)*N;

disti=@(i1,i2) min(abs(fi(i1)-fi(i2)),N-abs(fi(i1)-fi(i2)));
distj=@(j1,j2) min(abs(fj(j1)-fj(j2)),M-abs(fj(j1)-fj(j2)));


%Prepare the matrix A
SeyeL = sparse(1:L,1:L,1);
A = 4*SeyeL;

for i=1:N
    for j=1:M
        A(lij(i,j),lij(i+1,j)) = -1;
        A(lij(i,j),lij(i-1,j)) = -1;
        A(lij(i,j),lij(i,j+1)) = -1;
        A(lij(i,j),lij(i,j-1)) = -1;
    end
end


%Save the matrix A
fname = sprintf('A_%dx%d.mat',N,M);
save(fname,'A');

clear A

%Prepare the distance matrix
distM = zeros(L,L);

for l=1:L
    for k=l:L
        distM(l,k) = sqrt((disti(il(l),il(k))^2 + distj(jl(l),jl(k))^2))*dx;
        distM(k,l) = distM(l,k);
    end
end

%Save the matrix W
fname = sprintf('distM_%dx%d.mat',N,M);
save(fname,'distM');
