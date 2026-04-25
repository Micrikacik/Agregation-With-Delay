%calculate correlations between rho(t) and rho(t-tau)
%taking into account only values of rho <= Lval

Lval = 30;

%number of files
Nfiles = 15;
rho = zeros(1,Nfiles);

for K=1:Nfiles

    tau = (K-1)*30;

    %load data
    filename = sprintf('MCData_Reaction_2D_%d.mat', tau);
    load(filename);

    %some results are of size 401x401 with first row/column zero
    % - these we cut out
    L=Lval;
    if size(results{1}.thetaOccur) == 401
        L=Lval+1;
    end

%we need to have at least one nonzero element in each row/column
A=eye(L);

for k=1:size(results,1)
    A=A+results{k}.thetaOccur(1:L,1:L);
end

if size(results{1}.thetaOccur) == 401
    A = A(2:end,2:end);
end


P = A / sum(A(:));
N = size(A,1);
[u,v] = ndgrid(1:N, 1:N);

Eu = sum(u(:) .* P(:));
Ev = sum(v(:) .* P(:));

Covuv = sum((u(:)-Eu).*(v(:)-Ev).*P(:));
Varu = sum((u(:)-Eu).^2 .* P(:));
Varv = sum((v(:)-Ev).^2 .* P(:));

rho(K) = Covuv / sqrt(Varu * Varv);

end

plot(rho,'*');
