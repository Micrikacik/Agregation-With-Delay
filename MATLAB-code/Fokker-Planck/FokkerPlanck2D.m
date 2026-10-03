function FokkerPlanck2D(N,M,int_r,tau)
%Fokker-Planck equation for aggregation on the 2D torus
%finite volumes, semi-implicit discretization in time
%(N,M) = number of grid points
%int_r = interaction radius
%tau = delay measured in dt

%random seed
rng(1);

L = N*M;

T = 1e5;
dt = 1e-3;
dx = 1/N;
dtxx = dt/dx^2;

WhichPlot = 100;
snapshots=0;



%Load the matrix A
fname = sprintf('A_%dx%d.mat',N,M);
load(fname);

A = A*dtxx;


%Load the distance matrix and create matrix W
fname = sprintf('distM_%dx%d.mat',N,M);
load(fname);


W = sparse(distM <= int_r);
W = W*L/sum(sum(W));

%we don't need distM any more, so we can release the memory
clear distM;

%Prepare video
%fname = sprintf('exp%d/exp_%dx%d_%1.2f.avi',expno,N,M,int_r);
%vidObj = VideoWriter(fname);
%open(vidObj);


%Sparse unit matrix of size L
SeyeL = sparse(1:L,1:L,1);


%Initial condition
rho = rand(L,1);
rho = rho / (dx^2*sum(rho));

rhoBuf = zeros(tau,L);
rhoBuf(1,:) = rho;
for t=2:tau
    rhoBuf(t,:) = rho;
    %rhoBuf(t,:) = (SeyeL + A) \ rhoBuf(t-1,:)';
    %pcolor(reshape(rhoBuf(t,:),N,M)); shading interp; clim([0 2]); getframe;
end

rho = rhoBuf(tau,:)';

%pause;

Fcount = 1;
%Solve in time
for t=1:T
    
    %take rhodelay from the buffer
    ttau = mod(t-1-tau,tau)+1;
    rhodelay = rhoBuf(ttau,:)';

    %store present rho to the buffer
    rhoBuf(ttau,:,:) = rho;

    Wrho = W*rhodelay;
    %min(Wrho)
    %max(Wrho)    
    
    %Plot
    if(~mod(t-1,WhichPlot))
        
        subplot(1,2,1);
        pcolor(reshape(rhodelay,N,M)); shading interp;
        clim([0 2]);
        colorbar;
        
        subplot(1,2,2);
        pcolor(reshape(rho,N,M)); shading interp;
        clim([0 2]);
        colorbar;
        
        getframe;
        %F(Fcount) = getframe(gcf);
        Fcount = Fcount+1;
        %writeVideo(vidObj,currFrame);
        
        if(snapshots)
            fname = sprintf('exp%d/exp%d_%dx%d_%1.2f_%d.fig',expno,expno,N,M,int_r,t);
            saveas(gcf,fname)
            
            fname = sprintf('exp%d/exp%d_%dx%d_%1.2f_%d.eps',expno,expno,N,M,int_r,t);
            saveas(gcf,fname,'eps')        
        end
            
        %if (t==1) pause; end;
    end
    
    
    %Calculate G
    %G=2-(Wrho>1);  %step diffusions
    F = exp(-1*Wrho);    %exponential diffusion
    G = (F.^2)/2;
    
    %Make one timestep
    rho = (SeyeL + A*sparse(1:L,1:L,G)) \ rho;
    %[rho,flag] = bicg(SeyeL + A*sparse(1:L,1:L,G),rho,1e-8,100);    
    %if(flag) break; end;
            
    %sum(rho)*dx^2
        
end


% Close the video file
%close(vidObj);
%fname = sprintf('exp%d/exp_%dx%d_%1.2f.avi',expno,N,M,int_r);
%movie2avi(F, fname);

end
