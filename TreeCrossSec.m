function VV = TreeCrossSec(Sp, res_thres, Forder)

% TreeCrossSec: 2D cross-section reconstruction with free-form curves
% 
% Usage:
% VV = TreeCrossSec(Sp, res_thres, Forder)
% 
% Input:
% Sp: nxm point cloud (unit in m). The dimension m should be at least 2.
% optional:
% res_thres: residual tolerance e.g., 0.01 (cm)
% Forder: Fourier series order (1 - 8)
%
% Output:
% VV: a structure variable contains diameter estimations from various
% methods, and reconstructed points in ModelX&ModelY. Usually FouCirc gives
% the best diameter estimation.
% 
% 1. Arc length extrapolation
% VV.Fourier
% 2. Arc length + Circle
% VV.FouCirc
% VV.FouCirc_cent_x 
% VV.FouCirc_cent_y
% 3. Refined direct Circle fitting
% VV.RefinedCircle
% 4. Direct Circle fitting
% VV.DirectCircle
% 5. modeled points
% VV.ModelX
% VV.ModelY
% 6. parameters
% VV.Forder
% VV.Thres
%   
% If you use this code in your own work, please cite the following paper:
% Wang, D., Kankare, V., Puttonen, E., Hollaus, M., & Pfeifer, N. (2017).
% Reconstructing stem cross section shapes from terrestrial laser scanning.
% IEEE Geoscience and Remote Sensing Letters, 14(2), 272-276.
%
% Di Wang
% di.wang@aalto.fi


%% arguments
if nargin<1
    disp('function needs at least a point list (nx2) as input !')
    VV = [];
    return;
elseif (nargin < 2)
    res_thres = 0.02;
end
if (nargin <3)
    Forder = 5;
end

%%
Sp = Sp(:,1:2);
Ori_Sp = Sp;%backup
%% downsample to ~2k points to increase speed
if size(Sp,1)>2000
    rnum = 2000;
else
    rnum = size(Sp,1);
end
randindx = randperm(size(Sp,1),rnum);
Sp = Sp(randindx,:);

if size(Sp,1)>=10
    %% find initial parameters Tp for a rough circle with a fast inscribing method
    xx = min(Sp(:,1)):0.02:max(Sp(:,1));
    yy = min(Sp(:,2)):0.02:max(Sp(:,2));
    
    [xx,yy] = meshgrid(xx,yy);
    pot_loc = [xx(:),yy(:)];
    
    N = nan(size(pot_loc,1),1);
    for i = 1:size(pot_loc,1)
        dist = sqrt((Sp(:,1) - pot_loc(i,1)).^2 + (Sp(:,2) - pot_loc(i,2)).^2);
        N(i) = sum(abs(dist - prctile(dist,5))<0.02);
    end
    ia = find(N == max(N),1,'first');
    L1M = pot_loc(ia,:);
    
    r_init = median(sqrt(sum((Sp - repmat(L1M,size(Sp,1),1)).^2,2)));
    Tp = [L1M,r_init];
    
    %     visual debug
    %     figure
    %     scatter(Sp(:,1),Sp(:,2),2, 'fill')
    %     hold on
    %     th = 0:pi/50:2*pi;
    %     xunit = Tp(3) * cos(th) + Tp(1);
    %     yunit = Tp(3) * sin(th) + Tp(2);
    %     plot(xunit, yunit,'b','linewidth',2);
    %     hold off
    %     axis equal
    
    %% centered by rough center
    x = Sp(:,1) - Tp(1);
    y = Sp(:,2) - Tp(2);
    
    %% transform to polar corrdinates
    [THETA,RHO] = cart2pol(x,y);
    
    %     visual debug
    %     figure
    %     scatter(THETA,RHO,4,'fill')
    
    %% filter - 1
    T = [];
    for j = -3.2:0.02:3.2
        a = find(THETA>=j&THETA<j+0.02);
        if isempty(a)~=1
            c = RHO(a);
            d = c - min(c)<=res_thres*2;
            T = [T;a(d)];
        end
    end
    THETA = THETA(T);
    RHO = RHO(T);
    
    %% duplicate to 4*pi domain
    [~,bb] = sort(THETA);
    
    Sa = THETA(bb);
    Sb = RHO(bb);
    
    SbC = [Sb;Sb];
    SaC = [Sa;Sa+2*pi];
    
    %% Fit robust fourier series // Fourier approximation
    res=1;
    reqp=20;
    
    % Set up fittype and options.
    ftype = ['fourier' num2str(Forder)];
    ft = fittype( ftype );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    %opts.Display = 'Off';
    opts.Lower = [-Inf*ones(1,2*Forder+1) 1];
    opts.Robust = 'Bisquare';
    %     opts.Robust = 'off';
    opts.StartPoint = [zeros*ones(1,2*Forder+1) 1];
    opts.Upper = [Inf*ones(1,2*Forder+1) 1];
    
    while max(res)>res_thres && reqp>=19
        [xData, yData] = prepareCurveData( SaC, SbC );
        %% Fit model to data.
        [foo, ~] = fit( xData, yData, ft, opts );
        yy = feval(foo,SaC);
        res = SbC - yy;
        %% exclude outliers iteratively
        aa = find(res==max(res));
        SaC(aa) = [];
        SbC(aa) = [];
        
        reqp = length(SaC);
        
        %         hold on
        %         plot(foo)
        %         legend off
        %         pause(.01)
    end
    
    %% find approximate effective domain
    dis_thres = 0.1;
    
    dd = acos((2*(Tp(3).^2)-dis_thres^2)/2/(Tp(3).^2));
    [~, C, ~] = dbscan2105([SaC,SbC]', dd, 3);
    CL = cellfun(@length, C);
    Icl = find(CL==max(CL),1,'first');
    Ind_c = cell2mat(C(Icl));
    
    %% effective domain Sa_i
    Sa_i = SaC(Ind_c);
    Sb_i = SbC(Ind_c);
    yyy = feval(foo,Sa_i);
    
    %% transform fitted data -> Precise center location and radius
    [xa,ya] = pol2cart(Sa_i,yyy);
    Par = CircleFitByTaubin([xa,ya]);
    %% transform original domain data -> Direct circle fitting
    [~,ff] = sort(Sa_i);
    [xc,yc] = pol2cart(Sa_i(ff),Sb_i(ff));
    Par_d = CircleFitByTaubin([xc,yc]);
    
    a = max(Sa_i)-min(Sa_i);
    
    %% Arc length
    % calculate arc length in poalr CS, then convert to diameter
    % symbolic expression of foo
    sf = str2sym(formula( foo ));
    cn = coeffnames( foo);
    cv = coeffvalues( foo );
    for i = 1:length( cv )
        sf = subs( sf, cn{i}, cv(i) );
    end
    % Differentiate s
    syms x
    sf(x) = sf;
    ds = diff(sf,x);
    %syms x r
    arcL = sqrt(sf^2+ds^2);
    % to function handle
    g = matlabFunction(arcL);
    % perimeter
    q = integral(g,min(Sa_i),max(Sa_i));
    
    %%
    % ****************************************************************
    % Diamter results with various methods
    % ****************************************************************
    
    %% 1. Arc length extrapolation
    VV.Fourier = q*2*pi/a/2/pi;
    %% 2. Arc length + Circle
    VV.FouCirc = q/2/pi+Par(3)*(2*pi-a)/2/pi;
    VV.FouCirc_cent_x = Par(1) + Tp(1);
    VV.FouCirc_cent_y = Par(2) + Tp(2);
    %% 3. Refined direct Circle fitting
    VV.RefinedCircle = Par_d(3);
    %% 4. Direct Circle fitting
    direct_circle = CircleFitByTaubin(Ori_Sp);
    VV.DirectCircle = direct_circle(3);
    %% modeled points
    VV.ModelX = xa + Tp(1);
    VV.ModelY = ya + Tp(2);
    %% parameters
    VV.Forder = Forder;
    VV.Thres = res_thres;
    
    %% display results
    figure
    hold all
    scatter(Ori_Sp(:,1),Ori_Sp(:,2),4, 'fill','MarkerEdgeColor',[0.5 .5 .5],'MarkerFaceColor',[0.5 .5 .5])
    
    th = 0:pi/50:2*pi;
    
    % direct circle
    xunit = direct_circle(3) * cos(th) + direct_circle(1);
    yunit = direct_circle(3) * sin(th) + direct_circle(2);
    plot(xunit, yunit,'b','linewidth',2);
    % refined circle
    xunit = Par_d(3) * cos(th) + (Par_d(1)+Tp(1));
    yunit = Par_d(3) * sin(th) + (Par_d(2)+Tp(2));
    plot(xunit, yunit,'g','linewidth',2);
    % fouier fitting
    plot(VV.ModelX,VV.ModelY,'r','linewidth',2);
    plot([VV.ModelX(end);VV.ModelX(1)],[VV.ModelY(end);VV.ModelY(1)],'--r','linewidth',2);
    
    hold off
    
    lgd = legend({'Points','Direct Circle fit','Refined circle fit','Fourier fit'},'Location','northoutside');
    legend boxoff
    lgd.NumColumns = 4;
    box on
    axis equal
else
    %% 1. Arc length extrapolation
    VV.Fourier = nan;
    %% 2. Arc length + Circle
    VV.FouCirc = nan;
    VV.FouCirc_cent_x = nan;
    VV.FouCirc_cent_y = nan;
    %% 3. Refined direct Circle fitting
    VV.RefinedCircle = nan;
    %% 4. Direct Circle fitting
    VV.DirectCircle = nan;
    %% modeled points
    VV.ModelX = nan;
    VV.ModelY = nan;
    %% parameters
    VV.Forder = Forder;
    VV.Thres = res_thres;
end
end

% end of main function

%% ===========================================
function [W, C, ptsC] = dbscan2105(P, E, minPts)

[~, Npts] = size(P);

%neighbours=dis(P',E);
MdlKDT = KDTreeSearcher(P');
neighbours= rangesearch(MdlKDT,P',E);

W={};
ptsC  = zeros(Npts,1);
C     = {};
Nc    = 0;               % Cluster counter.
Pvisit = zeros(Npts,1);  % Array to keep track of points that have been visited.

for n = 1:Npts
    if ~Pvisit(n)                            % If this point not visited yet
        Pvisit(n) = 1;                       % mark as visited
        
        neighbourPts=cell2mat(neighbours(n));
        
        if length(neighbourPts) <= minPts-1  % Not enough points to form a cluster
            ptsC(n) = 0;                    % Mark point n as noise.
            
        else                % Form a cluster...
            Nc = Nc + 1;    % Increment number of clusters and process
            % neighbourhood.
            
            C{Nc} = [n];    % Initialise cluster Nc with point n
            ptsC(n) = Nc;   % and mark point n as being a member of cluster Nc.
            
            ind = 1;        % Initialise index into neighbourPts array.
            
            % For each point P' in neighbourPts ...
            while ind <= length(neighbourPts)
                
                nb = neighbourPts(ind);
                
                if ~Pvisit(nb)        % If this neighbour has not been visited
                    Pvisit(nb) = 1;   % mark it as visited.
                    
                    % Find the neighbours of this neighbour and if it has
                    % enough neighbours add them to the neighbourPts list
                    %neighbourPtsP = regionQuery(P, nb, E);
                    %nei2= rangesearch(KDtree,P(:,nb)',E);
                    neighbourPtsP = cell2mat(neighbours(nb));
                    
                    if length(neighbourPtsP) >= minPts
                        neighbourPts = [neighbourPts  neighbourPtsP];
                        %neighbourPts = unique(neighbourPts);
                    end
                end
                
                % If this neighbour nb not yet a member of any cluster add it
                % to this cluster.
                if ~ptsC(nb)
                    C{Nc} = [C{Nc} nb];
                    ptsC(nb) = Nc;
                end
                
                ind = ind + 1;  % Increment neighbour point index and process
                % next neighbour
            end
        end
    end
end

for i=1:length(C)
    
    W(i) ={P(:,cell2mat(C(i)))'};
end
end % of dbscan

%% ===================================================
function Par = CircleFitByTaubin(XY)

%--------------------------------------------------------------------------
%
%     Circle fit by Taubin
%      G. Taubin, "Estimation Of Planar Curves, Surfaces And Nonplanar
%                  Space Curves Defined By Implicit Equations, With
%                  Applications To Edge And Range Image Segmentation",
%      IEEE Trans. PAMI, Vol. 13, pages 1115-1138, (1991)
%
%     Input:  XY(n,2) is the array of coordinates of n points x(i)=XY(i,1), y(i)=XY(i,2)
%
%     Output: Par = [a b R] is the fitting circle:
%                           center (a,b) and radius R
%
%     Note: this fit does not use built-in matrix functions (except "mean"),
%           so it can be easily programmed in any programming language
%
%--------------------------------------------------------------------------

n = size(XY,1);      % number of data points

centroid = mean(XY);   % the centroid of the data set

%     computing moments (note: all moments will be normed, i.e. divided by n)

Mxx = 0; Myy = 0; Mxy = 0; Mxz = 0; Myz = 0; Mzz = 0;

for i=1:n
    Xi = XY(i,1) - centroid(1);  %  centering data
    Yi = XY(i,2) - centroid(2);  %  centering data
    Zi = Xi*Xi + Yi*Yi;
    Mxy = Mxy + Xi*Yi;
    Mxx = Mxx + Xi*Xi;
    Myy = Myy + Yi*Yi;
    Mxz = Mxz + Xi*Zi;
    Myz = Myz + Yi*Zi;
    Mzz = Mzz + Zi*Zi;
end

Mxx = Mxx/n;
Myy = Myy/n;
Mxy = Mxy/n;
Mxz = Mxz/n;
Myz = Myz/n;
Mzz = Mzz/n;

%    computing the coefficients of the characteristic polynomial

Mz = Mxx + Myy;
Cov_xy = Mxx*Myy - Mxy*Mxy;
A3 = 4*Mz;
A2 = -3*Mz*Mz - Mzz;
A1 = Mzz*Mz + 4*Cov_xy*Mz - Mxz*Mxz - Myz*Myz - Mz*Mz*Mz;
A0 = Mxz*Mxz*Myy + Myz*Myz*Mxx - Mzz*Cov_xy - 2*Mxz*Myz*Mxy + Mz*Mz*Cov_xy;
A22 = A2 + A2;
A33 = A3 + A3 + A3;

xnew = 0;
ynew = 1e+20;
epsilon = 1e-12;
IterMax = 20;

% Newton's method starting at x=0

for iter=1:IterMax
    yold = ynew;
    ynew = A0 + xnew*(A1 + xnew*(A2 + xnew*A3));
    if abs(ynew) > abs(yold)
        disp('Newton-Taubin goes wrong direction: |ynew| > |yold|');
        xnew = 0;
        break;
    end
    Dy = A1 + xnew*(A22 + xnew*A33);
    xold = xnew;
    xnew = xold - ynew/Dy;
    if (abs((xnew-xold)/xnew) < epsilon), break, end
    if (iter >= IterMax)
        disp('Newton-Taubin will not converge');
        xnew = 0;
    end
    if (xnew<0.)
        fprintf(1,'Newton-Taubin negative root:  x=%f\n',xnew);
        xnew = 0;
    end
end

%  computing the circle parameters

DET = xnew*xnew - xnew*Mz + Cov_xy;
Center = [Mxz*(Myy-xnew)-Myz*Mxy , Myz*(Mxx-xnew)-Mxz*Mxy]/DET/2;

Par = [Center+centroid , sqrt(Center*Center'+Mz)];

end    %    CircleFitByTaubin
