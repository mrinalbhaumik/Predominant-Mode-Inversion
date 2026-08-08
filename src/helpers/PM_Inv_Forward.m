function[v, Domi_v, amp_v1] = PM_Inv_Forward(vs, nu, rho,...,
    h, w, HTLM, Measured_component)
% Forward model of the layered elastic half-space: HTLM/PML discretization,
% per-frequency eigen-solution, and modal/dominant dispersion extraction.
% Author      : Mrinal Bhaumik
% Affiliation : Utah State University, 2026

%% Orient all inputs as row vectors --------------------------------------------------------------
if ~isrow(vs)
    vs = vs';
end
if ~isrow(nu)
    nu = nu';
end
if ~isrow(rho)
    rho = rho';
end
if ~isrow(h)
    h = h';
end
if ~isrow(w)
    w = w';
end
d = HTLM.d;
EPW = HTLM.EPW;
PML = HTLM.PML;

%% Convert Poisson's ratio to P-wave velocity ----------------------------------------------------
% nu < 0.5 is treated as Poisson's ratio; a real Vp (>= 0.5) can also be passed directly
if nu < 0.5
    vp = vs.*sqrt(2.*(1-nu)./(1-2.*nu));
else
    vp = nu;
end

%% Finite-layer discretization -------------------------------------------------------------------
cs      =   vs(1:end-1);
cp      =   vp(1:end-1);
roS     =   rho(1:end-1);


h(isnan(h))    = [];
Lamda_min      = min(real(vs))/max(w);  % Expected minimum wavelength
v_max1          =   1*max(vs);

% Non-uniform discretization
Q               = 1:100;
factor_k        = 1/(EPW*1.5);
Z_min           = 1*Lamda_min;
d_Q             = Z_min * factor_k * exp(Q*factor_k);
h_profile       = cumsum(h);
% adjust the layering
store_loc       = [];
for zz = 1 : length(h_profile)
    h_Q         = cumsum(d_Q);
    [~, loc]    = find(h_Q > h_profile(zz));
    store_loc(zz) = loc(1);
    diffr       = h_Q(loc(1)) - h_profile(zz);
    d_Q(loc(1)) = d_Q(loc(1)) - diffr;
end
h_Q             = cumsum(d_Q);
d_Q (h_Q  > (sum(h)+1)) = [];
nDivS           = [store_loc(1) diff(store_loc)];


%% Bottom half-space [using PMDLs] ---------------------------------------------------------------
csB     =   vs(end) ;
cpB     =   vp(end);
roB     =   rho(end);
kzTol   =   10^-4;

%% Stiffness matrix
% tic
[Kzz, Kzy, Kyy, M] = stiffness1(cs, cp, roS, nDivS, h, d, d_Q);
% toc

%% Half-space

[K0, K2, M] = PMDL(Kzz,Kzy,Kyy,M,csB,cpB,roB, PML);

%% Eigen value
v_min = min(cs);
[kz, evcR, evcL] = Eigen(w, M, K0, K2, v_max1, v_min);

[kz, amp_v] = filtering (w, kz, kzTol, evcR, evcL, K2, v_max1, v_min, Measured_component);


%% Dispersion curve ------------------------------------------------------------------------------
v = disp_curve(w,kz,cs,csB,cp,cpB);

%% Extract dominant mode -------------------------------------------------------------------------
amp_v1 = abs(real(amp_v)); amp_v1(amp_v1==0)=NaN; amp_v1 = amp_v1./max(amp_v1);
z = amp_v1(:,:);
Domi_v = zeros(length(w),1);

for ii         = 1 : length(w)
    v_domi     = v(ii,(z(:,ii)==1));
    
    if isempty(v_domi)
        % disp(['v_domi is empty at index ii = ', num2str(ii)]);
        v_domi = cs(1);
    end

    Domi_v(ii) = v_domi(1);
end


end


%% Sub-function of : Main_code.m
%  Sub functions : shape_fun.m

% Input :
%       Finite layer parameters (cp, cs, roS, h)
%       nDivs - Nomber of sub layer in each layer
%       GP - Gaus quadrature point and their weight
%       d -  Order of Polynomial
% Output :
%       Stiffness matrices for finite layer


function [Kzz, Kzy, Kyy, M] = stiffness1(cs, cp, roS, nDivS, h, d, d_Q)



% syms x
% [N1, B1] = shape_fun(d);
% % [N1,B1] = shape_fn(d);
n = d+1; % n=d+1 produce same result as kausel (d=1), while n=d produce better result, but create problem while modelling plates; for plate n=d+1 should be taken;
[GP] = lgwt(n,-1,1);
% [GP1] = lgwt(n,-1,1);

% Note : Order of shape function elements (d) - linear(1); quadratic(2);
%        cubic(3); quartic(4);
%        Ordert of stiffness matrix element (order of N'N) (p) - linear(2);
%        quadratic(4); cubic(6); quartic(8);
%        Required integration points p=2n-1

pt = d + 1; % before 07-12-22 it was d + 1
nL = length(cs);
NumL = sum(nDivS);
Kzz = zeros(2*pt*NumL-2*(NumL-1));
Kzy = Kzz; Kyy = Kzz; M = Kzz;

mu = roS.*cs.^2;
Lam =  roS.*cp.^2-2*mu;

e1 = 2*pt -1;
e2 = e1-1;
Ldis =[1  e1:e2:length(Kzz)];
Ldis1 = Ldis(2:end)+1;
j = 0;

for i = 1:nL                                    % starting for each layer
    Dzz = [Lam(i)+2*mu(i) 0 ; 0 mu(i)];
    Dzy = [0 Lam(i) ;  mu(i) 0];
    Dyy = [mu(i) 0 ; 0 Lam(i)+2*mu(i)];

    for k = 1 : nDivS(i)                        % for each sub layer
        j = j+1;
        for m = 1 : size(GP,2)                  % at each intrigation point of sub layers
            p_wt = GP(:,m);
            p = p_wt(1);
            wt = p_wt(2);
            [N1,B1] = shape_fn(d,p); % before 07-12-22 it was d
            N = double(N1');
            N = kron(N,eye(2));

            Bs = double(B1');

            L = d_Q(j);
            J = L/2;
            B = Bs/J; B = kron(B,eye(2));
            rc = Ldis(j):Ldis1(j);
            %% Calculate elements
            Kzz(rc,rc) = Kzz(rc,rc)   + (N.'*Dzz*N) * (J*wt);
            Kzy(rc,rc) = Kzy(rc,rc)   - (N.'*Dzy*B) * (J*wt);
            Kzy(rc,rc) = Kzy(rc,rc)   + (B.'*Dzy.'*N)*J*wt; % This Kzy = (Kzy^T-Kzy) in the equation (7)
            Kyy(rc,rc) = Kyy(rc,rc)   + (B.'*Dyy*B) * (J*wt);
            M(rc,rc)   = M(rc,rc)     + (N.'*eye(2)*N) * (J*wt*roS(i));
        end
    end
end


end


%% PMDL

% Adds the Perfectly-Matched Discrete Layer absorbing half-space; returns K0, K2, M.
function [K0, K2, M] = PMDL(Kzz,Kzy,Kyy,M,csB,cpB,roB,nDivB)
Dzz = [roB*cpB^2 0; 0 roB*csB^2];
Dzy = [0 roB*(cpB^2-2*csB^2); roB*csB^2 0];
Dyy = [roB*csB^2 0; 0 roB*cpB^2];

GP = [0 2];
Lpml = 1*2.^(0:nDivB-1);

sz = length(Kyy)+2*nDivB;
rc = length(Kyy)-1:2:sz;
rc1 = rc + 3;
Kzz(sz,sz)=0; Kzy(sz,sz)=0; Kyy(sz,sz) = 0; M(sz,sz)=0;

for k = 1:nDivB
    s = GP(1); wt = GP(2);
    N = [(1-s) (1+s)]./2;
    N=kron(N,eye(2));
    Bs = [-1 1]./2;
    L = Lpml(k);
    L = L/2;
    B = Bs/L;
    B=kron(B,eye(2));
    pp = rc(k):rc1(k);
    Kzz(pp,pp)=Kzz(pp,pp)    + (N.'*Dzz*N)*(L*wt);
    Kzy(pp,pp)=Kzy(pp,pp)    - (N.'*Dzy*B)*(L*wt);
    Kzy(pp,pp)=Kzy(pp,pp)    + (B.'*Dzy.'*N)*L*wt;
    Kyy(pp,pp)=Kyy(pp,pp)    + (B.'*Dyy*B)*(L*wt);
    M(pp,pp)=M(pp,pp)        + (N.'*eye(2)*N)*(L*wt*roB);

end
Kzz=Kzz(1:end-2,1:end-2); Kzy=Kzy(1:end-2,1:end-2); Kyy=Kyy(1:end-2,1:end-2); M=M(1:end-2,1:end-2);
%%
Z=zeros(size(Kzz,1)/2);
z=1:2:size(Kzz,1)-1;
y=2:2:size(Kzz,1);

K2=([Kzz(z,z) Z; -Kzy(y,z) Kzz(y,y)]); % z- odd points (Horizontal); y-even points (Vertical)
K0=([Kyy(z,z) Kzy(z,y);Z Kyy(y,y)]);
M=([M(z,z) Z;Z M(y,y)]);

end





%%
% Legendre-Gauss quadrature nodes and weights on the interval [a, b].
function [GP]=lgwt(N,a,b)
N=N-1;
N1=N+1; N2=N+2;

xu=linspace(-1,1,N1)';

% Initial guess
y=cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2);

% Legendre-Gauss Vandermonde Matrix
L=zeros(N1,N2);

% Derivative of LGVM
Lp=zeros(N1,N2);

% Compute the zeros of the N+1 Legendre Polynomial
% using the recursion relation and the Newton-Raphson method

y0=2;

% Iterate until new points are uniformly within epsilon of old points
while max(abs(y-y0))>eps


    L(:,1)=1;
    Lp(:,1)=0;

    L(:,2)=y;
    Lp(:,2)=1;

    for k=2:N1
        L(:,k+1)=( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k;
    end

    Lp=(N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);

    y0=y;
    y=y0-L(:,N2)./Lp;

end

% Linear map from[-1,1] to [a,b]
x=(a*(1-y)+b*(1+y))/2;

% Compute the weights
w=(b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2;

x = flip(x'); w = flip(w');
GP = [x;w];
end

%%



%%
% Solves the per-frequency quadratic eigenvalue problem for wavenumbers (kz) and eigenvectors.
function[kz, evcR, evcL] = Eigen(w, M, K0, K2, v_max, v_min)

% Only a few propagating modes are physical (the first 5 are used downstream).
% Computing every eigenvalue with eigs is slow, so instead we shift-invert
% around the physical wavenumber band [kmin, kmax] and extract nEig modes.
% Increase nEig if dispersion branches go missing; decrease for more speed.
nEig = min(20, size(M,2));
% nEig = max(20, size(M,2)); %for the slow option


evcR = zeros(size(M,1), nEig, length(w)); evcL = evcR;
kz   = zeros(nEig, length(w));

%% Eigen value problem

K0 = sparse(K0);
M  = sparse(M);
K2 = sparse(K2);

parfor ii = 1 : length(w)
    f = 2*pi* w(ii);

    % Shift into the middle of the physical wavenumber band (in kz^2 space)
    kmin  = f / v_max;
    kmax  = f / (0.8 * v_min);
    sigma = ((kmin + kmax) / 2)^2;

    % [evcRi, kzi, evcLi] = eigs(K0 - f^2*M, -K2, size(M,2)); kzi = sqrt(diag(kzi)); % old: ALL eigenvalues (slow)
    [evcRi, kzi, evcLi]  = eigs(K0 - f^2 * M, -K2, nEig, sigma);  kzi = sqrt(diag(kzi));

    evcR(:,:,ii) = (evcRi); evcL(:,:,ii) = (evcLi);
    kz(:,ii) = kzi;
end


end
%%


%% Filtering
% Selects and normalizes the physical propagating modes; returns wavenumbers and modal amplitudes.
function [kz, amp_v] = filtering (w, kz, kzTol, evcR, evcL, K2, v_max, v_min, Measured_component)
kz1 = zeros(size(kz));
amp_v = zeros(size(kz.')); 

for i=1:length(w)

    kzi = kz(:,i);

    %% Mode sorting

    % All -ve real parts propagates in +ve direction and -ve real parts
    % propagates opposite direction, therefore only +ve real parts will be
    % considered.
    % +ve real with 0 imaginary - Real mode
    % +ve real with -ve imaginary - decaying mode
    % +ve ral with +ve imaginary - amplitude exponentially increasing

    kzi(~isfinite(kzi))  = 0;
    kzi(real(kzi)<0)     = 0;           % eliminating all -ve real part
    kzi(imag(kzi)>kzTol) = 0;           % eliminating all +ve complex part

    f1                  = w(i)*2*pi;
    kmin                = f1/v_max;
    kmax                = f1/(0.8*v_min);
    kzi(real(kzi)>kmax) = 0;
    kzi(real(kzi)<kmin) = 0;

    kzi(abs(imag(kzi)) > kzTol) = 0;       % elastic: eliminate residual complex modes
    kz1(:,i) = kzi;                         % Purely propagating mode

    %%

    [r,~] = find(kz(:,i)==kzi);
    kzi2 = sort(kzi); kzi2 = flip(kzi2); kzi2=kzi2(1:length(r));

    r2 = zeros(1, length(r));
    for pp = 1 : length(r)
        [r1, ~] = find(kzi2(pp)==kzi);
        r2(pp) = r1;
    end
    r=r2;
    ampv = 0; amph = 0;
    evcs = evcR(:,:,i);

    % Normalization according to Kausel and Park

    evcL1 = evcL(:,:,i);
    evcL2 = zeros(size(evcL1));
    evcL2(1:end/2,:) = evcs(1:end/2,:) .* transpose((kz(:,i)));
    evcL2(end/2+1:end,:) = evcs(end/2+1:end,:) .* transpose(1./kz(:,i));
    fct = (transpose(evcL2) * K2 * evcs) .* transpose(kz(:,i));
    fct = diag(fct);
    evcs = evcs .* transpose(sqrt(1 ./ fct));  % Element-wise multiplication

    % %
    for j = 1:length(r)
        evc = evcs(:,r(j));
        fctr = 1;
        evc1 = abs(evc*fctr);
        ampv(j) = evc1(end/2+1);
        amph(j) = evc1(1);
    end
    if strcmpi(Measured_component, 'RCPM')
        amp_v(i,1:length(ampv)) = ampv.*amph; %%%% This multiplication effectively giving radial component due to a vertical load on the surface
    elseif strcmpi(Measured_component, 'VCPM')
        amp_v(i,1:length(ampv)) = ampv.*ampv; %%%% This multiplication effectively giving vertical component due to a vertical load on the surface
    end

    % clear amp
end
kz = kz1;
amp_v = transpose(amp_v);
end
%%

%% Dispersion curve

% Converts wavenumbers to modal phase-velocity curves.
function [v,v1] = disp_curve(w,kz,cs,csB,cp,cpB)

kz(kz==0)=NaN;
f = repmat(w,size(kz,1),1);
f = f.*2*pi;
v = f./kz;

v(v>max([cs csB])-1)=NaN;
vmin = 0.8*min(real([cs csB]));
v(v<vmin) = NaN;
v1 = v;
v = sort(v);
v = v';
v=real(v);


end
%%



% Lagrange shape functions N and their derivatives B for polynomial order d.
function [N1,B1] = shape_fn(d,x)

% syms x

% for 1st order (d=1)
if d==1
    N1 = [(1-x)/2;(1+x)/2];
    B1 = [-1/2;1/2];
end

% for 2nd order (d=2)
if d==2
    N1 = [(x*(x - 1))/2; (1 - x^2); (x*(x/2 + 1/2))];
    B1 = [(x - 1/2); -2*x; (x+1/2)];
end

% for 3rd order (d=3)
if d==3
    N1 = [(-(3*((3*x)/2 + 1/2)*(x - 1)*(x - 1/3))/8)
        ((9*((3*x)/2 + 3/2)*(x - 1)*(x - 1/3))/8)
        (-(9*((3*x)/4 + 3/4)*(x - 1)*(x + 1/3))/4)
        ((9*(x/2 + 1/2)*(x - 1/3)*(x + 1/3))/8)];
    B1 = [(9*x)/8 - (27*x^2)/16 + 1/16
        (81*x^2)/16 - (9*x)/8 - 27/16
        27/16 - (81*x^2)/16 - (9*x)/8
        (27*x^2)/16 + (9*x)/8 - 1/16];
end

% for 4th order (d==4)
if d==4
    N1 = [(x*(2*x + 1)*(x - 1)*(x - 1/2))/3
        -(4*x*(2*x + 2)*(x - 1)*(x - 1/2))/3
        4*x^4 - 5*x^2 + 1
        -4*x*((2*x)/3 + 2/3)*(x - 1)*(x + 1/2)
        (4*x*(x/2 + 1/2)*(x - 1/2)*(x + 1/2))/3];

    B1 = [(8*x^3)/3 - 2*x^2 - x/3 + 1/6
        - (32*x^3)/3 + 4*x^2 + (16*x)/3 - 4/3
        16*x^3 - 10*x
        - (32*x^3)/3 - 4*x^2 + (16*x)/3 + 4/3
        (8*x^3)/3 + 2*x^2 - x/3 - 1/6];
end

% for 5th order (d=5)
if d==5
    N1 =[-(125*((5*x)/2 + 3/2)*(x - 1)*(x - 1/5)*(x + 1/5)*(x - 3/5))/384
        (625*((5*x)/2 + 5/2)*(x - 1)*(x - 1/5)*(x + 1/5)*(x - 3/5))/384
        -(625*((5*x)/4 + 5/4)*(x - 1)*(x - 1/5)*(x - 3/5)*(x + 3/5))/96
        (625*((5*x)/6 + 5/6)*(x - 1)*(x + 1/5)*(x - 3/5)*(x + 3/5))/64
        -(625*((5*x)/8 + 5/8)*(x - 1)*(x - 1/5)*(x + 1/5)*(x + 3/5))/96
        (625*(x/2 + 1/2)*(x - 1/5)*(x + 1/5)*(x - 3/5)*(x + 3/5))/384];

    B1 =[- (3125*x^4)/768 + (625*x^3)/192 + (125*x^2)/128 - (125*x)/192 - 3/256
        (15625*x^4)/768 - (625*x^3)/64 - (1625*x^2)/128 + (325*x)/64 + 125/768
        - (15625*x^4)/384 + (625*x^3)/96 + (2125*x^2)/64 - (425*x)/96 - 375/128
        (15625*x^4)/384 + (625*x^3)/96 - (2125*x^2)/64 - (425*x)/96 + 375/128
        - (15625*x^4)/768 - (625*x^3)/64 + (1625*x^2)/128 + (325*x)/64 - 125/768
        (3125*x^4)/768 + (625*x^3)/192 - (125*x^2)/128 - (125*x)/192 + 3/256];
end

% for 6th order (d==6)
if d==6
    N1 = [(27*x*(3*x + 2)*(x - 1)*(x - 1/3)*(x + 1/3)*(x - 2/3))/80
        -(81*x*(3*x + 3)*(x - 1)*(x - 1/3)*(x + 1/3)*(x - 2/3))/40
        (81*x*((3*x)/2 + 3/2)*(x - 1)*(x - 1/3)*(x - 2/3)*(x + 2/3))/8
        - (81*x^6)/4 + (63*x^4)/2 - (49*x^2)/4 + 1
        (81*x*((3*x)/4 + 3/4)*(x - 1)*(x + 1/3)*(x - 2/3)*(x + 2/3))/4
        -(81*x*((3*x)/5 + 3/5)*(x - 1)*(x - 1/3)*(x + 1/3)*(x + 2/3))/8
        (81*x*(x/2 + 1/2)*(x - 1/3)*(x + 1/3)*(x - 2/3)*(x + 2/3))/40];
    B1 = [(243*x^5)/40 - (81*x^4)/16 - (9*x^3)/4 + (27*x^2)/16 + x/10 - 1/20
        - (729*x^5)/20 + (81*x^4)/4 + 27*x^3 - (27*x^2)/2 - (27*x)/20 + 9/20
        (729*x^5)/8 - (405*x^4)/16 - (351*x^3)/4 + (351*x^2)/16 + (27*x)/2 - 9/4
        -(x*(243*x^4 - 252*x^2 + 49))/2
        (729*x^5)/8 + (405*x^4)/16 - (351*x^3)/4 - (351*x^2)/16 + (27*x)/2 + 9/4
        - (729*x^5)/20 - (81*x^4)/4 + 27*x^3 + (27*x^2)/2 - (27*x)/20 - 9/20
        (243*x^5)/40 + (81*x^4)/16 - (9*x^3)/4 - (27*x^2)/16 + x/10 + 1/20];
end

% for 7th order (d=7)
if d==7
    N1 =[-(16807*((7*x)/2 + 5/2)*(x - 1)*(x - 1/7)*(x + 1/7)*(x - 3/7)*(x + 3/7)*(x - 5/7))/46080
        (117649*((7*x)/2 + 7/2)*(x - 1)*(x - 1/7)*(x + 1/7)*(x - 3/7)*(x + 3/7)*(x - 5/7))/46080
        -(117649*((7*x)/4 + 7/4)*(x - 1)*(x - 1/7)*(x + 1/7)*(x - 3/7)*(x - 5/7)*(x + 5/7))/7680
        (117649*((7*x)/6 + 7/6)*(x - 1)*(x - 1/7)*(x - 3/7)*(x + 3/7)*(x - 5/7)*(x + 5/7))/3072
        -(117649*((7*x)/8 + 7/8)*(x - 1)*(x + 1/7)*(x - 3/7)*(x + 3/7)*(x - 5/7)*(x + 5/7))/2304
        (117649*((7*x)/10 + 7/10)*(x - 1)*(x - 1/7)*(x + 1/7)*(x + 3/7)*(x - 5/7)*(x + 5/7))/3072
        -(117649*((7*x)/12 + 7/12)*(x - 1)*(x - 1/7)*(x + 1/7)*(x - 3/7)*(x + 3/7)*(x + 5/7))/7680
        (117649*(x/2 + 1/2)*(x - 1/7)*(x + 1/7)*(x - 3/7)*(x + 3/7)*(x - 5/7)*(x + 5/7))/46080];

    B1 =[- (823543*x^6)/92160 + (117649*x^5)/15360 + (84035*x^4)/18432 - (16807*x^3)/4608 - (12691*x^2)/30720 + (12691*x)/46080 + 5/2048
        (5764801*x^6)/92160 - (117649*x^5)/3072 - (991613*x^4)/18432 + (141659*x^3)/4608 + (171157*x^2)/30720 - (24451*x)/9216 - 343/10240
        - (5764801*x^6)/30720 + (352947*x^5)/5120 + (420175*x^4)/2048 - (36015*x^3)/512 - (445557*x^2)/10240 + (63651*x)/5120 + 1715/6144
        (5764801*x^6)/18432 - (117649*x^5)/3072 - (6974905*x^4)/18432 + (199283*x^3)/4608 + (648613*x^2)/6144 - (92659*x)/9216 - 8575/2048
        - (5764801*x^6)/18432 - (117649*x^5)/3072 + (6974905*x^4)/18432 + (199283*x^3)/4608 - (648613*x^2)/6144 - (92659*x)/9216 + 8575/2048
        (5764801*x^6)/30720 + (352947*x^5)/5120 - (420175*x^4)/2048 - (36015*x^3)/512 + (445557*x^2)/10240 + (63651*x)/5120 - 1715/6144
        - (5764801*x^6)/92160 - (117649*x^5)/3072 + (991613*x^4)/18432 + (141659*x^3)/4608 - (171157*x^2)/30720 - (24451*x)/9216 + 343/10240
        (823543*x^6)/92160 + (117649*x^5)/15360 - (84035*x^4)/18432 - (16807*x^3)/4608 + (12691*x^2)/30720 + (12691*x)/46080 - 5/2048];
end

% for 8th order (d==8)
if d==8
    N1 = [(128*x*(4*x + 3)*(x - 1)*(x - 1/2)*(x + 1/2)*(x - 1/4)*(x + 1/4)*(x - 3/4))/315
        -(1024*x*(4*x + 4)*(x - 1)*(x - 1/2)*(x + 1/2)*(x - 1/4)*(x + 1/4)*(x - 3/4))/315
        (1024*x*(2*x + 2)*(x - 1)*(x - 1/2)*(x - 1/4)*(x + 1/4)*(x - 3/4)*(x + 3/4))/45
        -(1024*x*((4*x)/3 + 4/3)*(x - 1)*(x - 1/2)*(x + 1/2)*(x - 1/4)*(x - 3/4)*(x + 3/4))/15
        (1024*(x - 1)*(x + 1)*(x - 1/2)*(x + 1/2)*(x - 1/4)*(x + 1/4)*(x - 3/4)*(x + 3/4))/9
        -(1024*x*((4*x)/5 + 4/5)*(x - 1)*(x - 1/2)*(x + 1/2)*(x + 1/4)*(x - 3/4)*(x + 3/4))/9
        (1024*x*((2*x)/3 + 2/3)*(x - 1)*(x + 1/2)*(x - 1/4)*(x + 1/4)*(x - 3/4)*(x + 3/4))/15
        -(1024*x*((4*x)/7 + 4/7)*(x - 1)*(x - 1/2)*(x + 1/2)*(x - 1/4)*(x + 1/4)*(x + 3/4))/45
        (1024*x*(x/2 + 1/2)*(x - 1/2)*(x + 1/2)*(x - 1/4)*(x + 1/4)*(x - 3/4)*(x + 3/4))/315];
    B1 = [(4096*x^7)/315 - (512*x^6)/45 - (128*x^5)/15 + (64*x^4)/9 + (56*x^3)/45 - (14*x^2)/15 - x/35 + 1/70
        - (32768*x^7)/315 + (1024*x^6)/15 + (512*x^5)/5 - 64*x^4 - (256*x^3)/15 + (48*x^2)/5 + (128*x)/315 - 16/105
        (16384*x^7)/45 - (7168*x^6)/45 - (6656*x^5)/15 + (1664*x^4)/9 + (5408*x^3)/45 - (676*x^2)/15 - (16*x)/5 + 4/5
        - (32768*x^7)/45 + (7168*x^6)/45 + (14848*x^5)/15 - (1856*x^4)/9 - (15616*x^3)/45 + (976*x^2)/15 + (128*x)/5 - 16/5
        (2*x*(4096*x^6 - 5760*x^4 + 2184*x^2 - 205))/9
        - (32768*x^7)/45 - (7168*x^6)/45 + (14848*x^5)/15 + (1856*x^4)/9 - (15616*x^3)/45 - (976*x^2)/15 + (128*x)/5 + 16/5
        (16384*x^7)/45 + (7168*x^6)/45 - (6656*x^5)/15 - (1664*x^4)/9 + (5408*x^3)/45 + (676*x^2)/15 - (16*x)/5 - 4/5
        - (32768*x^7)/315 - (1024*x^6)/15 + (512*x^5)/5 + 64*x^4 - (256*x^3)/15 - (48*x^2)/5 + (128*x)/315 + 16/105
        (4096*x^7)/315 + (512*x^6)/45 - (128*x^5)/15 - (64*x^4)/9 + (56*x^3)/45 + (14*x^2)/15 - x/35 - 1/70];
end

% for 9th order (d==9)
if d==9
    N1 =[-(531441*((9*x)/2 + 7/2)*(x - 1)*(x - 1/3)*(x + 1/3)*(x - 1/9)*(x + 1/9)*(x - 5/9)*(x + 5/9)*(x - 7/9))/1146880
        (4782969*((9*x)/2 + 9/2)*(x - 1)*(x - 1/3)*(x + 1/3)*(x - 1/9)*(x + 1/9)*(x - 5/9)*(x + 5/9)*(x - 7/9))/1146880
        -(4782969*((9*x)/4 + 9/4)*(x - 1)*(x - 1/3)*(x + 1/3)*(x - 1/9)*(x + 1/9)*(x - 5/9)*(x - 7/9)*(x + 7/9))/143360
        (4782969*((3*x)/2 + 3/2)*(x - 1)*(x - 1/3)*(x - 1/9)*(x + 1/9)*(x - 5/9)*(x + 5/9)*(x - 7/9)*(x + 7/9))/40960
        -(4782969*((9*x)/8 + 9/8)*(x - 1)*(x - 1/3)*(x + 1/3)*(x - 1/9)*(x - 5/9)*(x + 5/9)*(x - 7/9)*(x + 7/9))/20480
        (4782969*((9*x)/10 + 9/10)*(x - 1)*(x - 1/3)*(x + 1/3)*(x + 1/9)*(x - 5/9)*(x + 5/9)*(x - 7/9)*(x + 7/9))/16384
        -(4782969*((3*x)/4 + 3/4)*(x - 1)*(x + 1/3)*(x - 1/9)*(x + 1/9)*(x - 5/9)*(x + 5/9)*(x - 7/9)*(x + 7/9))/20480
        (4782969*((9*x)/14 + 9/14)*(x - 1)*(x - 1/3)*(x + 1/3)*(x - 1/9)*(x + 1/9)*(x + 5/9)*(x - 7/9)*(x + 7/9))/40960
        -(4782969*((9*x)/16 + 9/16)*(x - 1)*(x - 1/3)*(x + 1/3)*(x - 1/9)*(x + 1/9)*(x - 5/9)*(x + 5/9)*(x + 7/9))/143360
        (4782969*(x/2 + 1/2)*(x - 1/3)*(x + 1/3)*(x - 1/9)*(x + 1/9)*(x - 5/9)*(x + 5/9)*(x - 7/9)*(x + 7/9))/1146880];
    B1 = [- (43046721*x^8)/2293760 + (4782969*x^7)/286720 + (1240029*x^6)/81920 - (531441*x^5)/40960 - (102789*x^4)/32768 + (102789*x^3)/40960 + (87183*x^2)/573440 - (29061*x)/286720 - 35/65536
        (387420489*x^8)/2293760 - (4782969*x^7)/40960 - (15411789*x^6)/81920 + (5137263*x^5)/40960 + (1449981*x^4)/32768 - (1127763*x^3)/40960 - (1288143*x^2)/573440 + (47709*x)/40960 + 3645/458752
        - (387420489*x^8)/573440 + (4782969*x^7)/14336 + (3720087*x^6)/4096 - (885735*x^5)/2048 - (2473497*x^4)/8192 + (274833*x^3)/2048 + (496449*x^2)/28672 - (91935*x)/14336 - 5103/81920
        (129140163*x^8)/81920 - (4782969*x^7)/10240 - (48361131*x^6)/20480 + (6908733*x^5)/10240 + (8063469*x^4)/8192 - (2687823*x^3)/10240 - (2155491*x^2)/20480 + (239499*x)/10240 + 6615/16384
        - (387420489*x^8)/163840 + (4782969*x^7)/20480 + (152523567*x^6)/40960 - (7263027*x^5)/20480 - (28258227*x^4)/16384 + (3139803*x^3)/20480 + (9974907*x^2)/40960 - (369441*x)/20480 - 178605/32768
        (387420489*x^8)/163840 + (4782969*x^7)/20480 - (152523567*x^6)/40960 - (7263027*x^5)/20480 + (28258227*x^4)/16384 + (3139803*x^3)/20480 - (9974907*x^2)/40960 - (369441*x)/20480 + 178605/32768
        - (129140163*x^8)/81920 - (4782969*x^7)/10240 + (48361131*x^6)/20480 + (6908733*x^5)/10240 - (8063469*x^4)/8192 - (2687823*x^3)/10240 + (2155491*x^2)/20480 + (239499*x)/10240 - 6615/16384
        (387420489*x^8)/573440 + (4782969*x^7)/14336 - (3720087*x^6)/4096 - (885735*x^5)/2048 + (2473497*x^4)/8192 + (274833*x^3)/2048 - (496449*x^2)/28672 - (91935*x)/14336 + 5103/81920
        - (387420489*x^8)/2293760 - (4782969*x^7)/40960 + (15411789*x^6)/81920 + (5137263*x^5)/40960 - (1449981*x^4)/32768 - (1127763*x^3)/40960 + (1288143*x^2)/573440 + (47709*x)/40960 - 3645/458752
        (43046721*x^8)/2293760 + (4782969*x^7)/286720 - (1240029*x^6)/81920 - (531441*x^5)/40960 + (102789*x^4)/32768 + (102789*x^3)/40960 - (87183*x^2)/573440 - (29061*x)/286720 + 35/65536];

end

% for 10 th order (d==10)
if d==10
    N1 =[(78125*x*(5*x + 4)*(x - 1)*(x - 1/5)*(x + 1/5)*(x - 2/5)*(x + 2/5)*(x - 3/5)*(x + 3/5)*(x - 4/5))/145152
        -(390625*x*(5*x + 5)*(x - 1)*(x - 1/5)*(x + 1/5)*(x - 2/5)*(x + 2/5)*(x - 3/5)*(x + 3/5)*(x - 4/5))/72576
        (390625*x*((5*x)/2 + 5/2)*(x - 1)*(x - 1/5)*(x + 1/5)*(x - 2/5)*(x + 2/5)*(x - 3/5)*(x - 4/5)*(x + 4/5))/8064
        -(390625*x*((5*x)/3 + 5/3)*(x - 1)*(x - 1/5)*(x + 1/5)*(x - 2/5)*(x - 3/5)*(x + 3/5)*(x - 4/5)*(x + 4/5))/2016
        (390625*x*((5*x)/4 + 5/4)*(x - 1)*(x - 1/5)*(x - 2/5)*(x + 2/5)*(x - 3/5)*(x + 3/5)*(x - 4/5)*(x + 4/5))/864
        -(390625*(x - 1)*(x + 1)*(x - 1/5)*(x + 1/5)*(x - 2/5)*(x + 2/5)*(x - 3/5)*(x + 3/5)*(x - 4/5)*(x + 4/5))/576
        (390625*x*((5*x)/6 + 5/6)*(x - 1)*(x + 1/5)*(x - 2/5)*(x + 2/5)*(x - 3/5)*(x + 3/5)*(x - 4/5)*(x + 4/5))/576
        -(390625*x*((5*x)/7 + 5/7)*(x - 1)*(x - 1/5)*(x + 1/5)*(x + 2/5)*(x - 3/5)*(x + 3/5)*(x - 4/5)*(x + 4/5))/864
        (390625*x*((5*x)/8 + 5/8)*(x - 1)*(x - 1/5)*(x + 1/5)*(x - 2/5)*(x + 2/5)*(x + 3/5)*(x - 4/5)*(x + 4/5))/2016
        -(390625*x*((5*x)/9 + 5/9)*(x - 1)*(x - 1/5)*(x + 1/5)*(x - 2/5)*(x + 2/5)*(x - 3/5)*(x + 3/5)*(x + 4/5))/8064
        (390625*x*(x/2 + 1/2)*(x - 1/5)*(x + 1/5)*(x - 2/5)*(x + 2/5)*(x - 3/5)*(x + 3/5)*(x - 4/5)*(x + 4/5))/72576];
    B1 = [(1953125*x^9)/72576 - (390625*x^8)/16128 - (78125*x^7)/3024 + (78125*x^6)/3456 + (8125*x^5)/1152 - (40625*x^4)/6912 - (5125*x^3)/9072 + (5125*x^2)/12096 + x/126 - 1/252
        - (9765625*x^9)/36288 + (390625*x^8)/2016 + (1015625*x^7)/3024 - (203125*x^6)/864 - (59375*x^5)/576 + (59375*x^4)/864 + (157625*x^3)/18144 - (31525*x^2)/6048 - (125*x)/1008 + 25/504
        (9765625*x^9)/8064 - (1171875*x^8)/1792 - (1796875*x^7)/1008 + (359375*x^6)/384 + (90625*x^5)/128 - (90625*x^4)/256 - (67625*x^3)/1008 + (13525*x^2)/448 + (125*x)/126 - 25/84
        - (9765625*x^9)/3024 + (390625*x^8)/336 + (1328125*x^7)/252 - (265625*x^6)/144 - (40625*x^5)/16 + (40625*x^4)/48 + (546125*x^3)/1512 - (109225*x^2)/1008 - (125*x)/21 + 25/21
        (9765625*x^9)/1728 - (390625*x^8)/384 - (78125*x^7)/8 + (109375*x^6)/64 + (1009375*x^5)/192 - (1009375*x^4)/1152 - (208625*x^3)/216 + (41725*x^2)/288 + (125*x)/3 - 25/6
        -(x*(1953125*x^8 - 3437500*x^6 + 1918125*x^4 - 382250*x^2 + 21076))/288
        (9765625*x^9)/1728 + (390625*x^8)/384 - (78125*x^7)/8 - (109375*x^6)/64 + (1009375*x^5)/192 + (1009375*x^4)/1152 - (208625*x^3)/216 - (41725*x^2)/288 + (125*x)/3 + 25/6
        - (9765625*x^9)/3024 - (390625*x^8)/336 + (1328125*x^7)/252 + (265625*x^6)/144 - (40625*x^5)/16 - (40625*x^4)/48 + (546125*x^3)/1512 + (109225*x^2)/1008 - (125*x)/21 - 25/21
        (9765625*x^9)/8064 + (1171875*x^8)/1792 - (1796875*x^7)/1008 - (359375*x^6)/384 + (90625*x^5)/128 + (90625*x^4)/256 - (67625*x^3)/1008 - (13525*x^2)/448 + (125*x)/126 + 25/84
        - (9765625*x^9)/36288 - (390625*x^8)/2016 + (1015625*x^7)/3024 + (203125*x^6)/864 - (59375*x^5)/576 - (59375*x^4)/864 + (157625*x^3)/18144 + (31525*x^2)/6048 - (125*x)/1008 - 25/504
        (1953125*x^9)/72576 + (390625*x^8)/16128 - (78125*x^7)/3024 - (78125*x^6)/3456 + (8125*x^5)/1152 + (40625*x^4)/6912 - (5125*x^3)/9072 - (5125*x^2)/12096 + x/126 + 1/252];
end

if d > 10
    disp('shape function exceeds; modify the code')
    return
end

end

%%
