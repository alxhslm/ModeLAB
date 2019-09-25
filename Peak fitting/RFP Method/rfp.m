function [modal_par,alpha]=rfp(rec,omega,N,wt)
 
%RFP Modal parameter estimation from frequency response function using 
% rational fraction polynomial method.
%
% Syntax: [alpha,modal_par]=rfp(rec,omega,N)
%
% rec   = FRF measurement (receptance)
% omega = frequency range vector (rad/sec).
% N     = number of degrees of freedom.
% alpha = FRF generated (receptance).
% modal_par = Modal Parameters [freq,damp,Ci,Oi]: 
%             freq = Natural frequencies (rad/sec)
%             damp = Damping ratio
%             Ci   = Amplitude modal constant
%             Oi   = Phase modal constant (degrees)
%
% Reference: Mark H.Richardson & David L.Formenti "Parameter Estimation 
%           from Frequency Response Measurements Using Rational Fraction 
%           Polynomials", 1ºIMAC Conference, Orlando, FL. November, 1982.
%**********************************************************************
%Chile, March 2002, Cristian Andrés Gutiérrez Acuña, crguti@icqmail.com
%**********************************************************************

omega = omega(:);%omega is now a column
if size(rec,2) == length(omega)
	rec=rec.';     %rec is now a column
end

if nargin < 4
    wt = 0*omega+1;
end

nom_omega=max(omega);
omega=omega./nom_omega; %omega normalization

n=2*N;   % n is number of polynomial terms in denominator
m=2*N-1; % m is number of polynomial terms in numerator

%orthogonal function that calculates the orthogonal polynomials
[phi,coeff_A] = orthogonal(omega,ones(size(omega)),m);
[theta,coeff_B] = orthogonal(omega,(abs(rec)).^2,n);

P = phi.*wt;
T = repmat(rec,1,size(theta,2)-1).*(wt.*theta(:,1:end-1));
W = rec.*(wt.*theta(:,end));
X = -2*real(P'*T);
H = 2*real(P'*W);

I = eye(size(X,2));
D = -(I-X.'*X)\(X.'*H);
C = H - X*D;   %{C} orthogonal numerator    polynomial coefficients
D = [D;1];     %{D} orthogonal denominator  polynomial coefficients

%calculation of FRF (alpha)
numer = phi*C;
denom = theta*D;
alpha = numer./denom;

%convert orthogonal coefficients -> actual polynomial coefficients
A=coeff_A*C;
A=A(end:-1:1); %{A} standard numerator polynomial coefficients

B=coeff_B*D;
B=B(end:-1:1); %{B} standard denominator polynomial coefficients

%calculation of the poles and residues
[R,P]=residue(A,B);
residuals = R(1:2:end);
poles = P(1:2:end);

%scale appropriately
residuals = residuals(end:-1:1)*nom_omega; %residues
poles = poles(end:-1:1)*nom_omega;       %poles
freq = abs(poles);                 %Natural frequencies (rad/sec)
damp = -real(poles)./abs(poles);   %Damping ratios

%compute modal participation factors
Ai=-2*(real(residuals).*real(poles)+imag(residuals).*imag(poles));
Bi=2*real(residuals);
Ar=complex(Ai,abs(poles).*Bi);
    
modal_par.H = alpha;
modal_par.omega = freq;
modal_par.zeta = damp;
modal_par.A = Ar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%