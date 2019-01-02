function [modal_par,alpha]=grfp(rec,omega,N,wt)
 
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

Ndof = size(rec,2);
nom_omega=max(omega);
omega=omega./nom_omega; %omega normalization

m=2*N-1; %number of polynomial terms in numerator
n=2*N;   %number of polynomial terms in denominator

%orthogonal function that calculates the orthogonal polynomials
[phi,coeff_A] = orthogonal(omega,ones(size(omega)),m);
[theta,coeff_B] = orthogonal(omega,(mean(abs(rec),2)).^2,n);

P = kron(eye(Ndof),wt.*phi);
T = repmat(rec(:),1,size(theta,2)-1).*repmat(wt.*theta(:,1:end-1),Ndof,1);
W = rec(:).*repmat(wt.*theta(:,end),Ndof,1);
X = -2*real(P'*T);
H = 2*real(P'*W);

I = eye(size(X,2));
D = -(I-X.'*X)\(X.'*H);
C = H - X*D;   %{C} orthogonal numerator    polynomial coefficients
C = reshape(C,[],Ndof);
D = [D;1];     %{D} orthogonal denominator  polynomial coefficients

%calculation of FRF (alpha)
numer = phi*C;
denom = theta*D;
alpha = numer./denom;

%convert orthogonal coefficients -> actual polynomial coefficients
A=coeff_A*C;
A=A(end:-1:1,:); %{A} standard numerator polynomial coefficients

B=coeff_B*D;
B=B(end:-1:1); %{B} standard denominator polynomial coefficients

%calculation of the poles and residues
P = roots(B);
poles = P(1:2:end);
for i = 1:Ndof
    R=residue(A(:,i),B);
    residuals(:,i) = R(1:2:end);
end

%scale appropriately
residuals = residuals(end:-1:1,:)*nom_omega; %residues
poles = poles(end:-1:1)*nom_omega;       %poles
freq = abs(poles);                 %Natural frequencies (rad/sec)
damp = -real(poles)./abs(poles);   %Damping ratios

%compute modal participation factors
Ai=-2*(real(residuals).*real(poles)+imag(residuals).*imag(poles));
Bi=2*real(residuals);
Ar=complex(Ai,abs(poles).*Bi).';
    
modal_par.alpha = alpha;

modal_par.Frequency = freq;
modal_par.Damping = damp;

modal_par.Real = real(Ar);
modal_par.Imaginary = imag(Ar);
modal_par.Magnitude = abs(Ar);
modal_par.Phase = angle(Ar);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%