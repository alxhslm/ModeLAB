function [P,coeff]=orthogonal(omega,q,kmax)

%ORTHOGONAL Orthogonal polynomials required for rational fraction 
% polynomials method. (This code was written to be used with rfp.m)
%
% Syntax: [P,coeff]=orthogonal(rec,omega,phitheta,kmax)
%
% rec      = FRF measurement (receptance).
% omega    = frequency range vector (rad/sec).
% phitheta = weighting function (must be 1 for phi matrix or 2 for 
%            theta matrix).  
% kmax     = degree of the polynomial.
% P        = matrix of the orthogonal polynomials evaluated at the 
%            frequencies.
% coeff    = matrix used to transform between the orthogonal polynomial 
%            coefficients and the standard polynomial.
%
% Reference: Mark H.Richardson & David L.Formenti "Parameter Estimation 
%           from Frequency Response Measurements Using Rational Fraction 
%           Polynomials", 1ºIMAC Conference, Orlando, FL. November, 1982.
%**********************************************************************
%Chile, March 2002, Cristian Andrés Gutiérrez Acuña, crguti@icqmail.com
%**********************************************************************

%orthogonal polynomials matrix
R_m1=zeros(size(omega));
R_0=1/sqrt(2*sum(q)).*ones(size(omega));
R=[R_m1,R_0];    %polynomials -1 and 0.

%transform matrix (this bit is not in the paper)
coeff=zeros(kmax+1,kmax+2);
coeff(1,2)=1/sqrt(2*sum(q));

% generate the orthogonal polynomials matrix and transform matrix 
% using the simplifited Forsythe method
for k=1:kmax
	Vkm1 = 2*sum(omega.*R(:,end).*R(:,end-1).*q);
	Sk   = omega.*R(:,end)-Vkm1*R(:,end-1);
	Dk   = sqrt(2*sum((Sk.^2).*q));
	R(:,end+1) = (Sk/Dk); %#ok
	coeff(:,k+2) = -Vkm1*coeff(:,k);
	coeff(2:k+1,k+2) = coeff(2:k+1,k+2)+coeff(1:k,k+1);
    coeff(:,k+2) = coeff(:,k+2)/Dk;
end

R=R(:,2:kmax+2);         %orthogonal polynomials matrix
coeff=coeff(:,2:kmax+2); %transform matrix

%make complex by multiplying by i^k
P = 0*R;
jk = zeros(1,kmax+1);
for k=0:kmax
   P(:,k+1)=R(:,k+1)*1i^k; %complex orthogonal polynomials matrix
   jk(1,k+1)=1i^k;
end
coeff=(jk'*jk).*coeff;    %complex transform matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%