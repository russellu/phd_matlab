function f=image_registr_MI(x)

% This is a subroutine of opti_MI_scaling.m.
%
% Originally written by K.Artyushkova
% 10_2003
% Kateryna Artyushkova
% Postdoctoral Scientist
% Department of Chemical and Nuclear Engineering
% The University of New Mexico
% (505) 277-0750
% kartyush@unm.edu 

% Modified by Hosang Jin
% 2_2005
% 
% Graduate Student
% Department of Nuclear & Radiological Engineering
% The University of Florida
% hsjin@ufl.edu

load image


img1 = imread('C:\shared\utef\alex\rgbs.png') ; 
IM1=double(rgb2gray(img1));
img2 = imread('C:\shared\utef\biz\rgbs.png') ; 
IM2=double(rgb2gray(img2));


%IM1=double(IM1);
%IM2=double(IM2);
%IM2=imresize(IM2, x(4), 'bilinear');
J=imrotate(double(IM2), x(3),'bilinear'); %rotated cropped IMAGE2
J=abs(J)*255/max(max(J));

[n1 n2]=size(IM1);
[n3 n4]=size(J);

if false%n1>n3-x(1)/2
    %{
    f=1000;
    message=strvcat('The scaling factor is too small.', 'Press Ctrl+C to stop.',...
        'Increase x0(4) and restart.');
    disp('Press Ctrl+C to stop.')
    Errordlg(message)
    pause;
    %}
else
    if x(1)>n3-n1
        x(1)=n3-n1-1;
        IM1(1:n1, 1:n2)=255;
    end
    
    if x(2)>n4-n2
        x(2)=n4-n2-1;
        IM1(1:n1, 1:n2)=255;
    end
    
    if x(1)<0
        x(1)=0;
        IM1(1:n1, 1:n2)=255;
    end
    
    if x(2)<0
        x(2)=0;
        IM1(1:n1, 1:n2)=255;
    end
    
    xt=1:n1;
    yt=1:n2;
    
    xx=round(xt+x(1));
    yy=round(yt+x(2));
    
    IM2=round(J(xx, yy)); % selecting part of IMAGE2 matching the size of IMAHE1
    
    rows=size(IM1,1);
    cols=size(IM2,2);
    N=256;
    
    h=zeros(N,N);
    
    for ii=1:rows;    %  col 
        for jj=1:cols;   %   rows
            h(IM1(ii,jj)+1,IM2(ii,jj)+1)= h(IM1(ii,jj)+1,IM2(ii,jj)+1)+1;
        end
    end
    
    [r,c] = size(h);
    b= h./(r*c); % normalized joint histogram
    y_marg=sum(b); %sum of the rows of normalized joint histogram
    x_marg=sum(b');%sum of columns of normalized joint histogran
    
    Hy=0;
    for i=1:c;    %  col
        if( y_marg(i)==0 )
            %do nothing
        else
            Hy = Hy + -(y_marg(i)*(log2(y_marg(i)))); %marginal entropy for image 1
        end
    end
    
    Hx=0;
    for i=1:r;    %rows
        if( x_marg(i)==0 )
            %do nothing
        else
            Hx = Hx + -(x_marg(i)*(log2(x_marg(i)))); %marginal entropy for image 2
        end   
    end
    h_xy = -sum(sum(b.*(log2(b+(b==0))))); % joint entropy
    
    f=-(Hx+Hy-h_xy);% Mutual information
    %x
end