% Program No:6
% Write an M-function for performing local histogram equalization

function g = localhist(x,w,k)
f=x;
f=im2double(f);
%w=input('\nEnter the Neighborhood or Window size : ');
%k=input('\nEnter the value of the constant k (value should be between 0 and 1) : ');
%w = 31 ;
%wx = size(x,1) ; 
%wy = size(x,2) ; 
%k = 1 ; 

M=mean2(f);
z=colfilt(f,[w 1],'sliding',@std);
m=colfilt(f,[w 1],'sliding',@mean);
A=k*M./z;
g=A.*(f-m)+m;
%imshow(f), figure, 
%imagesc(g);
end