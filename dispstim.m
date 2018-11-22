cd c:/mscripts/finalstims ; 
ls 

rands = load('randarrs.mat') ; 
rands = rands.randarrs ;
[x,y] = meshgrid(-100:99,-100:99) ;
circ = sqrt(x.^2 + y.^2) < 90 ;
figure,
icount =1 ; 
for i=1:4:13
    subplot(1,4,icount) ;
    a = (squeeze(rands(1:200,1:200,i)).*circ) ; 
    a(circ==0) = 127.5 ; 
    imshow(uint8(a)) ;
    icount = icount + 1 ; 
end
    
plaid = load('plaid.mat') ; 
plaid = plaid.plaid ; 
p1 = squeeze(plaid(1:200,1:200,1)) ;
p2 = rot90(squeeze(plaid(1:200,1:200,1))) ; 
figure,
subplot(1,4,2) ;
a = ((p1/2 + p2/2).*circ) ; a(circ==0) = 127.5 ; 
imshow(uint8(a)) ; 
subplot(1,4,1) ; 
a = p1.*circ ; a(circ==0) = 127.5 ; 
imshow(uint8(a)) ; 

figure,
c1 = load('c1.mat') ; c1 = c1.c1 ; c1 = (c1(1:200,1:200).*circ) ; c1(circ==0) = 127 ; 
subplot(1,4,3) ; imshow(uint8(c1)) ; 
c2 = load('c2.mat') ; c2 = c2.c2 ; c2 = (c2(1:200,1:200).*circ) ; c2(circ==0) = 127 ; 
subplot(1,4,2) ; imshow(uint8(c2)) ;
c3 = load('c3.mat') ; c3 = c3.c3 ; c3 = (c3(1:200,1:200).*circ) ; c3(circ==0) = 127 ; 
subplot(1,4,1) ; imshow(uint8(c3)) ;



