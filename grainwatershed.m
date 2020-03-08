function [water] = grainwatershed(gambar,minthres)

dist2=-bwdist(~gambar); %menghitung jarak antar nilai sel (biner)
%dist2=-dist2;
dist2(~gambar)=-Inf; %yang 0 jadi infinit

%minimum threshold
dist2=imhmin(dist2,minthres);
water=watershed(dist2);
%QC pore
water(~gambar)=0;

end

