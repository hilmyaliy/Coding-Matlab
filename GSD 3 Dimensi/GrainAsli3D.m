function [the_bw,ukuran_gambar,volume_model,radius,luaspermukaan_bw,xx,yy,zz,partikel,bw_lama] = GrainAsli3D(grain,rmax,rmin,xmax,ymax,zmax,grain_persen)
%radius yg ri atas adalah radius model
%buat ruang final
the_bw=zeros(xmax,ymax,zmax);


%random 
igrain=1;
ulang=0;
while igrain<=grain

acak=normrnd(0.5,0.25);
rr=rmin+(rmax-rmin)*acak; %ukuran 3-4
rr=round(rr);

if rr>=rmin && rr<=rmax
    
radius(igrain)=rr;
    
volume(igrain)=(4/3)*pi*radius(igrain)^3;
[x,y,z] = meshgrid(-radius(igrain):radius(igrain));

bw_lama = sqrt((x).^2 + (y).^2 +(z).^2) <=radius(igrain);

xmove(igrain)=round((xmax-30)*rand(1)-radius(igrain)); %posisi x
ymove(igrain)=round((ymax-30)*rand(1)-radius(igrain)); %posisi y
zmove(igrain)=round((zmax-30)*rand(1)-radius(igrain)); %posisi z

    if (xmove(igrain)+radius(igrain))>xmax
        xmove(igrain)=xmove(igrain)-radius(igrain);
    elseif (xmove(igrain)-radius(igrain))< 0
        xmove(igrain)=xmove(igrain)+radius(igrain);
    end
    
    if (ymove(igrain)+radius(igrain))>ymax
        ymove(igrain)=ymove(igrain)-radius(igrain);
    elseif (ymove(igrain)-radius(igrain))<0
        ymove(igrain)=ymove(igrain)+radius(igrain);
    end
    
    if (zmove(igrain)+radius(igrain))>zmax
        zmove(igrain)=zmove(igrain)-radius(igrain);
    elseif (zmove(igrain)-radius(igrain))<0
        zmove(igrain)=zmove(igrain)+radius(igrain);
    end
    
    
    
%put smaller matrix into larger one
%hanya ruang yang diperbesar
bw_lama2=zeros(xmax,ymax,zmax);
for i=1:length(x)
    for j=1:length(y)
        for k=1:length(z)
            
            if bw_lama(i,j,k)==1
                bw_lama2(i,j,k)=1;
            end
            
        end
    end
end

bw_baru=zeros(xmax,ymax,zmax);

%perpindahan ke ruang matriks baru
x1=length(bw_lama2);

for i =1:length(x)
    for j=1:length(y)
        for k=1:length(z)
            
            if bw_lama2(i,j,k)==1
                %jika partikel dipindah ternyata masih di dalam frame
                if i+xmove(igrain)<=x1 && j+ymove(igrain)<=x1
                    if k+zmove(igrain)<=x1
                        
                    bw_baru(i+xmove(igrain),j+ymove(igrain),k+zmove(igrain))=bw_lama2(i,j,k);
                
                    end
                end 
            end
            
        end
    end
end
bw_baru=logical(bw_baru);
volume_model(igrain)=sum(sum(sum(bw_baru)));
daerah_bw               = regionprops3(bw_baru,'SurfaceArea');
luaspermukaan_bw(igrain)= round(daerah_bw.SurfaceArea);


%alami overlap
nabrak=0;
ntotal=0;
partikel1=zeros(xmax,ymax,zmax);
partikel1=bw_baru;

for i =1:1:xmax
    for j=1:1:ymax
        for k=1:1:zmax
            
            %terjadi kontak
            if bw_baru(i,j,k) ==1
                if bw_baru(i,j,k)==the_bw(i,j,k)
                 nabrak=nabrak+1;
                 partikel1(i,j,k)=0;
                end
                ntotal=ntotal+1;
            end
            
        end
    end
end
     

if (nabrak/ntotal*100) >= grain_persen
    ulang=ulang+1;
    if ulang==50
        grain=grain-1;
        ulang=0;
    end
    
else %jika di bawah tingkat overlap yg diminta
   
volumeHP(igrain)=ntotal-nabrak;
partikel(:,:,:,igrain)=logical(partikel1);

%plot !

%untuk membuat ruang
[xx,yy,zz] = meshgrid(1:xmax);

%figure(1), isosurface(xx,yy,zz,bw_lama2), axis equal, title('BW')
%xlabel x, ylabel y, zlabel z
%xlim([-50 200]), ylim([-50 200]), zlim([-50 200])


figure(2), isosurface(xx,yy,zz,partikel(:,:,:,igrain)), axis equal, title('Processing...')
xlabel x, ylabel y, zlabel z
xlim([-5 xmax+5]), ylim([-5 ymax+5]), zlim([-5 zmax+5])


hold on

the_bw= the_bw | partikel(:,:,:,igrain);


igrain=igrain+1;
ulang=0;

end


end


end

[ukuran_gambar(1), ukuran_gambar(2), ukuran_gambar(3)]=size(the_bw);
end

