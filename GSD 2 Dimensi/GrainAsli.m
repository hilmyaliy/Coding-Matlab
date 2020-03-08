function [gambar_asli,ukuran_gambar,x,y,r,luas] = GrainAsli(grain,rmax,rmin,xmax,ymax,grain_persen)
%UNTITLED14 Summary of this function goes here
%   Detailed explanation goes here
axis([1 xmax 1 ymax]), axis square, box on,
set(gca,'units','pixels','position',[50,50,300,300],... %[kiri bawah lebar tinggi]
    'XTick',[],'YTick',[])

%BUAT GRAIN

for i=1:grain
    
    %ukuran random
    r(i)=rmin+(rmax-rmin)*rand(1); %ukuran 3-4
    phi=360*rand();
    x(i)=xmax*rand(1); %posisi x
    y(i)=ymax*rand(1); %posisi y
    luas(i)=pi*r(i)^2;
    a(i)=rectangle('Position',[x(i),y(i),2*r(i),2*r(i)],...
              'Curvature',[1,1],...%kelengkungan
              'FaceColor','k'); %hitam
          
    %buat ruang
    axis([1 xmax 1 ymax]), axis square, box on
    set(gca,'units','pixels','XTick',[],'YTick',[]) %tick sebagai panjang sumbu
    H=gca;
    
    img2 = getframe(gca);
    imwrite(img2.cdata,'gambargrain.bmp','bmp');
    img=imread('gambargrain.bmp');
    gambar1(:,:,i)=im2bw(img); %ambil data gambar (cdata), didapat gambar biner
    gambar(:,:,i)=~gambar1(:,:,i);%grain-1-putih/pore-0-hitam
    
    delete(a(i));
    %delete(a);
end

axis off

ruang=zeros(length(gambar(:,:,1)));

nabrak=0;
ntotal=0;
ngrain=grain;


for k=1:ngrain  %grain
    %dimuali dari 2 karena ada bingkai yang nilainya 1
    for i=2:length(gambar(:,1,1)) %baris
        for j=2:length(gambar(1,:,1)) %kolom
            
            if gambar(i,j,k)==1 && ruang(i,j)==1  %
                nabrak=nabrak+1;
            
            end
            
            if gambar(i,j,k)==1
            ntotal=ntotal+1;
            end

        end
    end
    
    
    if (nabrak/ntotal*100) >= grain_persen
        gambar(:,:,k)=0;
        
    else
        ruang(:,:)=ruang(:,:)+gambar(:,:,k);
        
    end
    
    nabrak=0;
    ntotal=0;
    
end

%menghapus bingkai bagian atas kiri
ruang(1:length(ruang(:,1)),1)=0;
ruang(1,1:length(ruang(:,1)))=0;
gambar_asli=logical(ruang);


%DONE !

    
[ukuran_gambar(1), ukuran_gambar(2), ukuran_gambar(3)]=size(gambar_asli);

end

