function [gambar_asli,ukuran_gambar,r,luas_input,luas,partikel2,i1,d] = GrainAsli2(grain,rmax,rmin,xmax,ymax,grain_persen,dimensi)

%UNTITLED14 Summary of this function goes here
%   Detailed explanation goes here


axis([0 xmax 0 ymax]), axis square, box on,
set(gca,'units','pixels','outerposition',[20,20,xmax+10,ymax+10],'position',[30,30,xmax,ymax],... %[kiri bawah lebar tinggi]
    'XTick',[],'YTick',[])

%BUAT GRAIN

tgrain=grain;

i1=1;
ulang1=0;
ulang2=0;
izin=0;
ruang=zeros(xmax,ymax);

ntotal=0;
nabrak=0;
perlu=1;
partikel2=zeros(xmax,ymax);

while i1<=grain
    
    %acak normal
    acak=normrnd(0.5,0.25);
    rr=rmin+(rmax-rmin)*acak;
    
    if rr>=rmin && rr<=rmax
    r(i1)=rr;
    luas_input(i1)=round(pi*r(i1)^2);
    d(i1)=2*r(i1); %diameter
    phi=360*rand();
    
    x(i1)=round((xmax)*rand(1)); %posisi centroid x
    y(i1)=round((ymax)*rand(1)); %posisi centroid y
    
    if (x(i1)+r(i1)>xmax)
        x(i1)=x(i1)-r(i1);
    elseif (x(i1)-r(i1)<0)
        x(i1)=x(i1)+r(i1);
    end
    
    if (y(i1)+r(i1)>ymax)
        y(i1)=y(i1)-r(i1);
    elseif (y(i1)-r(i1)<0)
        y(i1)=y(i1)+r(i1);
    end

    %PERSEKITARAN
    coba=zeros(xmax,ymax);
    coba(round(x(i1)),round(y(i1)))=1; %centroid
    %buat persekitaran
    sekitarx=(round(x(i1)-1)):1:(round(x(i1)+1));
    sekitary=(round(y(i1)-1)):1:(round(y(i1)+1));
    coba(sekitarx,sekitary)=ones(3);
    
    kena=0;
    
    for ii=sekitarx(1):1:sekitarx(length(sekitarx))
        for jj=sekitary(1):1:sekitary(length(sekitary))
             if coba(ii,jj) == ruang(ii,jj) %jika centroid terkena padatan lain
                kena=kena+1;
            end
        end
    end
    
    if kena==0
       izin=1; 
    end
    
    if izin==1
    %luas(i1)=pi*r(i1)^2;
    a(i1)=rectangle('Position',[x(i1)-r(i1),y(i1)-r(i1),d(i1),d(i1)],...
              'Curvature',[1,1],...%kelengkungan
              'FaceColor','k'); %hitam
          
    %buat ruang
    axis([0 xmax 0 ymax]), axis square, box on
    set(gca,'units','pixels') %tick sebagai panjang sumbu
    
    H=gca;
    
    img2 = getframe(gca);
    %imwrite(img2.cdata,'gambargrain.bmp','bmp');
    %img=imread('gambargrain.bmp');
    img=img2.cdata;
    %perkecil gambar
    img=imresize(img,0.5);
    
    gambar1(:,:,i1)=im2bw(img); %ambil data gambar (cdata), didapat gambar biner
    gambar(:,:,i1)=~gambar1(:,:,i1);%grain-1-putih/pore-0-hitam
    
    gambar_p=gambar(:,:,i1);
    
    partikel1=zeros(xmax,ymax);
    
    partikel1=double(gambar_p);
    
    %menghapus frame pojok kiri atas bernilai 1
    %partikel1(1:length(partikel1(:,1)),1)=0;
    %partikel1(1,1:length(partikel1(:,1)))=0;
    
    luas1=sum(sum(partikel1)); %luas 1 grain
    
    delete(a(i1));
    
    %---------------------------------------------------
    %izin proses ?
    
    for i=1:length(gambar(:,1,1)) %baris
        for j=1:length(gambar(1,:,1)) %kolom
            
            if gambar(i,j,i1)==1 && ruang(i,j)==1  %terjadi kontak
                nabrak=nabrak+1;
                partikel1(i,j)=0;
                
            end
            
            if gambar(i,j,i1)==1
            ntotal=ntotal+1;
            end

        end
    end
    
    
    if (nabrak/ntotal*100) >= grain_persen %jika di atas toleransi overlap
        
        %gambar(:,:,k)=0;
        ulang2=ulang2+1;
        
        if ulang2==50 %jika 10 kali tetap overlap, maka grain dihilangkan
            
            gambar(:,:,i1)=0;
            
            grain=grain-1;
            ulang2=0;
            
            if perlu==2
                perlu=3;
            elseif perlu==3
                perlu=1;
            end
            
        end
        
    else %di bawah toleransi overlap
       
            ruang(:,:)=ruang(:,:)+gambar(:,:,i1);
            luas(i1)=ntotal-nabrak;
       
            %partikel1(sekitarx,sekitary)=ones(3); %untuk centroid
            partikel2(:,:,i1)=logical(partikel1);
            i1=i1+1;
            perlu=2;
            ulang2=0;
    end
    
    nabrak=0;
    ntotal=0;
    izin=0;
    
    end
    %----------------------------------------------
   
    end    
end

if perlu==2 || perlu==3
    i1=i1-1;
    %if perlu==3
       r=r(1:i1);
       d=d(1:i1);
    %end
    
 %mengurangi kelebihan butiran
end

%axis off



%nabrak=0;
%ntotal=0;
ngrain=grain;



%menghapus bingkai bagian atas kiri (karena mengandung nilai 1)
%ruang(1:length(ruang(:,1)),1)=0;
%ruang(1,1:length(ruang(:,1)))=0;
gambar_asli=logical(ruang);


%DONE !

    
[ukuran_gambar(1), ukuran_gambar(2), ukuran_gambar(3)]=size(gambar_asli);



%if tgrain>grain
%    i2=grain;
%    
%    while i2<tgrain
%     acak=normrnd(0.5,0.25);
%    rr=rmin+(rmax-rmin)*acak;
    
%        if rr>=rmin && rr<=rmax
%            i2=i2+1;
%            r(i2)=rr;
%        end
%    
%    end
    
%end
    





end

