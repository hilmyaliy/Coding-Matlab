clc;
clear all;
close all;

location = 'batupasir';       %  folder in which your images exists
ds = imageDatastore(location);         %  Creates a datastore for all images in your folder
i=1;
while hasdata(ds) 
    grainasli(:,:,i) = read(ds);             % read image from datastore
    i=i+1;
end
i=i-1;

grainasli=imbinarize(grainasli);
[ukuran_gambar(1), ukuran_gambar(2), ukuran_gambar(3)]=size(grainasli);

[xx,yy] = meshgrid(1:ukuran_gambar(1));
zz=meshgrid(-70:30);

bin_radius=2;
bin_volume=200;
bin_luaspermukaan=50;


%figure(1)
%imagesc(grainasli(:,:,1))
%colormap(flipud(gray))
%axis equal; xlim([0 ukuran_gambar(1)]); ylim([0 ukuran_gambar(2)]);



%--------------------------------------------------------
daerah_DA                   = regionprops3(grainasli, 'Centroid', 'Volume','PrincipalAxisLength','SurfaceArea');
unikbutiran_DA              = unique(grainasli);
databutiran_DA             = unikbutiran_DA(unikbutiran_DA ~= 0);
%databutiran_DA             = databutiran_DA(databutiran_DA ~= 1);
%databutiran_DA             = databutiran_DA(databutiran_DA ~= 2); % Noise in 1-2?
satuan_butiran_DA          = length(databutiran_DA);
centroid_butiran_DA         = round(daerah_DA.Centroid);
volume_DA                   = round(daerah_DA.Volume);
luaspermukaan_DA            = round(daerah_DA.SurfaceArea);   
deret_citra2                    = length(grainasli(:));

jumlah_butiran_DA = length(volume_DA);

diameter1_DA=daerah_DA.PrincipalAxisLength;
for kk=1:length(diameter1_DA)
    diameter_DA(kk)=mean(diameter1_DA(kk,:));
end
radius_DA=diameter_DA/2;

figure(2)
subplot(3,1,1)
hd1=histogram(radius_DA); 
hd1.BinWidth=bin_radius;
xlabel('Radius (voksel)')
ylabel('Frekuensi')
title('Radius butiran Metode DA')
subplot(3,1,2)
hd2=histogram(volume_DA); 
hd2.BinWidth=bin_volume;
xlabel('Volume (voksel)')
ylabel('Frekuensi')
title('Volume butiran Metode DA')
subplot(3,1,3)
hd3=histogram(luaspermukaan_DA); 
hd3.BinWidth=bin_luaspermukaan;
xlabel('Luas Permukaan (voksel)')
ylabel('Frekuensi')
title('Luas Permukaan butiran Metode DA')


%-------------------------------------------------------------------
%proses watershed
D = bwdist(~grainasli);
%figure(4), isosurface(xx,yy,zz,D,4), axis equal
%title('Isosurface of distance transform')
%xlabel x, ylabel y, zlabel z
%xlim([-5 xmax+5]), ylim([-5 ymax+5]), zlim([-5 zmax+5])
%view(3), camlight,colormap('gray')

hold off

D = -D;
D(~grainasli) = Inf;

%-----
minThres=1;
D = imhmin(D,minThres);
%-----

water = watershed(D);
water(~grainasli) = 0;
%figure(5)
%for igrain=1:100
%isosurface(xx,yy,zz,water==igrain,0.5)
%hold on
%end

%axis equal
%title('Citra Watershed')
%xlabel x, ylabel y, zlabel z
%xlim([-5 100+5]), ylim([-5 100+5]), zlim([-70 30])
%view(3), camlight,colormap('gray')

%hold off

daerah_proses                   = regionprops3(water, 'Centroid', 'Volume','SurfaceArea');
unikbutiran_proses              = unique(water);
databutiran_proses             = unikbutiran_proses(unikbutiran_proses ~= 0);
databutiran_proses             = databutiran_proses(databutiran_proses ~= 1);
databutiran_proses             = databutiran_proses(databutiran_proses ~= 2); % Noise in 1-2?
satuan_butiran_proses          = length(databutiran_proses); %jumlah butiran proses
allGrainCentroid                = round(daerah_proses.Centroid(databutiran_proses,:));
volume_proses                   = round(daerah_proses.Volume(databutiran_proses,:));
deret_citra3                    = length(water(:));

luaspermukaan_proses=round(daerah_proses.SurfaceArea);


jumlah_butiran_proses=satuan_butiran_proses;    
    
[radius_proses,idxSingleGrain] = ukuranGrain3D(jumlah_butiran_proses,water,ukuran_gambar,deret_citra3,databutiran_proses,allGrainCentroid);

for kk=1:length(radius_proses(:,1))
   radius_proses_rata2(kk)=(sum(radius_proses(kk,:))/6); 
end


figure(6)
subplot(3,1,1)
hs1=histogram(radius_proses_rata2); 
hs1.BinWidth=bin_radius;
xlabel('Radius (voksel)')
ylabel('Frekuensi')
title('Radius butiran Metode PSA')
subplot(3,1,2)
hs2=histogram(volume_proses); 
hs2.BinWidth=bin_volume;
xlabel('Volume (voksel)')
ylabel('Frekuensi')
title('Volume butiran Metode PSA')
subplot(3,1,3)
hs3=histogram(luaspermukaan_proses); 
hs3.BinWidth=bin_luaspermukaan;
xlabel('Luas Permukaan (voksel)')
ylabel('Frekuensi')
title('Luas Permukaan butiran Metode PSA')
hold off


grainasli2= ~grainasli;

mkdir(['Butiran 3D citra asli(x y)']);
    for n = 1:size(grainasli2,3)
        tmp = squeeze(grainasli2(:,:,n));
        imwrite(tmp,['Butiran 3D citra asli(x y)/gambar_',num2str(n,'%0.4i'),'.bmp'],'bmp')
    end
    
mkdir(['Butiran 3D citra asli(x z)']);
    for n = 1:size(grainasli2,3)
        tmp = squeeze(grainasli2(:,n,:));
        imwrite(tmp,['Butiran 3D citra asli(x z)/gambar_',num2str(n,'%0.4i'),'.bmp'],'bmp')
    end    
    
    mkdir(['Butiran 3D citra asli(y z)']);
    for n = 1:size(grainasli2,3)
        tmp = squeeze(grainasli2(n,:,:));
        imwrite(tmp,['Butiran 3D citra asli(y z)/gambar_',num2str(n,'%0.4i'),'.bmp'],'bmp')
    end
    
    

water2=logical(water);   
water2= ~water2;
    
    
    mkdir(['watershed 3D citra asli(x y)']);
    for n = 1:size(water2,3)
        tmp = squeeze(water2(:,:,n));
        imwrite(tmp,['watershed 3D citra asli(x y)/gambar_',num2str(n,'%0.4i'),'.bmp'],'bmp')
    end
    
    mkdir(['watershed 3D citra asli(x z)']);
    for n = 1:size(water2,3)
        tmp = squeeze(water2(:,n,:));
        imwrite(tmp,['watershed 3D citra asli(x z)/gambar_',num2str(n,'%0.4i'),'.bmp'],'bmp')
    end

        mkdir(['watershed 3D citra asli(y z)']);
    for n = 1:size(water2,3)
        tmp = squeeze(water2(n,:,:));
        imwrite(tmp,['watershed 3D citra asli(y z)/gambar_',num2str(n,'%0.4i'),'.bmp'],'bmp')
    end










