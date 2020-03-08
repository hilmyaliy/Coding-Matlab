clc;
clear all;
close all;

%3 dimensi

grain=100;
rmin=10;
rmax=20;
vmin=((4/3)*pi*rmin^3);
vmax=((4/3)*pi*rmax^3);
lpmin=4*pi*rmin^2;
lpmax=4*pi*rmax^2;


xmax=200;
ymax=200;
zmax=200;

grain_persen=10;

bin_radius=4;
bin_volume=3000;
bin_luaspermukaan=100;

the_bw=zeros(xmax,ymax,zmax);

[the_bw,ukuran_gambar,volume_model,radius_model,luaspermukaan_model,xx,yy,zz,partikel,bw_lama]=GrainAsli3D(grain,rmax,rmin,xmax,ymax,zmax,grain_persen);

%gabungan semua bola
figure(3), isosurface(xx,yy,zz,the_bw), title('Citra model 3 dimensi')
axis equal
xlim([-5 xmax+5]), ylim([-5 ymax+5]), zlim([-5 zmax+5])
view(3), camlight,colormap('gray')
camlight
xlabel x, ylabel y, zlabel z


figure(4)
subplot(3,1,1)
hm1=histogram(radius_model); 
hm1.BinWidth=bin_radius;
xlabel('Radius (voksel)')
ylabel('Frekuensi')
title('Histogram Radius butiran Parameter Model')
subplot(3,1,2)
hm2=histogram(volume_model); 
hm2.BinWidth=bin_volume;
xlabel('Volume (voksel)')
ylabel('Frekuensi')
title('Histogram Volume butiran Parameter Model')
subplot(3,1,3)
hm3=histogram(luaspermukaan_model); 
hm3.BinWidth=bin_luaspermukaan;
xlabel('Luas Permukaan (voksel)')
ylabel('Frekuensi')
title('Histogram Luas Permukaan butiran Parameter Model')


%Direct Analysis
%------------------------------------------------------------
daerah_DA                   = regionprops3(the_bw, 'Centroid', 'Volume','PrincipalAxisLength','SurfaceArea');
unikbutiran_DA              = unique(the_bw);
databutiran_DA              = unikbutiran_DA(unikbutiran_DA ~= 0);
%databutiran_DA             = databutiran_DA(databutiran_DA ~= 1);
%databutiran_DA             = databutiran_DA(databutiran_DA ~= 2); % Noise in 1-2?
satuan_butiran_DA           = length(databutiran_DA);
centroid_butiran_DA         = round(daerah_DA.Centroid); %centroid
volume_DA                   = round(daerah_DA.Volume);
luaspermukaan_DA            = round(daerah_DA.SurfaceArea);
deret_citra2                = length(the_bw(:));

jumlah_butiran_DA = length(volume_DA);

diameter1_DA=daerah_DA.PrincipalAxisLength;
for kk=1:length(diameter1_DA)
    diameter_DA(kk)=mean(diameter1_DA(kk,:));
end
radius_DA=diameter_DA/2;

figure(5)
subplot(3,1,1)
hd1=histogram(radius_DA); 
hd1.BinWidth=bin_radius;
xlabel('Radius (voksel)')
ylabel('Frekuensi')
title('Histogram Radius butiran Metode DA')
subplot(3,1,2)
hd2=histogram(volume_DA); 
hd2.BinWidth=bin_volume;
xlabel('Volume (voksel)')
ylabel('Frekuensi')
title('Histogram Volume butiran Metode DA')
subplot(3,1,3)
hd3=histogram(luaspermukaan_DA); 
hd3.BinWidth=bin_luaspermukaan;
xlabel('Luas Permukaan (voksel)')
ylabel('Frekuensi')
title('Histogram Luas Permukaan butiran Metode DA')



%HILMY PROCEDUR

for i=1:length(partikel(1,1,1,:))

  %tanpa watershed
  daerah_HP = regionprops3(partikel(:,:,:,i), 'Centroid', 'Volume','SurfaceArea');
  unikbutiran_model = unique(partikel(:,:,:,i)); 
  data_butiran_HP         = unikbutiran_model(unikbutiran_model ~= 0); 
  satuan_butiran_HP(i)    = length(data_butiran_HP);
  
  centroid_HP      = round(daerah_HP.Centroid(data_butiran_HP,:));
  volume_HP(i)     = round(daerah_HP.Volume(data_butiran_HP,:)); 
  luaspermukaan_HP(i) = round(daerah_HP.SurfaceArea);
  
  luas_citra1(:,:,:)    = partikel(:,:,:,i);
  deret_citra1          = length(luas_citra1(:)); 
  

    R_model           = zeros(satuan_butiran_HP(i),6); %jari2
    grainAzimuth2     = zeros(satuan_butiran_HP(i),6); 
    grainInclination2 = zeros(satuan_butiran_HP(i),6);  
  

jumlah_butiran_HP(i)=satuan_butiran_HP(i);
%    #3

[R_model] = ukuranGrain3D(jumlah_butiran_HP,partikel(:,:,:,i),ukuran_gambar,deret_citra1,data_butiran_HP,centroid_HP);

radius_HP(i,1:6)=R_model;


hold off

    
end

jumlah_butiran_HP=sum(satuan_butiran_HP(:));


for kk=1:length(radius_HP(:,1))
   radius_HP_rata2(kk)=(sum(radius_HP(kk,:))/6); 
end


figure(6)
subplot(3,1,1)
hp1=histogram(radius_HP_rata2);
hp1.BinWidth=bin_radius;
xlabel('Radius (voksel)')
ylabel('Frekuensi')
title('Histogram Radius Rata2 butiran Metode NP')

subplot(3,1,2)
hp2=histogram(volume_HP);
hp2.BinWidth=bin_volume;
xlabel('Volume (voksel)')
ylabel('Frekuensi')
title('Histogram Volume butiran Metode NP')

subplot(3,1,3)
hp3=histogram(luaspermukaan_HP); 
hp3.BinWidth=bin_luaspermukaan;
xlabel('Luas Permukaan (voksel)')
ylabel('Frekuensi')
title('Histogram Luas Permukaan butiran Metode NP')


%Error histogram HP dan Model
figure(7)
subplot(3,1,1)
[error_HPr]=Dua_Histogram2(radius_HP_rata2,radius_model,bin_radius,rmin,rmax);
xlabel('Radius (voksel)')
ylabel('Frekuensi')
title('Histogram Radius butiran antara Parameter Model dan Metode NP')
legend('Radius rata2 butiran Metode NP','Radius butiran Parameter Model')

subplot(3,1,2)
[error_HPv]=Dua_Histogram2(volume_HP,volume_model,bin_volume,vmin,vmax);
xlabel('Volume (voksel)')
ylabel('Frekuensi')
title('Histogram Volume butiran antara Parameter Model dan Metode NP')
legend('Volume butiran Metode NP','Volume butiran Model')

subplot(3,1,3)
[error_HPlp]=Dua_Histogram2(luaspermukaan_HP,luaspermukaan_model,bin_luaspermukaan,lpmin,lpmax);
xlabel('Luas Permukaan (voksel)')
ylabel('Frekuensi')
title('Histogram Luas Permukaan butiran antara Parameter Model dan Metode NP')
legend('Luas Permukaan butiran Metode NP','Luas Permukaan butiran Parameter Model')

%radius_DA = mean(diameter_DA())
%for kk=1:length(radius_DA)
%radius_DA_rata2(kk)=sum(radius_DA(kk,:))/2;
%end

%-------------------------------------------------------------------
%proses watershed
D = bwdist(~the_bw);

figure(8), isosurface(xx,yy,zz,D,4), axis equal
title('Isosurface of distance transform')
xlabel x, ylabel y, zlabel z
xlim([-5 xmax+5]), ylim([-5 ymax+5]), zlim([-5 zmax+5])
view(3), camlight,colormap('gray')

D = -D;
D(~the_bw) = Inf;

%-----
minThres=1;
D = imhmin(D,minThres);
%-----

water = watershed(D);
water(~the_bw) = 0;

figure(9)
for igrain=1:grain
isosurface(xx,yy,zz,water==igrain,0.5)
hold on
end

axis equal
title('Citra Watershed')
xlabel x, ylabel y, zlabel z
xlim([-5 xmax+5]), ylim([-5 ymax+5]), zlim([-5 zmax+5])
view(3), camlight,colormap('gray')

daerah_proses                   = regionprops3(water, 'Centroid', 'Volume','SurfaceArea');
unikbutiran_proses              = unique(water);
databutiran_proses             = unikbutiran_proses(unikbutiran_proses ~= 0);
databutiran_proses             = databutiran_proses(databutiran_proses ~= 1);
databutiran_proses             = databutiran_proses(databutiran_proses ~= 2); % Noise in 1-2?
satuan_butiran_proses          = length(databutiran_proses); %jumlah butiran proses
allGrainCentroid                = round(daerah_proses.Centroid(databutiran_proses,:));

volume_proses                   = round(daerah_proses.Volume(databutiran_proses,:));
luaspermukaan_proses           = round(daerah_proses.SurfaceArea);
deret_citra3                    = length(water(:));



jumlah_butiran_proses=satuan_butiran_proses;    
    
[radius_proses,idxSingleGrain] = ukuranGrain3D(jumlah_butiran_proses,water,ukuran_gambar,deret_citra3,databutiran_proses,allGrainCentroid);

for kk=1:length(radius_proses(:,1))
   radius_proses_rata2(kk)=(sum(radius_proses(kk,:))/6); 
end


figure(10)
subplot(3,1,1)
hs1=histogram(radius_proses_rata2); %kuning
hs1.BinWidth=bin_radius;
xlabel('Radius (voksel)')
ylabel('Frekuensi')
title('Histogram Radius butiran Metode PSA')
subplot(3,1,2)
hs2=histogram(volume_proses); %kuning
hs2.BinWidth=bin_volume;
xlabel('Volume (voksel)')
ylabel('Frekuensi')
title('Histogram Volume butiran Metode PSA')
subplot(3,1,3)
hs3=histogram(luaspermukaan_proses); 
hs3.BinWidth=bin_luaspermukaan;
xlabel('Luas Permukaan (voksel)')
ylabel('Frekuensi')
title('Histogram Luas Permukaan butiran Metode PSA')
hold off


%mencari error dua histogram
figure(11)
subplot(3,1,1)
[error_DAr]=Dua_Histogram2(radius_DA,radius_model,bin_radius,rmin,rmax);
xlabel('Radius (voksel)')
ylabel('Frekuensi')
title('Histogram Radius butiran antara Parameter Model dan Metode DA')
legend('Radius rata2 butiran Metode DA','radius butiran Parameter Model')

subplot(3,1,2)
[error_DAv]=Dua_Histogram2(volume_DA,volume_model,bin_volume,vmin,vmax);
xlabel('Volume (voksel)')
ylabel('Frekuensi')
title('Histogram Volume butiran antara Parameter Model dan Metode DA')
legend('Volume butiran Metode DA','Volume butiran Parameter Model')

subplot(3,1,3)
[error_DAlp]=Dua_Histogram2(luaspermukaan_DA,luaspermukaan_model,bin_luaspermukaan,lpmin,lpmax);
xlabel('Luas Permukaan (voksel)')
ylabel('Frekuensi')
title('Histogram Luas Permukaan butiran antara Parameter Model dan Metode DA')
legend('Luas Permukaan butiran Metode DA','Luas Permukaan butiran Parameter Model')


figure(12)
subplot(3,1,1)
[error_prosesr]=Dua_Histogram2(radius_proses_rata2,radius_model,bin_radius,rmin,rmax);
xlabel('Radius (voksel)')
ylabel('Frekuensi')
title('Histogram Radius butiran antara Parameter Model dan Metode PSA')
legend('Radius rata2 butiran Metode PSA','Radius butiran Parameter Model')

subplot(3,1,2)
[error_prosesv]=Dua_Histogram2(volume_proses,volume_model,bin_volume,vmin,vmax);
xlabel('Volume (voksel)')
ylabel('Frekuensi')
title('Histogram Volume butiran antara Parameter Model dan Metode PSA')
legend('Volume butiran Metode PSA','Volume butiran Parameter Model')

subplot(3,1,3)
[error_proseslp]=Dua_Histogram2(luaspermukaan_proses,luaspermukaan_model,bin_luaspermukaan,lpmin,lpmax);
xlabel('Luas Permukaan (voksel)')
ylabel('Frekuensi')
title('Histogram Luas Permukaan butiran antara Parameter Model dan Metode PSA')
legend('Luas Permukaan butiran Metode PSA','Luas Permukaan butiran Parameter Model')





%POROSITAS

%poros_model=((xmax*ymax*zmax-sum(volume_model))/(xmax*ymax*zmax)*100);
poros_DA=((xmax*ymax*zmax-sum(volume_DA))/(xmax*ymax*zmax)*100);
poros_HP=((xmax*ymax*zmax-sum(volume_HP))/(xmax*ymax*zmax)*100);
poros_proses=((xmax*ymax*zmax-sum(volume_proses))/(xmax*ymax*zmax)*100);


%Rata2 luas dan radius
rata2_rmodel=sum(radius_model)/length(radius_model);
rata_rDA=sum(radius_DA)/length(radius_DA);
rata2_rHP=sum(radius_HP_rata2)/length(radius_HP_rata2);
rata2_rPS=sum(radius_proses_rata2)/length(radius_proses_rata2);

rata2_vPS=sum(volume_proses)/length(volume_proses);
rata2_vmodel=sum(volume_model)/length(volume_model);
rata2_vDA=sum(volume_DA)/length(volume_DA);
rata2_vHP=sum(volume_HP)/length(volume_HP);

rata2_lpPS=sum(luaspermukaan_proses)/length(luaspermukaan_proses);
rata2_lpmodel=sum(luaspermukaan_model)/length(luaspermukaan_model);
rata2_lpDA=sum(luaspermukaan_DA)/length(luaspermukaan_DA);
rata2_lpHP=sum(luaspermukaan_HP)/length(luaspermukaan_HP);


%sorting
[sigma_model]=sorting(radius_model);
[sigma_DA]=sorting(radius_DA);
[sigma_HP]=sorting(radius_HP_rata2);
[sigma_PS]=sorting(radius_proses_rata2);




%%

the_bw2= ~the_bw;
mkdir(['Butiran 3D (x y)']);
    for n = 1:size(the_bw2,3)
        tmp = squeeze(the_bw2(:,:,n));
        imwrite(tmp,['Butiran 3D (x y)/gambar_',num2str(n,'%0.4i'),'.bmp'],'bmp')
    end
    
mkdir(['Butiran 3D (x z)']);
    for n = 1:size(the_bw2,3)
        tmp = squeeze(the_bw2(:,n,:));
        imwrite(tmp,['Butiran 3D (x z)/gambar_',num2str(n,'%0.4i'),'.bmp'],'bmp')
    end    
    
    mkdir(['Butiran 3D (y z)']);
    for n = 1:size(the_bw2,3)
        tmp = squeeze(the_bw2(n,:,:));
        imwrite(tmp,['Butiran 3D (y z)/gambar_',num2str(n,'%0.4i'),'.bmp'],'bmp')
    end
    
    

water2=logical(water);   
water2= ~water2;
    
    
    mkdir(['watershed 3D (x y)']);
    for n = 1:size(water2,3)
        tmp = squeeze(water2(:,:,n));
        imwrite(tmp,['watershed 3D (x y)/gambar_',num2str(n,'%0.4i'),'.bmp'],'bmp')
    end
    
    mkdir(['watershed 3D (x z)']);
    for n = 1:size(water2,3)
        tmp = squeeze(water2(:,n,:));
        imwrite(tmp,['watershed 3D (x z)/gambar_',num2str(n,'%0.4i'),'.bmp'],'bmp')
    end

        mkdir(['watershed 3D (y z)']);
    for n = 1:size(water2,3)
        tmp = squeeze(water2(n,:,:));
        imwrite(tmp,['watershed 3D (y z)/gambar_',num2str(n,'%0.4i'),'.bmp'],'bmp')
    end


