clc;
close all;
clear all;



%input banyaknya grain
butiran_input=100;
%batas radius
rmax=15;
rmin=5;
%bingkai
xmax=300; 
ymax=300;
dimensi=2.5;
tingkat_overlap=5; % toleransi overlap (jika mau 0 persen, maka diinput antara 0<x<1)
bin_radius=2;
bin_luas=100;

%--------------------------------------------------------------------------


tic

[citra_model,ukuran_citra,radius_model,volume_model,volume_HP,partikel_HP,butiran_model,diameter_model] = GrainAsli2(butiran_input,rmax,rmin,xmax,ymax,tingkat_overlap); %GRAIN ASLI
citra_model(1,1)=0;


for k=1:length(radius_model)
Phi(k)=(-log10(2*radius_model(k))/log10(2));
end
figure(2)
hp=histogram(Phi);
hp.BinWidth=1;


%porositas
poros_model=(sum(volume_model)/(xmax*ymax)*100);

vmin=min(volume_model);
vmax=max(volume_model);

%Direct Analyze Model
daerah_DA=regionprops('table',citra_model,'Centroid','Area','MajorAxisLength','MinorAxisLength');
centroid_DA=round(daerah_DA.Centroid);

%gambar dengan centroid
figure(5)
imagesc(citra_model)
colormap(flipud(gray))
axis equal; xlim([0 ukuran_citra(2)]); ylim([0 ukuran_citra(1)]);
hold on
plot(centroid_DA(:,1),centroid_DA(:,2),'*r')
title('Citra dengan Metode DA')
hold off

volume_DA=daerah_DA.Area;

radius_DA=[daerah_DA.MajorAxisLength./2 daerah_DA.MinorAxisLength./2];
for kk=1:length(radius_DA(:,1))
radius_DA_rata2(kk)=sum(radius_DA(kk,:))/2;
end


figure(6)
subplot(1,2,1)
hd1=histogram(radius_DA_rata2);
hd1.BinWidth=bin_radius;
nilaibin_DA1=hd1.Values;

nilaibin_DA2=hd1.BinEdges;
nilaibin_DA2=nilaibin_DA2(nilaibin_DA2 ~= 0);


xlabel('Radius (piksel)')
ylabel('Frekuensi')
title('Histogram Radius butiran Metode DA')


subplot(1,2,2)
hd2=histogram(volume_DA);
hd2.BinWidth=bin_luas;

xlabel('Luas (piksel)')
ylabel('Frekuensi')
title('Histogram Luas butiran Metode DA')





jumlah_butiran_DA=length(centroid_DA);

figure(7)
subplot(2,1,1)
[error_DA_model_r]=Dua_Histogram2(radius_DA_rata2,radius_model,bin_radius,rmin,rmax);
xlabel('Radius (piksel)')
ylabel('Frekuensi')
title('Histogram Radius butiran antara Parameter Model dan Metode DA')
legend('Radius rata2 butiran Metode DA','Radius butiran Parameter Model')


subplot(2,1,2)
[error_DA_model_v]=Dua_Histogram2(volume_DA,volume_model,bin_luas,vmin,vmax);
xlabel('Luas (piksel)')
ylabel('Frekuensi')
title('Histogram Luas butiran antara Parameter Model dan Metode DA')
legend('Luas butiran Metode DA','Luas butiran Parameter Model')


%   #1
figure(9)
subplot(1,2,1)
histogram(radius_model)
xlabel('Radius (piksel)')
ylabel('Frekuensi')
title('Histogram Radius butiran Parameter Model')
subplot(1,2,2)
histogram(volume_model)
xlabel('Luas (piksel)')
ylabel('Frekuensi')
title('Histogram Luas butiran Parameter Model')


figure(10)
imagesc(citra_model)
axis equal
xlim([0 ukuran_citra(2)]); ylim([0 ukuran_citra(1)]);
title('Citra model')
colormap(flipud(gray))
hold off

%PCA Hilmy PROCEDURE----------------------------------------------------
ii=0;
figure(11)
imagesc(citra_model)
colormap(flipud(gray))
axis equal; xlim([0 ukuran_citra(2)]); ylim([0 ukuran_citra(1)]);
title('Citra Metode NP')

partikel_HP=logical(partikel_HP);
for i=1:length(partikel_HP(1,1,:))
    ii=ii+1;
    partikel_HP(1,1,i)=0;
  daerah_model = regionprops('table',partikel_HP(:,:,i), 'Centroid', 'Area');
  unikbutiran_model = unique(partikel_HP(:,:,i)); 
  data_butiran_model         = unikbutiran_model(unikbutiran_model ~= 0); 
  satuan_butiran_HP(i)    = length(data_butiran_model);
  
  centroid_model      = round(daerah_model.Centroid(data_butiran_model,:));
  %volume_model(i)     = round(daerah_model.Area(data_butiran_model,:)); 
  
  luas_citra1(:,:)      = partikel_HP(:,:,i);
  deret_citra1          = length(luas_citra1(:)); 
  volume2_HP(i)         = round(daerah_model.Area);
if ukuran_citra(3) == 1     % 2D------------------------------------------
    R_model           = zeros(satuan_butiran_HP(i),4); %jari2
    grainAzimuth2     = zeros(satuan_butiran_HP(i),4); 
    grainInclination2 = zeros(satuan_butiran_HP(i),4);  
  
end
%    #3

[R_model] = ukuranGrain(satuan_butiran_HP(i),data_butiran_model,deret_citra1,ukuran_citra,centroid_model,partikel_HP(:,:,i));

radius_HP(i,1:4)=R_model;

hold off
end


jumlah_butiran_HP=sum(satuan_butiran_HP(:));


for kk=1:length(radius_HP(:,1))
   radius_HP_rata2(kk)=(sum(radius_HP(kk,:))/4); 
end


figure(12)
subplot(1,2,1)
histogram(radius_HP_rata2)
xlabel('Radius (piksel)')
ylabel('Frekuensi')
title('Histogram Radius Rata2 butiran Metode NP')
subplot(1,2,2)
histogram(volume2_HP)
xlabel('Luas (piksel)')
ylabel('Frekuensi')
title('Histogram Luas butiran Metode NP')


figure(13)
subplot(2,1,1)
[error_HP_model_r]=Dua_Histogram2(radius_HP_rata2,radius_model,bin_radius,rmin,rmax);
xlabel('Radius (piksel)')
ylabel('Frekuensi')
title('Histogram Radius butiran antara Paremeter Model dan Metode NP')
legend('Radius rata2 butiran Metode NP','Radius butiran Parameter Model')

subplot(2,1,2)
[error_HP_model_v]=Dua_Histogram2(volume_HP,volume_model,bin_luas,vmin,vmax);
xlabel('Luas (piksel)')
ylabel('Frekuensi')
title('Histogram Luas butiran antara Parameter Model dan Metode NP')
legend('Luas butiran Metode NP','Luas butiran Parameter Model')

%-------------------------------------------------------------------------



%----------------------------PROSES DIMULAI !!!---------------------------

minthres=linspace(1,9,3);
%minthres=2;

ket1=2;
ket2=2;
k=1;
plott=13;
for i=1:length(minthres)
%WATERSHED
                                                               %input batas minimum thresholding
[water_proses(:,:,i)]=grainwatershed(citra_model,minthres(i));

daerah_proses = regionprops('table', water_proses(:,:,i), 'Centroid', 'Area');             %menemukan luas dan pusat massa grain

unikbutiran_proses = unique(water_proses(:,:,i));                                            %mengambil 1 nilai (unik) dari setiap sel dan diurutkan
                                                                            %hasil = mendapatkan grain dari 1 nilai pada citra gray 8 bit                                                 
% karasteristik grain (tidak termasuk 0)                                    %hapus data yang bukan grain,didapat data grain tanpa pori dan noise
                                                                                                                                            
data_butiran_proses          = unikbutiran_proses(unikbutiran_proses ~= 0);                            % menghilangkan nilai 0
data_butiran_proses          = data_butiran_proses(data_butiran_proses ~= 1);                          % menghilangkan noise, nilai data grain minimal 3 atau 2 tiap sel
%data_grain          = data_grain(data_grain ~= 2);

jumlah_butiran_proses(i)     = length(data_butiran_proses);                                  %jumlah partikel

centroid_proses    = round(daerah_proses.Centroid(data_butiran_proses,:));                 %pusat tiap partikel
volume_proses      = round(daerah_proses.Area(data_butiran_proses,:));                     %volume tiap partikel

luas_citra2(:,:)  = water_proses(:,:,i);
deret_citra2      = length(luas_citra2(:));                               %luas/volume gambar

%membuat ruang kosong
if ukuran_citra(3) == 1     % 2D------------------------------------------
    radius_proses    = zeros(jumlah_butiran_proses(i),4); %jari2
    grainAzimuth     = zeros(jumlah_butiran_proses(i),4); 
    grainInclination = zeros(jumlah_butiran_proses(i),4);
elseif ukuran_citra(3) > 1  % 3D------------------------------------------
    radius_proses      = zeros(jumlah_butiran_proses(i),6);
    grainAzimuth     = zeros(jumlah_butiran_proses(i),6);
    grainInclination = zeros(jumlah_butiran_proses(i),6);
end

%porositas
poros_proses(i)=((xmax*ymax-sum(volume_proses))/(xmax*ymax)*100);



figure(15)
subplot(1,length(minthres),i)
imagesc(water_proses(:,:,i))
colormap(flipud(gray))

title(['h-min = ' num2str(minthres(i))])
axis equal; xlim([0 ukuran_citra(2)]); ylim([0 ukuran_citra(1)]);


figure(16)
subplot(1,length(minthres),i)
imagesc(water_proses(:,:,i))
colormap(flipud(gray))

title(['h-min = ' num2str(minthres(i))])

axis equal; xlim([0 ukuran_citra(2)]); ylim([0 ukuran_citra(1)]);

[radius_proses] = ukuranGrain(jumlah_butiran_proses(i),data_butiran_proses,deret_citra2,ukuran_citra,centroid_proses,water_proses(:,:,i));
hold off


%mencari rata-rata grain dari radius analisis
for kk=1:length(radius_proses(:,1))
   radius_proses_rata2(kk)=(sum(radius_proses(kk,:))/4); 
end

[sigma_PS(i)]=sorting(radius_proses_rata2);


figure(17)
subplot(2,length(minthres),i)
[error_PSA_model_r(i)]=Dua_Histogram2(radius_proses_rata2,radius_model,bin_radius,rmin,rmax);
xlabel('Radius (piksel)')
ylabel('Frekuensi')

if i==ket1
title({'Histogram Radius butiran antara Parameter Model dan Metode PSA',['h-min = ' num2str(minthres(i))]})
legend('Radius rata2 butiran Metode PSA','Radius butiran Parameter Model')

else
    title(['h-min = ' num2str(minthres(i))])
end



%--------------------------------------------------------------------------

%VOLUME/LUAS butiran

%figure(18)
subplot(2,length(minthres),length(minthres)+i)
[error_PSA_model_v(i)]=Dua_Histogram2(volume_proses,volume_model,bin_luas,vmin,vmax);
xlabel('Luas (piksel)')
ylabel('Frekuensi')


if i==ket1
title({'Histogram Luas butiran antara Parameter Model dan Metode PSA',['h-min = ' num2str(minthres(i))]})
legend('Luas butiran Metode PSA','Luas butiran Parameter Model')

else
    title(['h-min = ' num2str(minthres(i))])
end


%histogram radius dan volume PSA
figure(19)
subplot(2,length(minthres),i)
histogram(radius_proses_rata2)
xlabel('Radius (piksel)')
ylabel('Frekuensi')
if i==ket1
title('Histogram Radius rata2 butiran Metode PSA')

end


subplot(2,length(minthres),length(minthres)+i)
histogram(volume_proses)
xlabel('Luas (piksel)')
ylabel('Frekuensi')
if i==ket1
title('Histogram Luas butiran Metode PSA')

end




%------------------------------------
% memperbaiki masalah jumlah grain telalu banyak
if jumlah_butiran_proses(i)>butiran_input
    jumlah_butiran_proses(i)=butiran_input;
end


k=k+1;

rata2_rPS(i)=sum(radius_proses_rata2)/length(radius_proses);
rata2_vPS(i)=sum(volume_proses)/length(volume_proses);

end

    %8
    %error radius
figure(20)
subplot(2,2,1)
plot(minthres,error_PSA_model_r,'*')
title('Error Radius butiran')
xlabel('H-minimum')
ylabel('Error (%)')

%error volume
subplot(2,2,2)
plot(minthres,error_PSA_model_v,'*')
title('Error Luas butiran')
xlabel('H-minimum')
ylabel('Error (%)')


%Error Selisih Partikel
for i=1:length(minthres)
   GSR(i)=(abs(butiran_input-jumlah_butiran_proses(i))/butiran_input)*100; 
end

subplot(2,2,[3 4])
plot(minthres,GSR,'*')
title('Nilai Error selisih jumlah Butiran')
xlim([0 minthres(length(minthres))])
xlabel('H-minimum')
ylabel('%')

hold off


%POROSITAS

poros_model=((xmax*ymax-sum(volume_model))/(xmax*ymax)*100);
poros_DA=((xmax*ymax-sum(volume_DA))/(xmax*ymax)*100);
poros_HP=((xmax*ymax-sum(volume2_HP))/(xmax*ymax)*100);

%Rata2 luas dan radius

rata2_rmodel=sum(radius_model)/length(radius_model);
rata_rDA=sum(radius_DA_rata2)/length(radius_DA_rata2);
rata2_rHP=sum(radius_HP_rata2)/length(radius_HP_rata2);


rata2_vmodel=sum(volume_model)/length(volume_model);
rata2_vDA=sum(volume_DA)/length(volume_DA);
rata2_vHP=sum(volume2_HP)/length(volume2_HP);



%sorting

%sorting
[sigma_model]=sorting(radius_model);
[sigma_DA]=sorting(radius_DA_rata2);
[sigma_HP]=sorting(radius_HP_rata2);



toc



%mencari rata-rata grain dari radius analisis





%ddata_grain=double(data_grain);

%Metode


% Plot Cutoff







%[histFB, binCenter, statGSD ] = hist_freq( grainDiameter );

%figure(121)
%[hf] = GrainSizeDistribution(binCenter, histFB);
%title('Metode Laser Diffraction')


%[histVB, binCenter, statGSD] = hist_vol( grainDiameter, allGrainVolume );

%figure(122)
%[hv] = GrainSizeDistribution(binCenter, histVB);
%title('Metode Sieve')

%[histPC, binCenter,maxPixel, statGSD] = point_count( grainDiameter, allGrainCentroid );
%figure(133)

%[hp] = GrainSizeDistribution(binCenter, histPC);
%title('Metode Point-Count')
%toc



%imwrite(gambar,'cobagambar.bmp','bmp');


%if there's something wrong with the frame (bingkai), use this
%cla reset;

%subplot(1,2,2)
%hist(grainRadius)
%title('Histrogram Radius')

%figure(10)

%subplot(1,2,2)
%hist(allGrainVolume)
%title('Setelah Watershed')

%figure(11)
%hist(r)

%10

%figure(20)
%[Dgrain]=DiameterGrain(radius_proses,"mean");




%grainRadiuss=radius_proses*0.1;
%Dgrain=Dgrain*0.01;
%hist(grainRadiuss)

%cutoffNum = [2 1 0.5 0.25 0.125 0.0625 0.0039 0];
%cutoffTxt = {'Kersik','Pasir Sangat Kasar','Pasir Kasar','Pasir Medium','Pasir Halus','Pasir Sangat Halus','Lanau','Tanah Liat'};
%cutoffTxtLoc = [2.4 1.2 0.6 0.3 0.15 0.075 0.0048 0.003];
%nTxt = length(cutoffNum);


%box on;
%hold on


%for iTxt = 1:nTxt
%    plot([log10(cutoffNum(iTxt)) log10(cutoffNum(iTxt))],[0 1],'--', 'color', [0.6350 0.0780 0.1840])
%    h = text(log10(cutoffTxtLoc(iTxt)),0.25, cutoffTxt{iTxt});
%    set(h, 'rotation', 90); h(1).Color = [0.6350 0.0780 0.1840];
    %hasilakhir = getframe(h);
    %hold on
%end
%hold off


     
  
   
   
   
   



