function [Error] = Dua_Histogram2(data1,data2,lebar_bin,rmin,rmax)
%UNTITLED Summary of this function goes here

%data 2 disarankan dengan data model

hr1=histogram(data1);
hr1.BinWidth=lebar_bin;
data1_1 = hr1.BinEdges;
data1_2 = hr1.BinCounts;
%data3_r1 = h1.Values;
data1_1=data1_1(1:length(data1_1)-1); %mengurangi data yg kelebihan


hold on
hr2=histogram(data2);
hr2.BinWidth=lebar_bin;
data2_1 = hr2.BinEdges;
data2_2 = hr2.BinCounts;
data2_3 = hr2.Values;

hold off

data2_1=data2_1(1:length(data2_1)-1); %mengurangi data yg kelebihan


%error histogram

titik1=rmin;
titik2=rmin+lebar_bin;
x1=0;x2=0;
ke=1;

%untuk yang didalam ruang lingkup
while titik1>=rmin && titik2<=rmax
    
    %hitung banyak data1 pada bin ke-titik
    
    for i=1:length(data1)
        if data1(i)>=titik1 && data1(i)<=titik2
            x1=x1+1;
        end
    end
    
    for j=1:length(data2)
        if data2(j)>=titik1 && data2(j)<=titik2
           x2=x2+1; 
        end
    end
    
    %menghitung error
    rata2=(x1+x2)/2;
    xx1=(rata2-x1)^2;
    xx2=(rata2-x2)^2;
    sigma=sqrt((xx1+xx2)/1);
    err(ke)=sigma/sqrt(2);
   
    
    ke=ke+1;
    titik1=titik1+lebar_bin;
    titik2=titik2+lebar_bin;
    x1=0;x2=0;
end





%untuk yang diluar lingkup pada data 1
%untuk mendeteksi berapa banyak bin
x3=0;x4=0;
for i=1:length(data1)
   if  data1(i)<rmin
      x3=x3+1; 
      
   elseif data1(i)>rmax
      x4=x4+1; 
   end
end


titik1=rmin;
titik2=rmax;

%jika ada bin di luar lingkup


if x3>0 || x4>0
    while x3>0 || x4>0
        
    titik1=titik1-lebar_bin;
    titik2=titik2+lebar_bin;
    x31=0;
    x41=0;
    
    
    for i=1:length(data1)
        if data1(i)>=titik1 && data1(i)<titik1+lebar_bin %bagian kiri
           x31=x31+1; %ada binnya
           x3=x3-1;
        end
        
        if data1(i)<=titik2 && data1(i)>titik2-lebar_bin %bagian kanan
           x41=x41+1;
           x4=x4-1;
        end
        
    end
    
    if x31>0
        %menghitung error
    rata2=(x31+0)/2;
    xx1=(rata2-x31)^2;
    xx2=(rata2-0)^2;
    sigma=sqrt((xx1+xx2)/1);
    err(ke)=sigma/sqrt(2);
    ke=ke+1;
    end
    
    if x41>0
    rata2=(x41+0)/2;
    xx1=(rata2-x41)^2;
    xx2=(rata2-0)^2;
    sigma=sqrt((xx1+xx2)/1);
    err(ke)=sigma/sqrt(2);
    ke=ke+1;
    end
        

    end
end

Error=sum(err)/(length(err));





end

