function [sigmai] = sorting(data1)

cumu_freq=0;
total_data=sum(data1);

%dari terbesar ke terkecil
urutan_data=sort(data1,'descend');

for i=1:length(data1)
   
    jumlah_data(i)=cumu_freq+urutan_data(i);
    cumu_freq=jumlah_data(i);
    
    cumu_percent(i)=(cumu_freq/total_data)*100;
    
end

x5=1;
x16=1;
x84=1;
x95=1;


for i=1:length(cumu_percent)
   
    %phi 5%
    if x5==1 && cumu_percent(i)>=5
        pi5=(-log10(2*urutan_data(i))/log10(2));
        x5=0;
    end
    
    %phi 16%
    if x16==1 && cumu_percent(i)>=16
        pi16=(-log10(2*urutan_data(i))/log10(2));
        x16=0;
    end
    
    %phi 84%
    if x84==1 && cumu_percent(i)>=84
        pi84=(-log10(2*urutan_data(i))/log10(2));
        x84=0;
    end
    
     %phi 95%
    if x95==1 && cumu_percent(i)>=95
        pi95=(-log10(2*urutan_data(i))/log10(2));
        x95=0;
    end
    
    
end




%Inclusive graphic standart deviation

sigmai=((pi84-pi16)/4)+(pi95-pi5)/6.6;





end

