

clear all;
close all;
clc;
%% Read Original Image
Original_Image = imread('SunnyLake.bmp');
figure(1)
subplot(1,2,1)
imshow(Original_Image)
title('Original Image');

%% Manual GrayScale Generate

[w,h,z] = size(Original_Image);
% GrayScale_Image=zeros(w,h);
for(i=1:w)
    for(j=1:h)
        GrayScale_Image(i,j) = (Original_Image(i,j,1)+Original_Image(i,j,2)+Original_Image(i,j,1))/3;
    end
end
figure(1)
subplot(1,2,2)
imshow(GrayScale_Image);
title('GrayScale Image');

%% Manual Histogram
Histogram_Vector=zeros(1,256);
for(i=1:w)
    for(j=1:h)
Histogram_Vector((GrayScale_Image(i,j))+1) = Histogram_Vector((GrayScale_Image(i,j))+1)+1;
    end
end

for(i=1:256)
    if(Histogram_Vector(i)>5000)
        Histogram_Vector(i)=5000;
    end
end
        
y = linspace(1,256,256);

figure(2)
subplot(1,2,1)
imshow(GrayScale_Image);
title('GrayScale Image');
figure(2)
subplot(1,2,2)
plot(Histogram_Vector)
title('Histogram Graph');

%% OTSU thresholding

[counts,x] = imhist(GrayScale_Image,16);
stem(x,counts)
T = otsuthresh(counts);

BW = imbinarize(GrayScale_Image,T);
figure(3)
imshow(BW)
title('Binari Image');

%%  Add Gaussian Noise

NoisyImage_1 = imnoise(Original_Image,'gaussian',0,1);
figure(4)
subplot(2,2,1)
imshow(NoisyImage_1);
title('Noisy Image with Sigma 1');

NoisyImage_5 = imnoise(Original_Image,'gaussian',0,5);
figure(4)
subplot(2,2,2)
imshow(NoisyImage_5);
title('Noisy Image with Sigma 5');

NoisyImage_10 = imnoise(Original_Image,'gaussian',0,10);
figure(4)
subplot(2,2,3)
imshow(NoisyImage_10);
title('Noisy Image with Sigma 10');

NoisyImage_20 = imnoise(Original_Image,'gaussian',0,20);
figure(4)
subplot(2,2,4)
imshow(NoisyImage_20);
title('Noisy Image with Sigma 20');

%% Create GrayScale Image

for(i=1:w)
    for(j=1:h)
        GrayScale_Image_Noise1(i,j) = (NoisyImage_1(i,j,1)+NoisyImage_1(i,j,2)+NoisyImage_1(i,j,1))/3;
    end
end
figure(5) 
subplot(2,2,1)
imshow(GrayScale_Image_Noise1);
title('GrayScale Image with 1 Sigma Noisy');

for(i=1:w)
    for(j=1:h)
        GrayScale_Image_Noise2(i,j) = (NoisyImage_5(i,j,1)+NoisyImage_5(i,j,2)+NoisyImage_5(i,j,1))/3;
    end
end
figure(5) 
subplot(2,2,2)
imshow(GrayScale_Image_Noise2);
title('GrayScale Image with 5 Sigma Noisy');

for(i=1:w)
    for(j=1:h)
        GrayScale_Image_Noise3(i,j) = (NoisyImage_10(i,j,1)+NoisyImage_10(i,j,2)+NoisyImage_10(i,j,1))/3;
    end
end
figure(5) 
subplot(2,2,3)
imshow(GrayScale_Image_Noise3);
title('GrayScale Image with 10 Sigma Noisy');

for(i=1:w)
    for(j=1:h)
        GrayScale_Image_Noise4(i,j) = (NoisyImage_20(i,j,1)+NoisyImage_20(i,j,2)+NoisyImage_20(i,j,1))/3;
    end
end
figure(5) 
subplot(2,2,4) 
imshow(GrayScale_Image_Noise4);
title('GrayScale Image with 20 Sigma Noisy');

%% Create Low-Pass Filters

LP_Kernel1 = ones(3,3)/9;

LP_Kernel2 = ones(5,5)/25;

%% Denoising Noisy Image 1 (Sigma 1) wtih kernel1 (3x3)
t=1;
k=1;
 LP_OUT_Noise1_Kernel1 = zeros(w,h);
for(i=1:w-2)
    for(j=1:h-2)
 
        for(k=1:3)
            for(l=1:3)
            
                LP_OUT_Noise1_Kernel1(i+1,j+1) =  LP_OUT_Noise1_Kernel1(i+1,j+1) + GrayScale_Image_Noise1(k+i-1,l+j-1).*LP_Kernel1(k,l);
            end 
        end
        
    end 
end

figure(6) 
subplot(2,2,1) 
imshow(uint8(LP_OUT_Noise1_Kernel1))
title('Denoising Image For Sigma 1 Noise with LP 3x3 Kernel ');

%% Denoising Noisy Image 1 (Sigma 1) wtih kernel2 (5x5)
t=1;
k=1;
 LP_OUT_Noise1_Kernel2 = zeros(w,h);
for(i=1:w-4)
    for(j=1:h-4)
 
        for(k=1:5)
            for(l=1:5)
            
                LP_OUT_Noise1_Kernel2(i+1,j+1) =  LP_OUT_Noise1_Kernel2(i+1,j+1) + GrayScale_Image_Noise1(k+i-1,l+j-1).*LP_Kernel2(k,l);
            end 
        end
        
    end 
end

figure(6) 
subplot(2,2,2) 
imshow(uint8(LP_OUT_Noise1_Kernel2))
title('Denoising Image For Sigma 1 Noise with LP 5x5 Kernel ');



%% Denoising Noisy Image 2 (Sigma 5) wtih kernel1 (3x3)
t=1;
k=1;
 LP_OUT_Noise2_Kernel1 = zeros(w,h);
for(i=1:w-2)
    for(j=1:h-2)
 
        for(k=1:3)
            for(l=1:3)
            
                LP_OUT_Noise2_Kernel1(i+1,j+1) =  LP_OUT_Noise2_Kernel1(i+1,j+1) + GrayScale_Image_Noise2(k+i-1,l+j-1).*LP_Kernel1(k,l);
            end 
        end
        
    end 
end

figure(6) 
subplot(2,2,3) 
imshow(uint8(LP_OUT_Noise2_Kernel1))
title('Denoising Image For Sigma 5 Noise with LP 3x3 Kernel ');

%% Denoising Noisy Image 2 (Sigma 5) wtih kernel2 (5x5)
t=1;
k=1;
 LP_OUT_Noise2_Kernel2 = zeros(w,h);
for(i=1:w-4)
    for(j=1:h-4)
 
        for(k=1:5)
            for(l=1:5)
            
                LP_OUT_Noise2_Kernel2(i+1,j+1) =  LP_OUT_Noise2_Kernel2(i+1,j+1) + GrayScale_Image_Noise2(k+i-1,l+j-1).*LP_Kernel2(k,l);
            end 
        end
        
    end 
end

figure(6) 
subplot(2,2,4) 
imshow(uint8(LP_OUT_Noise2_Kernel2))
title('Denoising Image For Sigma 5 Noise with LP 5x5 Kernel ');




%% Denoising Noisy Image 3 (Sigma 10) wtih kernel1 (3x3)
t=1;
k=1;
 LP_OUT_Noise3_Kernel1 = zeros(w,h);
for(i=1:w-2)
    for(j=1:h-2)
 
        for(k=1:3)
            for(l=1:3)
            
                LP_OUT_Noise3_Kernel1(i+1,j+1) =  LP_OUT_Noise3_Kernel1(i+1,j+1) + GrayScale_Image_Noise3(k+i-1,l+j-1).*LP_Kernel1(k,l);
            end 
        end
        
    end 
end

figure(7) 
subplot(2,2,1) 
imshow(uint8(LP_OUT_Noise3_Kernel1))
title('Denoising Image For Sigma10 Noise with 3x3 LP Kernel ');

%% Denoising Noisy Image 3 (Sigma 10) wtih kernel2 (5x5)
t=1;
k=1;
 LP_OUT_Noise3_Kernel2 = zeros(w,h);
for(i=1:w-4)
    for(j=1:h-4)
 
        for(k=1:5)
            for(l=1:5)
            
                LP_OUT_Noise3_Kernel2(i+1,j+1) =  LP_OUT_Noise3_Kernel2(i+1,j+1) + GrayScale_Image_Noise3(k+i-1,l+j-1).*LP_Kernel2(k,l);
            end 
        end
        
    end 
end

figure(7) 
subplot(2,2,2) 
imshow(uint8(LP_OUT_Noise3_Kernel2))
title('Denoising Image For Sigma 10 Noise with 5x5 LP Kernel ');



%% Denoising Noisy Image 4 (Sigma 20) wtih kernel1 (3x3)
t=1;
k=1;
 LP_OUT_Noise4_Kernel1 = zeros(w,h);
for(i=1:w-2)
    for(j=1:h-2)
 
        for(k=1:3)
            for(l=1:3)
            
                LP_OUT_Noise4_Kernel1(i+1,j+1) =  LP_OUT_Noise4_Kernel1(i+1,j+1) + GrayScale_Image_Noise4(k+i-1,l+j-1).*LP_Kernel1(k,l);
            end 
        end
        
    end 
end

figure(7) 
subplot(2,2,3) 
imshow(uint8(LP_OUT_Noise4_Kernel1))
title('Denoising Image For Sigma20 Noise with 3x3 LP Kernel ');

%% Denoising Noisy Image4 (Sigma 20) wtih kernel2 (5x5)
t=1;
k=1;
 LP_OUT_Noise4_Kernel2 = zeros(w,h);
for(i=1:w-4)
    for(j=1:h-4)
 
        for(k=1:5)
            for(l=1:5)
            
                LP_OUT_Noise4_Kernel2(i+1,j+1) =  LP_OUT_Noise4_Kernel2(i+1,j+1) + GrayScale_Image_Noise4(k+i-1,l+j-1).*LP_Kernel2(k,l);
            end 
        end
        
    end 
end

figure(7) 
subplot(2,2,4) 
imshow(uint8(LP_OUT_Noise4_Kernel2))
title('Denoising Image For Sigma20 Noise with 5x5 LP Kernel ');


%% Create High-Pass Filters

HP_Kernel1 = [-1 -1 -1; -1 8 -1; -1 -1 -1];

HP_Kernel2 = [0.17 0.67 0.17; 0.67 -3.33 0.67; 0.17 0.67 0.17 ];

%% Denoising Noisy Image1 (Sigma 1) wtih HP kernel1 
t=1;
k=1;
 HP_OUT_Noise1_Kernel1 = zeros(w,h);
 HP_OUT_Noise1_Kernel2 = zeros(w,h);
for(i=1:w-2)
    for(j=1:h-2)
 
        for(k=1:3)
            for(l=1:3)
            
                HP_OUT_Noise1_Kernel1(i+1,j+1) =  HP_OUT_Noise1_Kernel1(i+1,j+1) + GrayScale_Image_Noise1(k+i-1,l+j-1).*HP_Kernel1(k,l);
                HP_OUT_Noise1_Kernel2(i+1,j+1) =  HP_OUT_Noise1_Kernel2(i+1,j+1) + GrayScale_Image_Noise1(k+i-1,l+j-1).*HP_Kernel2(k,l);
            
            end 
        end
        
    end 
end

figure(8)
subplot(2,1,1)
imshow(uint8(HP_OUT_Noise1_Kernel1))
title('Denoising Image For Sigma 1 Noise with HP Kernel 1 ');

subplot(2,1,2)
imshow(uint8(HP_OUT_Noise1_Kernel2))
title('Denoising Image For Sigma 1 Noise with HP Kernel 2 ');

%% Denoising Noisy Image2 (Sigma 5) wtih HP Kernel 1 and Kernel 2
t=1;
k=1;
 HP_OUT_Noise2_Kernel1 = zeros(w,h);
 HP_OUT_Noise2_Kernel2 = zeros(w,h);
for(i=1:w-2)
    for(j=1:h-2)
 
        for(k=1:3)
            for(l=1:3)
            
                HP_OUT_Noise2_Kernel1(i+1,j+1) =  HP_OUT_Noise2_Kernel1(i+1,j+1) + GrayScale_Image_Noise2(k+i-1,l+j-1).*HP_Kernel1(k,l);
                HP_OUT_Noise2_Kernel2(i+1,j+1) =  HP_OUT_Noise2_Kernel2(i+1,j+1) + GrayScale_Image_Noise2(k+i-1,l+j-1).*HP_Kernel2(k,l);
            
            end 
        end
        
    end 
end

figure(9)
subplot(2,1,1)
imshow(uint8(HP_OUT_Noise2_Kernel1))
title('Denoising Image For Sigma 5 Noise with HP Kernel 1 ');

subplot(2,1,2)
imshow(uint8(HP_OUT_Noise2_Kernel2))
title('Denoising Image For Sigma 5 Noise with HP Kernel 2 ');

%% Denoising Noisy Image3 (Sigma 10) wtih HP Kernel 1 and Kernel 2
t=1;
k=1;
 HP_OUT_Noise3_Kernel1 = zeros(w,h);
 HP_OUT_Noise3_Kernel2 = zeros(w,h);
for(i=1:w-2)
    for(j=1:h-2)
 
        for(k=1:3)
            for(l=1:3)
            
                HP_OUT_Noise3_Kernel1(i+1,j+1) =  HP_OUT_Noise3_Kernel1(i+1,j+1) + GrayScale_Image_Noise3(k+i-1,l+j-1).*HP_Kernel1(k,l);
                HP_OUT_Noise3_Kernel2(i+1,j+1) =  HP_OUT_Noise3_Kernel2(i+1,j+1) + GrayScale_Image_Noise3(k+i-1,l+j-1).*HP_Kernel2(k,l);
            
            end 
        end
        
    end 
end

figure(10)
subplot(2,1,1)
imshow(uint8(HP_OUT_Noise3_Kernel1))
title('Denoising Image For Sigma 10 Noise with HP Kernel 1 ');

subplot(2,1,2)
imshow(uint8(HP_OUT_Noise3_Kernel2))
title('Denoising Image For Sigma 10 Noise with HP Kernel 2 ');


%% Denoising Noisy Image4 (Sigma 20) wtih HP Kernel 1 and Kernel 2
t=1;
k=1;
 HP_OUT_Noise4_Kernel1 = zeros(w,h);
 HP_OUT_Noise4_Kernel2 = zeros(w,h);
for(i=1:w-2)
    for(j=1:h-2)
 
        for(k=1:3)
            for(l=1:3)
            
                HP_OUT_Noise4_Kernel1(i+1,j+1) =  HP_OUT_Noise4_Kernel1(i+1,j+1) + GrayScale_Image_Noise4(k+i-1,l+j-1).*HP_Kernel1(k,l);
                HP_OUT_Noise4_Kernel2(i+1,j+1) =  HP_OUT_Noise4_Kernel2(i+1,j+1) + GrayScale_Image_Noise4(k+i-1,l+j-1).*HP_Kernel2(k,l);
            
            end 
        end
        
    end 
end

figure(11)
subplot(2,1,1)
imshow(uint8(HP_OUT_Noise4_Kernel1))
title('Denoising Image For Sigma 20 Noise with HP Kernel 1 ');

subplot(2,1,2)
imshow(uint8(HP_OUT_Noise4_Kernel2))
title('Denoising Image For Sigma 20 Noise with HP Kernel 2 ');

%% Denoise Salt and Pepper Image

Noisy_SunnyLake = imread('SP_Noisy_SunnyLake.png');
RedChannel = Noisy_SunnyLake(:,:,1);
GreenChannel = Noisy_SunnyLake(:,:,2);
BlueChannel = Noisy_SunnyLake(:,:,3);

Denoise_R = medfilt2(RedChannel, [3 3]);
Denoise_G = medfilt2(GreenChannel, [3 3]);
Denoise_B = medfilt2(BlueChannel, [3 3]);

 %Reconstruct the noise free RGB%
Fin_denoise = cat(3,Denoise_R,Denoise_G,Denoise_B);

figure(12)
subplot(2,1,1)
imshow(Noisy_SunnyLake)
title('Noisy Image');

figure(12)
subplot(2,1,2)
imshow(Fin_denoise)
title('Denoised Image');


