clear all
Picture = imread('camera.pgm');

Picture_Gray = Picture;%灰度处理

Picture_FFT = fft2(Picture_Gray);%傅里叶变换
Picture_FFT_Shift = fftshift(Picture_FFT);%对频谱进行移动，是0频率点在中心
Picture_AM_Spectrum = log(abs(Picture_FFT_Shift));%获得傅里叶变换的幅度谱
Picture_Phase_Specture = log(angle(Picture_FFT_Shift)*180/pi);%获得傅里叶变换的相位谱
Picture_Restructure = ifft2(abs(Picture_FFT).*exp(j*(angle(Picture_FFT))));%双谱重构
% figure(1)
% subplot(221)
% imshow(Picture_Gray)
% title('原图像')
% subplot(222)
% imshow(Picture_AM_Spectrum,[])%显示图像的幅度谱，参数'[]'是为了将其值线性拉伸
% title('图像幅度谱')
% subplot(223)
% imshow(Picture_Phase_Specture,[]);
% title('图像相位谱')
% subplot(224)
% imshow(Picture_Restructure,[]);
% title('双谱重构图')
% Picture_Restructure = ifft2(abs(Picture_FFT));%幅度谱重构
% imshow(Picture_Restructure,[]);
Picture_Restructure = ifft2(exp(j*(angle(Picture_FFT))));%相位谱重构
imshow(Picture_Restructure,[]);
imwrite(Picture_Restructure,'camaraIDFTPhaseAngle.pgm')