function countobjects()

clear

imgPathes = {'Control', '1-DUIZHAO1_0001'; ...
    'CuSO4', '1-CuSO4-1';...
    'SBB', '1-SBB1'};

for i = 1:size(imgPathes,1)  % for each 3D image
    imgPathOri = ['./' imgPathes{i,1} '/' imgPathes{i,2} '.tif.frames'];
    [counto, zz, cc] = countImgObjects(imgPathes{i,2}, imgPathOri);
    % counto:bacterial load
    % zz:Proportional volume
    % cc:Percentage of average fluorescence intensity
end


function [counto, cc, zz] = countImgObjects(imgPath_i, imgPathOri)

% Calculating binarization threshold coefficient J 
M=0;
imgPathAll = dir([imgPathOri '/' imgPath_i '_C003Z*']);
for n=1:length(imgPathAll)
    imgPath = [imgPathOri '/' imgPathAll(n).name];
    I = imread(imgPath); 
    gray = rgb2gray(I);
    m=mean2(gray);
    M=M+m;
end
mean1=M/double(length(imgPathAll));
J=0.022*mean1+0.3861;%The fitted curve based on the three groups of average gray values and the debugged threshold
%Setting the length, width and height of the voxel
l=0.5101;w=0.5101;h=1;
%Setting Maximum diameter and minimum diameter
Md=5;md=0.5;
%% Loading Images
for n=1:length(imgPathAll)
    imgPath = [imgPathOri '/' imgPathAll(n).name];
    I = imread(imgPath); 
    s=size(I);
    gray = rgb2gray(I);

    %% Threshold segmentation and region growing method
    thresh = double(max(max(gray))) *J/255;
    bw = imbinarize(gray,thresh);
    L1 = logical(bw);
    STATS = regionprops(L1,'all');
    R=double(zeros(s(1),s(2)));
    for d =1:length(STATS)
        coor = STATS(d).Centroid;
        x0=coor(1);y0=coor(2);
        x0 = uint32(x0);y0 = uint32(y0);
        zhan = zeros(s(1)*s(2),2);
        pzhan = 1;
        zhan(pzhan,1)=x0; zhan(pzhan,2)=y0;
        R(y0,x0)=1;
        while pzhan > 0
           x1 = zhan(pzhan,1);y1 = zhan(pzhan,2);
            pzhan = pzhan - 1; 
            for a = -1 : 1
                for b = -1 : 1
                    if x1+a > 0 && x1+a <= s(2) &&y1+b > 0 && y1+b <=s(1) && gray(y1+b,x1+a)>gray(y0,x0)-15&& R(y1+b,x1+a) ~= 1
                        R(y1+b,x1+a) = 1;
                        pzhan = pzhan + 1;     
                        zhan(pzhan,1) = x1 + a;
                        zhan(pzhan,2) = y1 + b;
                    end
                end
            end
        end   
    end

    %% removing and filling 
    L2 = logical(R);
    STATS1 = regionprops(L2,'all');
    for d1 =1:length(STATS1)
        dia=STATS1(d1).MaxFeretDiameter;
        if dia*l>=100
            A=STATS1(d1).PixelList;
            for i=1:STATS1(d1).Area
                for j=1:2
                    B(j)=A(i,j);
                end
                R(B(2),B(1))=0;
            end
        end
    end
    r=imfill(R);

    %% creating volume
      Q=double(r).*double(gray);  
      model1(:,:,n)=Q;
      model2=model1;
end

%% Screening bacteria 
CC = bwconncomp(model1);
L3 = labelmatrix(CC);
stats = regionprops3(L3,model1,'all');
c=1;
z=0;
for v1=1:height(stats)
    if stats.BoundingBox(v1,4)*l<=Md&&stats.BoundingBox(v1,4)*l>=md...
            &&stats.BoundingBox(v1,5)*w<=Md&&stats.BoundingBox(v1,5)*w>=md...
            &&stats.BoundingBox(v1,6)*h<=Md&&stats.BoundingBox(v1,6)*h>=md
        V(c)=(stats.Volume(v1))*(l*w*h); 
         z=z+V(c);
        G(c)=stats.MeanIntensity(v1);
        c=c+1;
    else
        D=stats.VoxelList(v1);
        D1=cell2mat(D);
        for o=1:stats.Volume(v1)
            for p=1:3
                B(p)=D1(o,p);
            end
            model1(B(2),B(1),B(3))=0;
        end
    end
end
c=c-1;
%% counting mean fluorescence intensities and volume
C=zeros(1,9);
Z=zeros(1,4);
for k=1:c
    if V(k)>0&&V(k)<=4
        Z(1)=Z(1)+V(k);
    elseif V(k)>4&&V(k)<=8
        Z(2)=Z(2)+V(k);
    elseif V(k)>8&&V(k)<=12
        Z(3)=Z(3)+V(k);
    else
        Z(4)=Z(4)+V(k);
    end
end
zz=Z./z;  
for k=1:c
    if G(k)<=20
        C(1)=C(1)+1;
    elseif G(k)>20&&G(k)<=30
        C(2)=C(2)+1;
    elseif G(k)>30&&G(k)<=40
        C(3)=C(3)+1;
    elseif G(k)>40&&G(k)<=50
        C(4)=C(4)+1; 
    elseif G(k)>50&&G(k)<=60
        C(5)=C(5)+1; 
    elseif G(k)>60&&G(k)<=70
        C(6)=C(6)+1; 
    elseif G(k)>70&&G(k)<=80
        C(7)=C(7)+1; 
    elseif G(k)>80&&G(k)<=90
        C(8)=C(8)+1; 
    else
        C(9)=C(9)+1; 
     end
end
cc=C./c;
%% Counting bacterial load          
for o=1:length(imgPathAll)
    E=model1(:,:,o);
    bw1 = imbinarize(E,0);
    L4 = logical(bw1);
    STATS2 = regionprops(L4,'all');        
    counto(o)=length(STATS2);
end
                
              
              
        






