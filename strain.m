
%***********************************************************************************************%
% Author: Shuhui Liu
% E-mail: liushuhuibpt@126.com
%***********************************************************************************************%
%compute the strain

clear;
clc;
close all;
file='300O2';   
pathname = 'D:\github\input\';
filename = strcat(file,'.xlsx');%file of atom coordinates
pathname1 = 'D:\github\input\';%raw image
filename1 = strcat(file,'.tif');
coords = readmatrix(strcat(pathname,filename));
f = imread(strcat(pathname1,filename1));
coefficients=zeros(26,2);
figure,imshow(f,[]);
selected_indices = coords(:, 1) > 718 ;
coords = coords(selected_indices, :);%Remove the recognition coordinates that are not in the crystal area.
hold on
x=coords(:,2);
y=1350-coords(:,1);
scatter(x, 1350-y, 'filled');
hold on;
%Calculate the reference distance of this image
selected_indices = coords(:, 1) > 1151 & coords(:, 1) < 1297 & coords(:, 2) > 120 & coords(:, 2) < 1294;
coords_ref = coords(selected_indices, :);%Remove the recognition coordinates that are not in the crystal region
x=coords_ref(:,2);
y=1350-coords_ref(:,1);
scatter(x, 1350-y, 'filled');
hold on;
distance_H = [];%Store the distance between the atom and its lateral surrounding atoms
distance_V = [];%Store the distance between the atom and its vertical surrounding atoms
for i=1:length(coords_ref(:,1))
    x=coords_ref(i,1);%x
    y=coords_ref(i,2);%y
    %Find the coordinates of the laterally adjacent points
    selected_indices_H = coords_ref(:, 2) > y-37 & coords_ref(:, 2) < y+37 & coords_ref(:, 1) > x-5 & coords_ref(:, 1) < x+5 & coords_ref(:, 2) ~=y;
    coords_H = coords_ref(selected_indices_H,:);
    x1=coords_H(:,2);
    y1=1350-coords_H(:,1);
    scatter(x1, 1350-y1, 'filled');
    length_H = length(coords_H(:,1));
  
    for jH=1:length_H
       dis = pdist([coords_ref(i,:); coords_H(jH,:)]);
       if dis>19
        distance_H=[distance_H,dis];
       end
    end
    %Find the coordinates of the vertical adjacent points
    selected_indices_V = coords_ref(:, 1) > x-37 & coords_ref(:, 1) < x+37 & coords_ref(:, 2) > y-5 & coords_ref(:, 2) < y+5 & coords_ref(:, 2) ~=y;
    coords_V = coords_ref(selected_indices_V,:);
    length_V = length(coords_V(:,1));
    x2=coords_V(:,2);
    y2=1350-coords_V(:,1);
    scatter(x2, 1350-y2, 'filled');
    for jV=1:length_V
       dis = pdist([coords_ref(i,:); coords_V(jV,:)]);
       if dis>19
        distance_V=[distance_V, dis];
       end
    end
end
hold off
distance_H_average = 22.3534;
distance_V_average = 24.2287;
distance_H = [];
distance_V = [];

X_all = [];
Y_all = [];
U_all = [];
V_all = [];
hydrostatic_strain_all=[];
figure,imshow(f,[]);
hold on
result_all = {};
heat=zeros(size(f));
%sorted_points = sortrows(coords, [1, 2]);
for i=1:length(coords(:,1))
    if coords(i,1)<1012

        x=coords(i,1);%x
        y=coords(i,2);%y
        
        %Find the coordinates of the laterally adjacent points
        selected_indices_H = coords(:, 2) > y-30 & coords(:, 2) < y+30 & coords(:, 1) > x-10 & coords(:, 1) < x+10 & coords(:, 2) ~=y;
        coords_H = coords(selected_indices_H,:);
        length_H = length(coords_H(:,1));
        %ind the coordinates of the vertical adjacent points
        selected_indices_V = coords(:, 1) > x-30 & coords(:, 1) < x+30 & coords(:, 2) > y-10 & coords(:, 2) < y+10 & coords(:, 2) ~=y;
        coords_V = coords(selected_indices_V,:);
        length_V = length(coords_V(:,1));
        if length_V ~=2 || length_H ~=2 
            continue
        end
        if abs(coords_V(1,1)-coords_V(2,1))<20 || abs(coords_H(1,2)-coords_H(2,2))<20
            continue
        end
        %Calculate the four vectors
        vectorH1 = coords_H(1,:)-coords(i,:);
        vectorH2 = coords_H(2,:)-coords(i,:);
        vectorV1 = coords_V(1,:)-coords(i,:);
        vectorV2 = coords_V(2,:)-coords(i,:);
        result1 = vectorH1'*vectorH1;
        result2 = vectorH2'*vectorH2;
        result3 = vectorV1'*vectorV1;
        result4 = vectorV2'*vectorV2;
        result_H = result1+result2;
        result_V = result3+result4;
        result = result1+result2+result3+result4;
        trace_sum_H = trace(result_H);
        trace_sum_V = trace(result_V);
        trace_sum = trace(result);
        d0_H=2*(distance_H_average).^2; 
        d0_V=2*(distance_V_average).^2;
        d0 = mean([d0_V,d0_H]);
        hydrostatic_strain_H = (1/2).*((1/d0_H)*trace_sum_H-1);
        hydrostatic_strain_V = (1/2).*((1/d0_V)*trace_sum_V-1);

     
        hydrostatic_strain = (hydrostatic_strain_H+hydrostatic_strain_V)./2;
       
        unitmatix = [[1,0];[0,1]];
       
        Result_all(i).coords=coords(i,:);
        Result_all(i).hydrostatic_strain_H = hydrostatic_strain_H;
      

        x_point=coords(i,2);
        y_point=coords(i,1);
        magnitude=300;
        X_all = [X_all;x_point];
        Y_all = [Y_all;y_point];
        U_all = [U_all;hydrostatic_strain_H];
        V_all = [V_all;hydrostatic_strain_V];
        hydrostatic_strain_all = [hydrostatic_strain_all;hydrostatic_strain];
       
    end
end


hydrostatic_strain_all=hydrostatic_strain_all-0.008;
color_values = hydrostatic_strain_all;


% 
% color_min = min(color_values(:));
% color_max = max(color_values(:));


color_min = -0.08;
color_max = 0.08;


n = 256; 


custom_colormap = turbo;

colorm=custom_colormap;



M=sqrt(U_all.^2+V_all.^2);
Mdown=min(M(:));
Mup=max(M(:));
Mup=0.1;
scaler1=Mup./M;
U=U_all.*scaler1;
V=V_all.*scaler1;

for i = 1:size(X_all, 1)
    
    color_value_normalized = (color_values(i) - color_min) / (color_max - color_min);
    
    if color_value_normalized<=0
        color=colorm(1,:);
    elseif color_value_normalized>1
        color=colorm(256,:);
        
    else
        color = colorm(floor(color_value_normalized *(n-1) )+1 , :);
    end
   
    scatter(X_all(i), Y_all(i), 25,color,'filled');
    
    quiver(X_all(i), Y_all(i), 100*U(i), 100*V(i),1, 'Color', [1,1,1],'LineWidth', 1.5,'MaxHeadSize', 5.0);
    hold on

end
   
    quiver(739.37, 750.59, (-9.81), (-2.05),1, 'Color', 'w','LineWidth', 1.5,'MaxHeadSize', 5.0);

    hold on;

c = colorbar;
c.Colormap = custom_colormap;

c.Ticks =linspace(0, 255, 5);

c.TickLabels = {'-0.08', '-0.04', '0.0','0.04', '0.08'}; 

c.Label.String = 'Colorbar';

% save result
exportgraphics(gca, strcat('D:\github\output\',file,'.png'),'Resolution',1000) 