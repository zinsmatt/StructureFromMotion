% Dans ce script nous utilisons deux images successives issues d'une
% seule caméra couleur. 
% Il faut rajouter le dossier MatlabFns qui contient les différentes
% fonctions de Peter Kovesi

clear all;
close all;

addpath('MatlabFns\Spatial');
addpath('MatlabFns\Match');
addpath('MatlabFns\Misc');
addpath('MatlabFns\Robust');
addpath('MatlabFns\Projective');
addpath('devkit\matlab');

% chargement de la calibration de la camera
calibration = loadCalibrationCamToCam('calib_cam_to_cam.txt');
K = calibration.K{3};

% Chargement des images
im1 = imread('image_02\0000000040.png');
im2 = imread('image_02\0000000041.png');

% Transformation en niveaux de gris
im1_init = im1;
im2_init = im2;
if size(im1,3) == 3,
    im1 = rgb2gray(im1);
end
if size(im2,3)==3,
    im2 = rgb2gray(im2);
end


% Détection des points d'intérêts
fprintf('Détection des points d intérêt\n');
sigma = 1;
k = 0.04;
[cim1, r1, c1] = harris(im1,sigma,k,'N',10000,'radius',1);
[cim2, r2, c2] = harris(im2,sigma,k,'N',10000,'radius',1);


% Matching des points d'intérêts entre les deux images
fprintf('Matching des points\n');
w = 9;
dmax = 80;  % plus long mais les résultats sont meilleurs
[m1, m2, p1ind, p2ind, cormat] = matchbycorrelation(im1, [r1 c1]', im2, [r2 c2]', w, dmax);
matches = [m1(2,:);m1(1,:);m2(2,:);m2(1,:)];
x1 = matches(1:2,:);
x2 = matches(3:4,:);


% Estimation de la matrice fondamentale et des inliers
[F, inliers] = ransacfitfundmatrix(x1, x2, 0.001);
matches_init = matches;
matches = matches_init(:,inliers);


% Affichage des du matching des inliers
% show((double(im1_init)+double(im2_init))/255,5);hold on;
% plot(matches(1,:),matches(2,:),'.r');hold on;
% plot(matches(3,:),matches(4,:),'.b');
% line([matches(1,:); matches(3,:)],[matches(2,:); matches(4,:)],'color',[0 1 0]);


% Calcul de la matrice essentielle
E = K'*F*K;
[P1,P2,r1,t1,rot_axis,rot_angle,g] = torr_linear_EtoPX(E,matches,K,1);

% Triangulation des inliers matchés
fprintf('Triangulation\n');
matches2 = matches;
X = torr_triangulate(matches',1,P1,P2);

% Passage des coordonnées homogènes en coordonnées cartésiennes
x = []; y = []; z = [];
for i=1:size(X,2),
    if abs(X(4,i))>0.000001,
        x = [x X(1,i)/X(4,i)];
        y = [y X(2,i)/X(4,i)];
        z = [z X(3,i)/X(4,i)];
    else
        x = [x X(1,i)];
        y = [y X(2,i)];
        z = [z X(3,i)];
    end
end

% Suprression des points en z positifs. (par essais on voit que les négatifs doivent être gardés
mat = z>0;
idx = find(mat);
x(idx) = [];
y(idx) = [];
z(idx) = [];
matches2(:,idx) = [];

% Suppression des points trop éloignés du centroïde des points
mx = mean(x);
my = mean(y);
mz = mean(z);
dm = sum(sqrt((x-mx).^2+(y-my).^2+(z-mz).^2),2)/size(x,2);
mat = sqrt((x-mx).^2+(y-my).^2+(z-mz).^2) > dm;
idx = find(mat);
x(idx) = [];
y(idx) = [];
z(idx) = [];
matches2(:,idx) = [];

% Récupération des couleurs des points
coul = [];
for i=1:size(matches2,2),
    v1 = im1_init(matches2(2,i),matches2(1,i),:);
    v2 = im2_init(matches2(2,i),matches2(1,i),:);
    coul = [coul (double(v1) + double(v2))/(2*255)];
end
coul = squeeze(coul);
coul_2 = uint8(coul*255);

% Affichage des points 3D
figure();
scatter3(-x, y,z,10,coul,'filled');

% Création et enregistrement d'un nuage du nuage de points au format .ply
% pour pouvoir le visualiser plus facilement avec MeshLab
ptCloud = pointCloud([-x',y',z'],'Color',coul_2);
pcwrite(ptCloud,'points.ply');


