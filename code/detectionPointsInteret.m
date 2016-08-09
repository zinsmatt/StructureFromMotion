function [ row, col ] = detectionPointsInteret( img, N, w, rayon , methode)
% Renvoie les coordonnées de points d'intérêt
% Paramètre :
% - img : image à traiter
% - N : nombre maximal de points détectés
% - w : taille de la fenêtre d'analyse
% - rayon : rayon considéré lors de la suppression des non maxima
% - méthode : définit le masque utilisé pour calculer les dérivées directionnelles
%             choix : Prewitt, Sobel ou Gradient

% Transformation de l'image de niveaux de gris si nécessaire
if size(img,3) == 3,
    img = rgb2gray(img);
end
img = double(img);

w2 = (w-1)/2;

% Choix du masque pour les dérivées directionnelles
if strcmp(methode,'prewitt') == 1
    Hx = [-1,0,1;-1,0,1;-1,0,1];
elseif strcmp(methode,'sobel') == 1,
    Hx = [-1,0,1;-2,0,2;-1,0,1];
elseif strcmp(methode,'gradient') == 1,
    Hx = [0,0,0;-1,0,1;0,0,0];
end  
Hy = Hx';

% Calcul des dérivées directionnelles
dimgx = filter2(Hx,img);
dimgy = filter2(Hy,img);

% Calcul du masque Gaussien centré
xc = (w+1)/2;
yc = (w+1)/2;
sigma = 1;
masque_gauss = zeros(w,w);
for j=1:w,
    for i=1:w,
        masque_gauss(j,i) = -((i-xc)^2+(j-yc)^2)/(2*sigma^2);
    end 
end
masque_gauss = (1/(2*pi*sigma^2))*exp(masque_gauss);
    
R = zeros(size(img));
k = 0.04;

% Calcul de la réponse en chaque pixel
for j = 1+w2:size(img,1)-w2,
    for i = 1+w2:size(img,2)-w2,
        
        Ix = dimgx(max(1,j-w2):min(size(img,1),j+w2),max(1,i-w2):min(size(img,2),i+w2));
        Iy = dimgy(max(1,j-w2):min(size(img,1),j+w2),max(1,i-w2):min(size(img,2),i+w2));

        % Calcul de la matrice d'autocorrélation 
        ta = masque_gauss.*(Ix.^2);
        tb = masque_gauss.*(Iy.^2);
        tc = masque_gauss.*(Ix.*Iy);

        A = sum(ta(:));
        B = sum(tb(:));
        C = sum(tc(:));

        % det(M)-k*trace(M)^2
        R(j,i) = (A*B-C^2)-k*(A+B)^2;
    end
end


% Suppression des non maxima
[sortedX,sortingIndices] = sort(R(:),'descend');
[j,i] = ind2sub(size(R),sortingIndices);
zones = zeros(size(img));
save_i = [];
save_j = [];
k=1;
while k<=size(j,1) && size(save_i,2)<=N,
   vi = i(k);
   vj = j(k);
   non_maxi = 0;
   if zones(vj,vi) == 0,
      zones(max(1,vj-rayon):min(size(img,1),vj+rayon),max(1,vi-rayon):min(size(img,2),vi+rayon)) = 1;
   else
       non_maxi = 1;
   end
   if non_maxi == 0,
       save_i = [save_i vi];
       save_j = [save_j vj];
   end
    k=k+1;
end
row = save_j';
col = save_i';
end