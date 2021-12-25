clear all 
close all


%% variables de temps
dt=0.1;
T=200;
n=T/dt;

%% initialisation de A1 et A2


A1=-4*eye(400,400); %la matrice pour la propagation des températures pour l'état stable
A2=-4*eye(400,400); %la matrice pour la propagation des températures pour l'état évolutif

%remplissage de A1
for i = 1:400
    if (i<400) 
        A1(i,i+1)=1;
    end
    if (i>1) 
        A1(i,i-1)=1;
    end
    if((i>20)&&(i<=380))
        A1(i,i-20)=1;
        A1(i,i+20)=1;
    end
    if((i>=2)&&(i<=10))
        A1(i,i+20)=1;
        A1(i,i+300)=1;
    end
    if((i>=302)&&(i<=310))
        A1(i,i)=-5;
        A1(i,i-300)=1;
    end
    %les bords et les points à 300° ont une valeur de température fixe
    if ((mod((i-1),20)==0)||(mod(i,20)==0)||(i>380)||((mod((i-10),20)==0)&&(i<=210))||((mod((i-11),20)==0)&&(i<=211))||(i==27)||(i==28)||(i==76)||(i==97)||((i<=20)&&(i>10)))
        A1(i,:)=0;
        A1(i,i)=1;
    end
end

%remplissage de A2
for i = 1:400
    if (i<400) 
        A2(i,i+1)=1;
    end
    if (i>1) 
        A2(i,i-1)=1;
    end
    if((i>20)&&(i<=380))
        A2(i,i-20)=1;
        A2(i,i+20)=1;
    end
    if((i>=2)&&(i<=10))
        A2(i,i+20)=1;
        A2(i,i+300)=1;
    end
    if((i>=302)&&(i<=310))
        A2(i,i)=-5;
        A2(i,i-300)=1;
    end
    if (((mod((i-1),20)==0)&&(i~=1)&&(i~=381))||((mod(i-11,20)==0)&&(i~=11)&&(i<211)))
        A2(i,i)=-3;
        A2(i,i-20)=1;
        A2(i,i+20)=1;
        A2(i,i+1)=1;
        A2(i,i-1)=0;
    end %bords gauches
    if (((mod(i,20)==0)&&(i~=20)&&(i~=400))||((mod(i-10,20)==0)&&(i~=10)&&(i<210)))
        A2(i,i)=-3;
        A2(i,i-20)=1;
        A2(i,i+20)=1;
        A2(i,i-1)=1;
        A2(i,i+1)=0;
    end %bords droits
    
    if ((i>=382)&&(i~=400))
        A2(i,i)=-3;
        A2(i,i-20)=1;
        A2(i,i+1)=1;
        A2(i,i-1)=1;
    end %dessous
    
    if ((i>11)&&(i<20))
        A2(i,i)=-3;
        A2(i,i+20)=1;
        A2(i,i+1)=1;
        A2(i,i-1)=1;
    end %dessus droit
    
    if (i==11)
        A2(i,i)=-2;
        A2(i,i+1)=1;
        A2(i,i+20)=1;
        A2(i,i-1)=0;
    end %coin 11
    
    if (i==20)
        A2(i,i)=-2;
        A2(i,i-1)=1;
        A2(i,i+20)=1;
        A2(i,i+1)=0;
    end %coin 20
    
    if (i==381)
        A2(i,i)=-2;
        A2(i,i+1)=1;
        A2(i,i-20)=1;
        A2(i,i-1)=0;
    end %coin 381
    
    if (i==400)
        A2(i,i)=-2;
        A2(i,i-1)=1;
        A2(i,i-20)=1;
    end %coin 400
    
    if (i==1)
        A2(i,i)=-3;
        A2(i,i+20)=1;
        A2(i,i+1)=1;
        A2(i,i+300)=1;
    end %coin 1
    
    if (i==10)
        A2(i,i)=-3;
        A2(i,i+20)=1;
        A2(i,i-1)=1;
        A2(i,i+300)=1;
    end %coin 10
    
    
        
    
end

%% initialisation de X et B

% B est un vecteur composé des valeurs de X à t=0
B = zeros(400, 1);

% Initialisation du vecteur de température X
X = zeros(400, 1);

% Points sur le bord à 50°C
for i=1:20:381
    B(i) = 50;
end

for i=20:20:400
    B(i) = 50;
end

for i=381:400
    B(i) = 50;
end

for i=10:20
    B(i) = 50;
end

for i=30:20:210
    B(i) = 50;
end

for i=31:20:211
    B(i) = 50;
end

% Points à 300°C
B(27) = 300;
B(28) = 300;
B(76) = 300;
B(97) = 300;

%% Calcul des températures à l'état stable
[X,L]=gauss(A1,B,0.01);

Xaff1=zeros(20,20); %initialisation du tableau de températures pour l'affichage
for i=1:20
    for j=1:20
        Xaff1(i,j)=X(i+(j-1)*20);
    end
end
figure(1)

hold on
surf(Xaff1)
title('Carte des températures, à température stable')
hold off


%% Calcul des températures évoluantes
Xev=zeros(400,n);   %Variable stockant les températures pour chaque moment
Xev(:,1)=X;
expA=expm(A2*dt); %pour empêcher de recalculer cette valeur à chaque instance de la boucle suivante
for h=2:n
    Xev(:,h)=expA*Xev(:,h-1); 
    h %on affiche l'avancement
end
Xaff2=zeros(20,20,n);  %les valeurs de Xev pour l'affichage, les températures ne sont plus en colonne 400*1 mais en matrice 20*20  
for i=1:20
    for j=1:20
        Xaff2(i,j,:)=Xev(i+(j-1)*20,:);
    end
end

%% Tracé des températures évoluantes
figure(2)
for k2 =1:n                               % Tracé évolutif
        surf(Xaff2(:,:,k2))                  
        axis([0 20    0 20    0 300])          
        drawnow                                    
        pause(dt)
        dt*k2                             % Affichage du temps passé à chaque nouveau tracé
end