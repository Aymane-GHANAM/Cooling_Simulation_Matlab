function [Xs,k]=gauss(A,B,precision)
    %% initialisations
    taille=size(A,1)
    Xa=zeros(taille,1);
    Xb=zeros(taille,1);
    k=1;
    Xbs=zeros(taille,1000);
    %% Calculs
    while(((max(abs(A*Xb-B)))>precision)&&(k<1000))
        for i=1:taille
            somme1=0;
            somme2=0;
            for j=1:i-1
                somme1=somme1+A(i,j)*Xb(j);
            end
            for j=i+1:taille
                somme2=somme2+A(i,j)*Xa(j);
            end
            Xb(i)=(B(i)-somme1-somme2)/A(i,i);
            Xa(i)=Xb(i);
        end
        k=k+1;
    end
    Xs=Xb;
    fprintf("On trouve la matrice Xs en %d itérations",k)
end