function I = P1Z42_AUL_threeeighths(f, m, n)
% Projekt 1, zadanie 42
% Arkadiusz Ułanowski, 320747
%
% Obliczanie całek podwójnych po obszarze D = {(x,y) € R^2, |x| + |y| <= 1}
% z zastąpieniem funkcji całkowanej f funkcją g, której
% przeskalowana całka po K = [-1; 1] x [-1; 1] jest równa całce z f po D,
% a następnie przybliżaniem całki g po K rozbijając na całki iterowane
% i stosując złożoną kwadraturę Newtona "3/8" ze względu na obydwie zmienne
%
% Wejście:
%   f -     uchwyt do funkcji dwóch zmiennych rzeczywistych całkowalnej
%           w sensie Riemanna na D
%   m -     niezerowa liczba naturalna równych części, na które podzielono
%           [-1; 1] przy przybliżaniu całki iterowanej po x (domyślnie 50,
%           kierunek OX w kwadracie K wyznacza w obszarze D prosta y = x)        
%   n -     niezerowa liczba naturalna równych części, na które podzielono
%           [-1; 1] przy przybliżaniu całki iterowanej po y (domyślnie 50,
%           kierunek OY w kwadracie K wyznacza w obszarze D prosta y = -x)
% Wyjście:
%   I -     otrzymane przybliżenie całki z f po D

arguments
    f function_handle
    m = 50
    n = 50
end

g = @(x, y) f((x + y)/2, (x - y)/2); % zmiana obszaru całkowania
H_m = 2/m; % długość jednego przedziału po x
H_n = 2/n; % długość jednego przedziału po y;

x = -1 : H_m/3 : 1; % współrzędne x węzłów
y = -1: H_n/3 : 1; % współrzędne y węzłów
sizex = size(x, 2);
sizey = size(y, 2);
I = 0;
comf = H_m*H_n/128; % czynnik skalujący zależny od kroków całkowania
                    % oraz jakobianu przekształcenia f -> g
for i = 1:sizex
    for j = 1:sizey
        I = I + g(x(i), y(j))...
            *get_threeeights_coef(i, j, sizex, sizey);
    end
end 
I = comf*I; % przeskalowanie i otrzymanie kwadratury