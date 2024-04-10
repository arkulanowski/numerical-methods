function coef = get_threeeights_coef(i, j, sizex, sizey)
% Projekt 1, zadanie 42
% Arkadiusz Ułanowski, 320747
%
% Funkcja pomocnicza. Obliczanie kwadratury podwójnej ma postać sumy
% po i, po j z A_ij * f(x_i, y_j). Funkcja zwraca A_ij, czyli współczynnik
% (wagę) przed wartością w podanym węźle, podwójnej kwadratury Newtona 3/8.
% Przykład graficzny: gdy sizex = 10, sizey = 7, macierz współczynników A
% ma rozmiar 10x7 i wygląda następująco:
% 1     3     3     2     3     3     1
% 3     9     9     6     9     9     3
% 3     9     9     6     9     9     3
% 2     6     6     4     6     6     2
% 3     9     9     6     9     9     3
% 3     9     9     6     9     9     3
% 2     6     6     4     6     6     2
% 3     9     9     6     9     9     3
% 3     9     9     6     9     9     3
% 1     3     3     2     3     3     1
% wtedy A_ij = get_threeeights_coef(i, j, 10, 7) dla każdego i, j
%
% Wejście:
%   i -     liczba porządkowa określająca, którą z kolei w ciągu
%           posortowanych rosnąco wszystkich możliwych współrzędnych
%           x węzłów kwadratury jest współrzędna x węzła, dla którego
%           chcemy znaleźć wagę
%   j -     liczba porządkowa określająca, którą z kolei w ciągu
%           posortowanych rosnąco wszystkich możliwych współrzędnych
%           y węzłów kwadratury jest współrzędna y węzła, dla którego
%           chcemy znaleźć wagę        
%   sizex - liczba wszystkich różnych możliwych współrzędnych x
%           występujących w węzłach
%   sizey - liczba wszystkich różnych możliwych współrzędnych y
%           występujących w węzłach
% Wyjście:
%   coef  - waga, z jaką wartość węzła liczy się do sumy przybliżającej
%           całkę

tri_row = mod(i, 3) == 1; % zapamiętujemy, czy element A_ij leży
                          % w co trzecim (1, 4, 7, ...) wierszu
                          % macierzy współczynników
tri_col = mod(j, 3) == 1; % zapamiętujemy, czy element A_ij leży
                          % w co trzeciej (1, 4, 7, ...) kolumnie
                          % macierzy współczynników
vrim = j == 1 || j == sizey; % zapamiętujemy, czy element A_ij leży
                             % w pierwszej lub ostatniej kolumnie
                             % macierzy współczynników
hrim = i == 1 || i == sizex; % zapamiętujemy, czy element A_ij leży
                             % w pierwszym lub ostatnim wierszu
                             % macierzy współczynników

if(vrim || hrim) % element jest na krańcu macierzy
     if(vrim && hrim) % jeśli jest w rogu, ma współczynnik 1
         coef = 1;
         return;
     elseif(vrim && tri_row || hrim && tri_col) % jeśli jest w co trzecim
         coef = 2;                              % wierszu/co trzeciej kol.,
         return;                                % ma współczynnik 2
     else % jeśli nie żadne z powyższych, ma współczynnik 3
         coef = 3;
         return;
     end
else % element nie jest na krańcu macierzy
    if(tri_row || tri_col) % jeśli jest w co trzecim w./co trzeciej kol.
        if(tri_row && tri_col) % jeśli jest jednocześnie
            coef = 4;          % i w co trzecim wierszu i w co trzeciej
            return;            % kolumnie, ma współczynnik 4
        else % jeśli jest w dokładnie jednym z: {co trzecim wierszu,
            coef = 6; % co trzeciej kolumnie}, ma współczynnik 6
            return;
        end
    else % jeśli nie żadne z powyższych, ma współczynnik 9
        coef = 9;
        return;
    end
end