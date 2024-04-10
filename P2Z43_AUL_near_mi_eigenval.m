function [itNum, precision, eigenval] = P2Z43_AUL_near_mi_eigenval(A, mi, itMax, desiredPrecision, TESTING)
% Projekt 2, zadanie 43
% Arkadiusz Ułanowski, 320747
%
% Obliczanie wartości własnej macierzy A leżącej najbliżej podanej wartości
% mi poprzez obliczanie odwrotną metodą potęgową z normowaniem
% najmniejszej co do modułu wartości własnej B = A - mi*I.
% Występujące w każdej iteracji metody dzielenie przez B zastąpiono
% rozwiązywaniem układu równań liniowych z wykorzystaniem
% uprzedniej dekompozycji PBQ = LU.
%
% Wejście:
%   A                - macierz kwadratowa o elementach będących
%                      liczbami zespolonymi
%   mi               - liczba zespolona, dla której szukamy
%                      najbliższej jej wartości własnej macierzy A
%   itMax            - maksymalna dopuszczalna przez użytkownika liczba
%                      wykonanych iteracji metody (brak wartości domyślnej,
%                      użytkownik podaje itMax jawnie, aby uniknąć
%                      nieporozumień w interpretacji wyjściowego itNum)
%   desiredPrecision - żądana w warunku stopu dokładność wyniku
%                      (domyślnie 1e-7)
%   TESTING          - flaga wymuszająca wzięcie za przybliżenie początkowe
%                      wektora samych jedynek, domyślnie fałsz,
%                      przeznaczona wyłącznie do testów numerycznych
%                      funkcji
% Wyjście:
%   itNum            - liczba wykonanych iteracji metody (w szczególnym
%                      przypadku, gdy itNum = itMax + 1,
%                      metoda wykonała się itMax razy i nie znalazła
%                      rozwiązania na zadaną przez użytkownika dokładność)
%   precision        - dokładność zwróconego wyniku (użytkownik może
%                      samodzielnie sprawdzić, czy jest zadowalająca)
%   eigenval         - wyliczone przybliżenie wartości własnej A leżącej
%                      najbliżej wartości mi

arguments
    A
    mi
    itMax
    desiredPrecision = 1e-7
    TESTING = false;
end

itNum = 0; % przypisanie wartości początkowych
precision = 0;
n = length(A);
B = A - mi*eye(n);
if(rcond(B) <= eps("double")) % jeśli macierz B jest niepełnego wymiaru,
    eigenval = mi; % mi jest wartością własną
    precision = eps("double"); % z dokładnością do epsilona maszynowego
    return;
end
[p, q, L, U] = paqlu_decomp(B); % rozkład macierzy B

if(TESTING) % jeśli wywołano w trybie testu funkcji
    v = ones(n, 1); % początkowe przybliżenie wektora własnego jest
                    % wektorem samych jedynek
else
    v = rand(n, 1); % początkowe przybliżenie wektora własnego jest losowe
end
v = v/norm(v); % unormowanie początkowego przybliżenia
found = false;
LT.LT = true; % utworzenie struktur przekazywanych dalej funkcji linsolve,
UT.UT = true; % celem bezpośredniego wskazania, że układy równań są
              % dolno- i górnotrójkątne

while(itNum < itMax)
    itNum = itNum + 1;
    v_prev = v;
    y = linsolve(L, v(p), LT); % rozwiązanie układów równań z macierzą
    z = linsolve(U, y, UT); % dolno- i górnotrójkątną
    v(q) = z;
    v = v/norm(v); % unormowanie rozwiązania
    
    [~, ind] = max(abs(v_prev)); % znalezienie indeksu największego
                              % modułem elementu poprzedniego przybliżenia
    precision = norm((1/v(ind))*abs(v(ind))*v... % wyliczenie normy
    - (1/v_prev(ind))*abs(v_prev(ind))*v_prev);  % różnicy przeskalowanych
                                                 % kolejnych przybliżeń
    if(precision <= desiredPrecision) % jeśli osiągnięto precyzję:
        found = true; % znaleziono rozwiązanie, zakończ obliczenia
        break;
    end
end

if(~found) % jeśli nie znaleziono rozwiązania
    itNum = itMax + 1; % liczba iteracji jest równa maksymalnej + 1
end

eigenval = v'*A*v; % wylicz w. własną z ułamka Rayleigha