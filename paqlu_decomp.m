function [p, q, L, U] = paqlu_decomp(A)
% Projekt 2, zadanie 43
% Arkadiusz Ułanowski, 320747
%
% Funkcja przeprowadza rozkład macierzy A na macierze dolnotrójkątną
% i górnotrójkątną z pełnym wyborem elementu głównego.
% Mówiąc językiem MATLABa, spełnione jest A(p, q) = L*U.
% Również, aby rozwiązać układ Ax = b, wykonać można następujące:
% y = L\b(p)
% z = U\y
% x(q) = z
%
% Wejście:
%   A - macierz kwadratowa o elementach będących liczbami zespolonymi
% Wyjście:
%   p - wektor permutacji elementów w wierszach macierzy A
%   q - wektor permutacji elementów w kolumnach macierzy A
%   L - macierz dolnotrójkątna z jedynkami na diagonali
%   U - macierz górnotrójkątna

n = length(A);
p = 1:n;
q = p;
L = zeros(n);
U = A;

for i = 1:n-1
    [~, M] = max(abs(U(i:end, i:end)), [], 'all', 'linear');
    [row, column] = ind2sub(size(U(i:end, i:end)), M);
    row = row + i - 1;       % znalezienie współrzędnych największego
    column = column + i - 1; % modułem elementu pod główną przekątną

    U([i, row], :) = U([row, i], :);
    L([i, row], :) = L([row, i], :); % zamiana wierszy,
    p([i, row]) = p([row, i]); % zapisanie zamiany w w. permutacji wierszy
    U(:, [i, column]) = U(:, [column, i]);
    L(:, [i, column]) = L(:, [column, i]); % zamiana kolumn,
    q([i, column]) = q([column, i]); % zapisanie zamiany
                                     % w wektorze permutacji kolumn

    pivots = U(i+1:end, i)./U(i, i);
    L(i+1:end, i) = pivots; % zapisanie współczynników mówiących o tym,
                            % ile razy redukujemy dany wiersz
    U(i+1:end, :) = U(i+1:end, :) - U(i, :).*pivots; % eliminacja Gaussa
end

L(logical(eye(n))) = 1; % wyjedynkowanie diagonali macierzy L