clear; clc; close all;

syms x real 

Px = 4.5;
Py = 8;
epsilon = 1e-6;
f = 7 - x^(5/2);
D = (x - Px)^2 + (f - Py)^2;


fprintf('1. Funkcja celu (Sformułowana symbolicznie):\n');
pretty(D) 

d_sym = diff(D, x);

fprintf('\n2. Pochodna funkcji celu D''(x):\n');
pretty(simplify(d_sym))

rozwiazanie_vpasolve = vpasolve(d_sym == 0, x, [0 Inf]);
x_opt_vpasolve = double(rozwiazanie_vpasolve);

% Obliczenie wartości y oraz minimalnej odległości
y_opt_vpasolve = double(subs(f, x, x_opt_vpasolve));
min_dist_sq_vpasolve = double(subs(D, x, x_opt_vpasolve));
min_dist_vpasolve = sqrt(min_dist_sq_vpasolve);

fprintf('Znalezione x (vpasolve): %.6f\n', x_opt_vpasolve);
fprintf('Znalezione y (vpasolve): %.6f\n', y_opt_vpasolve);
fprintf('Minimalna odległość: %.6f\n\n', min_dist_vpasolve);

%% 2. Metoda numeryczna
f_num = @(x) 7 - x.^(5/2);
D_num = @(x) (x - Px).^2 + (f_num(x) - Py).^2;

options = optimset('Display','iter','TolX',epsilon);

fprintf('--- ROZWIĄZANIE NUMERYCZNE (fminbnd) ---\n');
[x_opt_num, fval_num, exitflag, output] = fminbnd(D_num, 0, 5, options);

y_opt_num = f_num(x_opt_num);
min_dist_num = sqrt(fval_num);

fprintf('\nZnalezione x (numerycznie): %.6f\n', x_opt_num);
fprintf('Znalezione y (numerycznie): %.6f\n', y_opt_num);
fprintf('Minimalna odległość: %.6f\n', min_dist_num);
fprintf('Liczba iteracji: %d\n', output.iterations);

%% 3. INTERPRETACJA GRAFICZNA
figure('Name', 'Projekt Optymalizacja', 'NumberTitle', 'off', 'Color', 'w');

% Wykres 1: Geometria problemu (Krzywa i najkrótszy odcinek)
subplot(1,2,1);
hold on; grid on; axis equal;
x_plot = 0:0.01:5;
y_plot = f_num(x_plot);

plot(x_plot, y_plot, 'b-', 'LineWidth', 1.5, 'DisplayName', 'f(x) = 7 - x^{5/2}');
plot(Px, Py, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8, 'DisplayName', 'Punkt P(4.5, 8)');
plot(x_opt_num, y_opt_num, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8, 'DisplayName', 'Punkt Optymalny');
plot([Px x_opt_num], [Py y_opt_num], 'k--', 'LineWidth', 1.5, 'DisplayName', 'Najkrótsza odległość');

xlabel('x'); ylabel('y');
title('Interpretacja Geometryczna');
legend('Location', 'southwest');
ylim([-10 10]);

% Wykres 2: Funkcja celu (Parabola minimalizacji)
subplot(1,2,2);
g_plot = D_num(x_plot);
plot(x_plot, g_plot, 'm-', 'LineWidth', 1.5);
hold on; grid on;
plot(x_opt_num, fval_num, 'ko', 'MarkerFaceColor', 'k');

xlabel('x'); ylabel('Kwadrat odległości g(x)');
title('Wykres Funkcji Celu g(x)');
text(x_opt_num, fval_num+5, sprintf('Min: x=%.4f', x_opt_num), 'HorizontalAlignment', 'center');

sgtitle('Minimalizacja odległości punktu od funkcji');

%% --- CZĘŚĆ 4: METODA BISEKCJI ---
fprintf('\n--- METODA BISEKCJI ---\n');
a = 0;
b = 5;
iter = 0;

fileID = fopen('wyniki_bisekcja.txt', 'w');
fprintf(fileID, 'Iteracja\t x_srodek\t Blad_bezwzgledny\n');

fprintf('Iter\t x_approx\t Blad\n');

while (b - a) >= (2 * epsilon)
    iter = iter + 1;
    
    x1 = a + (b - a)/4;
    xm = (a + b)/2;
    x2 = b - (b - a)/4;
    
    fx1 = D_num(x1);
    fxm = D_num(xm);
    fx2 = D_num(x2);

    if fx1 < fxm
        b = xm;
    elseif fx2 < fxm
        a = xm;
    else
        a = x1;
        b = x2;
    end
    x_approx = (a + b)/2;
    
    blad = abs(x_opt_vpasolve - x_approx);
    
    % Zapis do pliku i wyświetlenie w konsoli
    fprintf(fileID, '%d\t %.8f\t %.8f\n', iter, x_approx, blad);
    fprintf('%d\t %.6f\t %.8f\n', iter, x_approx, blad);
end

x_bisekcja = (a + b)/2;
y_bisekcja = f_num(x_bisekcja);
dist_bisekcja = sqrt(D_num(x_bisekcja));

fprintf('\n--- PODSUMOWANIE BISEKCJI ---\n');
fprintf('Wynik zapisano w pliku: wyniki_bisekcja.txt\n');
fprintf('Znalezione x (bisekcja): %.6f\n', x_bisekcja);
fprintf('Liczba iteracji: %d\n', iter);

%% ZESTAWIENIE REZULTATÓW  
fprintf('\n==============================================\n');
fprintf('     ZESTAWIENIE REZULTATÓW KOŃCOWYCH\n');
fprintf('==============================================\n');
fprintf('%-15s | %-12s | %-12s\n', 'Metoda', 'x*', 'Min. Odleglosc');
fprintf('----------------------------------------------\n');
fprintf('%-15s | %.6f     | %.6f\n', 'Matlab vpasolve', x_opt_vpasolve, min_dist_vpasolve);
fprintf('%-15s | %.6f     | %.6f\n', 'Matlab fminbnd', x_opt_num, min_dist_num);
fprintf('%-15s | %.6f     | %.6f\n', 'Własna Bisekcja', x_bisekcja, dist_bisekcja);
fprintf('----------------------------------------------\n');

max_iter = 100;
%% Walec
fprintf('Zadanie: Minimalizacja powierzchni walca przy stałej objętości V = %.6f\n', x_opt_vpasolve);
TargetFunc = @(r) 2*pi*r.^2 + 2*x_opt_vpasolve./r;
dTarget = @(r) 4*pi*r - 2*x_opt_vpasolve./(r.^2);
ddTarget = @(r) 4*pi + 4*x_opt_vpasolve./(r.^3);

fileID = fopen('wyniki_newton.txt', 'w');
fprintf(fileID, 'Iteracja\t Promien_r\t Pole_P(r)\t Pierwsza_Pochodna\n');
fprintf('Iter\t r_approx\t P(r)\t\t P''(r)\n');
fprintf('-----------------------------------------------------\n');

r_current = 1.0;
iter = 0;

while true
    iter = iter + 1;
    
    f_val = TargetFunc(r_current);
    df_val = dTarget(r_current);
    ddf_val = ddTarget(r_current);
    
    fprintf(fileID, '%d\t %.8f\t %.8f\t %.8f\n', iter, r_current, f_val, df_val);
    fprintf('%d\t %.6f\t %.6f\t %.8f\n', iter, r_current, f_val, df_val);
    
    r_next = r_current - (df_val / ddf_val);
    
    if abs(r_next - r_current) < epsilon || iter >= max_iter
        r_opt = r_next;
        break;
    end
    
    r_current = r_next;
end

fclose(fileID);

h_opt = x_opt_vpasolve / (pi * r_opt^2);
S_opt = TargetFunc(r_opt);

fprintf('-----------------------------------------------------\n');
fprintf('\nWYNIKI KOŃCOWE:\n');
fprintf('Plik z historią iteracji: wyniki_newton.txt\n');
fprintf('Optymalny promień (r):    %.6f [m]\n', r_opt);
fprintf('Optymalna wysokość (h):   %.6f [m]\n', h_opt);
fprintf('Minimalna powierzchnia:   %.6f [m^2]\n', S_opt);
fprintf('Liczba iteracji:          %d\n', iter);

fprintf('Stosunek h/r:             %.4f\n', h_opt/r_opt);

