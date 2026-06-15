
    %% ---------- VSTUPNÍ MATICE (uprav dle potřeby) ----------
    % Příklad: 2D konstantní rychlost (constant velocity model)
Ts = 0.1;
    
    
    A=kron([1 Ts; 0 1], eye(2));
    B=kron([Ts^2/2; Ts], eye(2));
    C=[eye(2) zeros(2)];
    
    q = 1e-3;
    Q = q*kron([Ts^3/3 Ts^2/2;Ts^2/2 Ts], eye(2))

       r = 5;
    R_st = r * eye(2);
    nu = 5;
    %normalni
    R = R_st * nu /(nu-2);

    
    P0 = eye(4);              % prior kovariance P_{0|-1} (nebo P_{0|0})

    %% ---------- PARAMETRY EXPERIMENTU ----------
    N     = 1000;   % fixní čas, jehož odhad sledujeme
    K_max = 2000;   % poslední čas, do kterého přidáváme měření

    n = size(A,1);

    %% ---------- HLAVNÍ SMYČKA PŘES HORIZONT k ----------
    k_values   = N:K_max;                 % 100, 101, ..., 200
    trace_PNk  = zeros(size(k_values));   % výsledky tr(P_{N|k})

    for idx = 1:numel(k_values)
        k = k_values(idx);

        % --- 1) Dopředná filtrace (kovariance only) od 0 do k ---
        % Uložíme prediktivní P_{t|t-1} a filtrované P_{t|t} pro RTS pass.
        Pp = cell(k+1,1);   % P_{t|t-1}  (prediktivní)
        Pf = cell(k+1,1);   % P_{t|t}    (filtrované)

        % Inicializace v čase t = 0
        P_pred = P0;        % P_{0|-1}
        for t = 0:k
            % --- Update (měření v čase t) ---
            S = C*P_pred*C' + R;          % inovační kovariance
            K_gain = P_pred*C'/S;          % Kalmanův zisk
            P_filt = (eye(n) - K_gain*C)*P_pred;
            % symetrizace (numerická stabilita)
            P_filt = 0.5*(P_filt + P_filt');

            Pp{t+1} = P_pred;
            Pf{t+1} = P_filt;

            % --- Predikce do času t+1 ---
            P_pred = A*P_filt*A' + Q;
            P_pred = 0.5*(P_pred + P_pred');
        end

        % --- 2) RTS backward pass z času k zpět do času N ---
        % Smoothovaná kovariance P_{t|k}. Začneme P_{k|k} = Pf{k+1}.
        P_smooth = Pf{k+1};               % P_{k|k}
        for t = (k-1):-1:N
            % G_t = P_{t|t} A' (P_{t+1|t})^{-1}
            G = Pf{t+1} * A' / Pp{t+2};
            P_smooth = Pf{t+1} + G*(P_smooth - Pp{t+2})*G';
            P_smooth = 0.5*(P_smooth + P_smooth');
        end
        % Po skončení smyčky je P_smooth == P_{N|k}

        trace_PNk(idx) = trace(P_smooth);
    end

    %% ---------- VYKRESLENÍ ----------
    figure('Color','w'); hold on; grid on; box on;
    plot(k_values, trace_PNk, '-o', 'LineWidth', 1.6, ...
         'MarkerSize', 4, 'MarkerFaceColor', [0 0.45 0.74]);

    % Zvýraznění bodu filtrace (k = N)
    plot(N, trace_PNk(1), 'rs', 'MarkerSize', 10, ...
         'LineWidth', 2, 'MarkerFaceColor', 'r');

    xlabel('horizont měření $k$', 'Interpreter','latex');
    ylabel('$\mathrm{tr}(P_{N|k}) = \mathbb{E}\,\|\hat{x}_{N|k} - x_N\|_2^2$', ...
           'Interpreter','latex');
    title(sprintf('Přínos vyhlazování: chyba odhadu v čase N=%d', N), ...
          'Interpreter','latex');
    legend({'$\mathrm{tr}(P_{N|k})$ (smoothing)', ...
            'filtrace $k=N$'}, ...
           'Interpreter','latex', 'Location','northeast');

    %% ---------- TEXTOVÝ VÝPIS ----------
    fprintf('Filtrace   tr(P_{%d|%d})  = %.4f\n', N, N, trace_PNk(1));
    fprintf('Smoothing  tr(P_{%d|%d})  = %.4f\n', N, K_max, trace_PNk(end));
    fprintf('Redukce chyby vyhlazováním: %.1f %%\n', ...
        100*(1 - trace_PNk(end)/trace_PNk(1)));
    