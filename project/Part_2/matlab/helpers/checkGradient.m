function checkGradient(f, df, x0)

    eps = 1e-7;
    x0  = x0(:);

    % Analytical gradient
    analytic = df(x0);

    % Preallocate
    numeric = zeros(length(x0), 1);
    error   = zeros(length(x0), 1);

    % === Finite-Difference Loop ===
    for k = 1:length(x0)
        e = zeros(length(x0),1);
        e(k) = 1;

        numeric(k) = (f(x0 + eps*e) - f(x0 - eps*e)) / (2*eps);
        error(k)   = abs(numeric(k) - analytic(k));
    end

    % === Print summary statistics ===
    fprintf('\n--- NUMERICAL GRADIENT SUMMARY ---\n');
    fprintf(' Max error : %.3e\n', max(error));
    fprintf(' Mean error: %.3e\n', mean(error));
    fprintf(' Min error : %.3e\n', min(error));
    fprintf('--------------------------------\n\n');

    % === Plot results ===
    figure('Name','Gradient Check','Color','w');
    subplot(3,1,1);
    plot(real(numeric), 'bo-', 'LineWidth', 1.2); hold on;
    plot(real(analytic), 'r.--', 'LineWidth', 1.2);
    ylabel('Gradient');
    title('Real Part: Numerical vs Analytical');
    legend('Numerical','Analytical');
    grid on;

    subplot(3,1,2);
    plot(imag(numeric), 'bo-', 'LineWidth', 1.2); hold on;
    plot(imag(analytic), 'r.--', 'LineWidth', 1.2);
    ylabel('Gradient');
    title('Imag Part: Numerical vs Analytical');
    legend('Numerical','Analytical');
    grid on;

    subplot(3,1,3);
    semilogy(error, 'k', 'LineWidth', 1.2);
    xlabel('Index');
    ylabel('Absolute Error');
    title('Absolute Gradient Error (log scale)');
    grid on;

end
