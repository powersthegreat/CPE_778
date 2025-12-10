function plotConvergence(hist)

    K = length(hist.obj);

    % === Normalization (only if the initial values are nonzero) ===
    obj_norm   = hist.obj   / max(abs(hist.obj(1)),   1e-12);
    hnorm_norm = hist.hnorm / max(abs(hist.hnorm(1)), 1e-12);
    gval_norm  = hist.gval  / max(abs(hist.gval(1)),  1e-12);

    % Gradient norm normalization
    grad_norm  = hist.gradNorm;
    grad_normN = grad_norm / max(abs(grad_norm(1)), 1e-12);

    % Penalty parameters (no normalization needed)
    rho_vals = hist.rho;
    eta_vals = hist.eta;

    % === FIGURE ===
    figure('Name','ALM Convergence Overview','Color','w');
    tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

    % -------------------------------------------------------------
    % Objective
    % -------------------------------------------------------------
    nexttile;
    semilogy(obj_norm,'LineWidth',1.5);
    grid on;
    ylabel('Norm. Objective');
    title('Objective Convergence');

    % -------------------------------------------------------------
    % Equality Constraint Norm
    % -------------------------------------------------------------
    nexttile;
    semilogy(hnorm_norm,'LineWidth',1.5);
    grid on;
    ylabel('Norm. ||h||');
    title('Equality Constraint Convergence');

    % -------------------------------------------------------------
    % Inequality Constraint Value
    % -------------------------------------------------------------
    nexttile;
    semilogy(gval_norm,'LineWidth',1.5);
    grid on;
    ylabel('Norm. g(x)');
    title('Inequality Constraint Convergence');

    % -------------------------------------------------------------
    % Gradient Norm
    % -------------------------------------------------------------
    nexttile;
    semilogy(grad_normN,'LineWidth',1.5);
    grid on;
    ylabel('Norm. ||∇L||');
    xlabel('Outer Iteration');
    title('Gradient Norm Convergence');

    % -------------------------------------------------------------
    % Penalty ρ
    % -------------------------------------------------------------
    nexttile;
    plot(rho_vals,'LineWidth',1.5);
    grid on;
    ylabel('\rho_k');
    xlabel('Outer Iteration');
    title('Equality Penalty Parameter');

    % -------------------------------------------------------------
    % Penalty η
    % -------------------------------------------------------------
    nexttile;
    plot(eta_vals,'LineWidth',1.5);
    grid on;
    ylabel('\eta_k');
    xlabel('Outer Iteration');
    title('Inequality Penalty Parameter');

end
