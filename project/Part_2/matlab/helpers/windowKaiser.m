function G = windowKaiser(nfft, over, dbDown)

    % Generate a tapered window using Kaiser-based beta approximation 
    % for a desired sidelobe level.
    %
    %   G = WindowSlepianKaiser(nfft, over, dbDown)
    %
    %   Inputs:
    %       nfft   - Total window length (FFT size)
    %       over   - Oversampling factor; determines effective taper length.
    %                Effective length L = floor(nfft / over)
    %       dbDown - Desired sidelobe attenuation in dB (positive scalar)
    %
    %   Output:
    %       G      - Window of length nfft normalized such that: sum(G) = nfft

    % Compute non-zero window length L based on oversampling factor
    L = floor(nfft / over);

    % Ensure L is odd (DPSS behaves better for odd symmetric lengths)
    if mod(L, 2) == 0
        L = L - 1;
    end

    A = abs(dbDown);
    if A > 50
        % Higher attenuation → Increase A slightly ("twiddle")
        A = A + 7.5;
        beta = 0.1102 * (A - 8.7);    % Kaiser approximation
    elseif (A >= 21) && (A <= 50)
        % Mid-range attenuation
        A = A + 15;
        beta = 0.5842*(A - 21)^0.4 + 0.07886*(A - 21);
    else
        % Low sidelobe levels → rectangular-like window
        beta = 0;
    end

    % Form Kaiser window
    G = kaiser(L, beta);

    % Zero pad filter
    G = [zeros(ceil((nfft-L)/2),1);G;zeros(floor((nfft-L)/2),1)];

    % Normalize filter
    G = G./sum(G)*length(G);

end