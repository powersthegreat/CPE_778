function mask = computeCorrelationMask(s, K)

    % Builds an autocorrelation-domain mask whose support corresponds to
    % the reciprocal bandwidth of s, based on the -6 dB PSD width.

    % PSD computation
    S = fft(s(:) ./ sqrt(length(s)), K);
    psd_db = fftshift(10*log10(abs(S).^2 / max(abs(S).^2)));

    % Identify -6 dB bandwidth
    idx = find(psd_db >= -6);
    psd_mask = zeros(K,1);
    psd_mask(idx(1):idx(end)) = 1;

    % Autocorrelation of binary PSD mask
    r = abs(ifft(psd_mask));

    % Find first null (monotonic decrease -> turn point)
    d1 = diff(r);
    firstNull = find(d1 > 0, 1, 'first');  % rising → past minimum
    if isempty(firstNull), firstNull = 1; end

    % Find last null (search in reversed sequence)
    d2 = diff(flipud(r));
    lastNullFromEnd = find(d2 > 0, 1, 'first');
    if isempty(lastNullFromEnd), lastNullFromEnd = 1; end
    lastNull = K - lastNullFromEnd + 1;

    % Build final ACF mask
    mask = zeros(K,1);
    mask(firstNull:lastNull) = 1;

return