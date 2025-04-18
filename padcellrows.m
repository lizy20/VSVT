function M = padcellrows(C)
    % C: a cell array where each row is a numeric vector (e.g., element node list)
    % Returns a numeric matrix padded with NaNs for unequal lengths

    maxlen = max(cellfun(@numel, C));      % Find the maximum row length
    M = nan(length(C), maxlen);            % Preallocate output with NaNs
    for i = 1:length(C)
        M(i, 1:numel(C{i})) = C{i};        % Fill each row with values
    end
end