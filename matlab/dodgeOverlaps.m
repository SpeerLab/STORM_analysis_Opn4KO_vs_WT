function offs = dodgeOverlaps(y, minGap, step)
% DODGEOVERLAPS  Horizontal offsets that fan out points colliding in y.
%
%   offs = dodgeOverlaps(y, minGap, step)
%
% Points whose y-values fall within MINGAP of a neighbour are treated as a
% cluster and spread symmetrically about x = 0 with spacing STEP
% (..., -step, 0, +step, ... for odd counts; +/-step/2, +/-3*step/2, ... for
% even counts). Isolated points stay centred at 0. Used to keep the N = 6
% per-animal dots from overlapping on the centre line.

    y = y(:);
    n = numel(y);
    offs = zeros(n, 1);
    if n < 2, return; end

    [ys, order] = sort(y);

    % group consecutive points that are within minGap of the previous one
    groupId = zeros(n, 1);
    gid = 0;
    for i = 2:n
        if ys(i) - ys(i-1) < minGap
            groupId(i) = gid;
        else
            gid = gid + 1;
            groupId(i) = gid;
        end
    end

    for g = unique(groupId).'
        idx = find(groupId == g);
        m = numel(idx);
        if m == 1, continue; end
        pos = ((0:m-1) - (m-1)/2) * step;   % symmetric about 0
        offs(order(idx)) = pos;
    end
end
