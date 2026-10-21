function check_euler_step(D, which)
% CHECK_EULER_STEP  Exact prerequisites for the explicit or implicit Euler bound.
% Performed before a numerical search; the independent replays repeat their own checks.
if ~isscalar(D.hq) || ~isAlways(D.hq > 0)
    error('check_euler_step:step', 'h must be positive');
end
for v = 1:numel(D.ver)
    for i = 1:D.m
        A = D.ver(v).A{i};
        if any(double(sign(A(~eye(D.d)))) < 0)
            error('check_euler_step:metzler', 'A_i must be Metzler');
        end
        if strcmp(which, 'explicit')
            if any(double(sign(D.hq * abs(diag(A)) - 1)) >= 0)
                error('check_euler_step:explicit', 'explicit Euler requires h max|a_jj| < 1 in each reading');
            end
        elseif strcmp(which, 'implicit')
            M = sym(eye(D.d)) - D.hq * A;
            if isAlways(det(M) == 0)
                error('check_euler_step:implicit', 'I - h A_i is singular');
            end
            z = M \ sym(ones(D.d, 1));
            if ~all(isAlways(z > 0)) || ~all(isAlways(M * z > 0))
                error('check_euler_step:implicit', 'implicit Euler requires h eta(A_i) < 1 in each reading');
            end
        else
            error('check_euler_step:kind', 'unknown Euler family');
        end
    end
end
end
