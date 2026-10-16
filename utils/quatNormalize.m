function qn = quatNormalize(q)
%QUATNORMALIZE Normalize a quaternion to unit length (guards against
%numeric drift accumulated by fixed-step integration).

n = norm(q);
if n < eps
    qn = [1; 0; 0; 0];
else
    qn = q / n;
end

end
