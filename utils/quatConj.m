function qc = quatConj(q)
%QUATCONJ Conjugate of quaternion [w;x;y;z].

qc = [q(1); -q(2); -q(3); -q(4)];

end
