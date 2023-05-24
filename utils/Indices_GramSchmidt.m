function [nj,mj] = Indices_GramSchmidt(n_max)
    nj = [];
    mj = [];

    for k = 0:n_max
        for j = 0:k
            n = k-j;
            m = j;
            nj = [nj, n];
            mj = [mj, m];
        end
    end
end

%{
function [nj,mj] = Indices_GramSchmidt(n_max)
    nj = [];
    mj = [];

    for n = 0:n_max
        for m = 0:n_max
            nj = [nj, n];
            mj = [mj, m];
        end
    end
end
%}