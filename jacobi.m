function A = jacobi(A,idx,para)
    if para == 1 % one-side for SVD
        % [~,~,V] = svd(A(:,idx));
        AA = A(:,idx)'*A(:,idx);
        tau = (AA(1,1)-AA(2,2))/(2*AA(1,2));
        t = sign(tau)/(abs(tau)+sqrt(1+tau^2));
        c = 1/sqrt(1+t^2);
        s = t*c;
        V = [c,-s;s,c];
        A(:,idx) = A(:,idx)*V;
    else % two-side for EVP
        AA = A(idx,idx);
        tau = (AA(1,1)-AA(2,2))/(2*AA(1,2));
        t = sign(tau)/(abs(tau)+sqrt(1+tau^2));
        c = 1/sqrt(1+t^2);
        s = t*c;
        V = [c,-s;s,c];
        A(:,idx) = A(:,idx)*V;
        A(idx,:) = V'*A(idx,:);
        A(idx,idx) = (A(idx,idx)+A(idx,idx)')/2;
    end
end
