rng(1);
close all
n = 1000;
Q = orth(randn(n,n));
alpha = 1;
iterMax = 300000;
iterplt = 100000;
err = zeros(2,iterMax);

eval = sort([1+exp(-alpha*(1:n/2)),linspace(2,3,n/2)])';
A = Q*diag(eval)*Q';
A = (A+A')/2;
[Vs,~] = eig(single(A));
[Vs,~] = qr(double(Vs));
A = Vs'*A*Vs;
A0 = (A+A')/2;
[~,idx] = sort(diag(A0));
A0 = A0(idx,idx);
nA = max(eval);
A = A0;
for iter = 1:iterMax
    err(1,iter) = norm(sort(diag(A))-eval)/nA;
    offA = abs(A-diag(diag(A)));
    [val, linearIndex] = max(offA(:));   
    [i, j] = ind2sub(size(A), linearIndex);
    A = jacobi(A, [i,j],2);
    if mod(iter,iterplt)==0
        iter
        figure
        imagesc(log(abs(A))./log(10));
        clim([-12,-6]);
        colorbar
        set(gcf, 'Color', 'w');
        title(['Classical: iteration=',num2str(iter)],'FontSize',14);
        export_fig (['fig\histc',num2str(iter/iterplt),'.eps']);
        export_fig (['fig\histc',num2str(iter/iterplt),'.pdf']);
    end
end

A = A0;
for iter = 1:iterMax
    err(2,iter) = norm(sort(diag(A))-eval)/nA;
    offA = abs(A-diag(diag(A)));
    den = abs(diag(A)-diag(A)');
    den = sqrt(den.^2+4*offA.^2)+den;
    offA = 2*offA.^2./den;
    [val, linearIndex] = max(offA(:));   
    [i, j] = ind2sub(size(A), linearIndex);
    A = jacobi(A, [i,j],2);
    if mod(iter,iterplt)==0
        iter
        figure
        imagesc(log(abs(A))./log(10));
        clim([-12,-6]);
        colorbar
        set(gcf, 'Color', 'w');
        title(['New: iteration=',num2str(iter)],'FontSize',14);
        export_fig (['fig\histn',num2str(iter/iterplt),'.eps']);
        export_fig (['fig\histn',num2str(iter/iterplt),'.pdf']);
    end
end
