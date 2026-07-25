rng(1);

iterMaxSet = [500000,1000000];
alphaSet = [0.01,1];
len = length(alphaSet);
errSet = cell(1,len);
erroffSet = cell(1,len);
n = 1000;
Q = orth(randn(n,n));
for ii = 1:len

iterMax = iterMaxSet(ii);
alpha = alphaSet(ii);
err = zeros(2,iterMax);
erroff = zeros(3,iterMax);
eval = sort([1+exp(-alpha*(1:n/2)),linspace(2,3,n/2)])';

A = Q*diag(eval)*Q';
A = (A+A')/2;
[Vs,~] = eig(single(A));
[Vs,~] = qr(double(Vs));
A = Vs'*A*Vs;
A0 = (A+A')/2;
nA = max(eval);
A = A0;
for iter = 1:iterMax
    if mod(iter,iterMax/100) ==0
        iter
    end
    err(1,iter) = norm(sort(diag(A))-eval)/nA;
    offA = abs(A-diag(diag(A)));
    erroff(1,iter) = norm(offA,'fro');
    [val, linearIndex] = max(offA(:));   
    [i, j] = ind2sub(size(A), linearIndex);
    A = jacobi(A, [i,j],2);
end

A = A0;
for iter = 1:iterMax    
    if mod(iter,iterMax/100) ==0
        iter
    end
    err(2,iter) = norm(sort(diag(A))-eval)/nA;
    offA = abs(A-diag(diag(A)));
    erroff(2,iter) = norm(offA,'fro');
    den = abs(diag(A)-diag(A)');
    den = sqrt(den.^2+4*offA.^2)+den;
    offA = 2*offA.^2./den;
    [erroff(3,iter), linearIndex] = max(offA(:));   
    [i, j] = ind2sub(size(A), linearIndex);
    A = jacobi(A, [i,j],2);
end
errSet{ii} = err;
erroffSet{ii} = erroff;

end
save('data')
%%

for ii = 1:len
err = errSet{ii};
erroff = erroffSet{ii};
    figure 
hold on
plot(err(1,:),'b-','LineWidth',2,'DisplayName','Classical')
plot(erroff(1,:),'b--','LineWidth',2,'DisplayName','Classical-off')
plot(err(2,:),'r-','LineWidth',2,'DisplayName','New')
plot(erroff(2,:),'r--','LineWidth',2,'DisplayName','New-off')
plot(erroff(3,:),'k--','LineWidth',2,'DisplayName','New-$L_{ij}$')
hold off
legend('FontSize',14,'box','off','Interpreter','latex')
set(gca,'yscale','log')
set(gcf, 'Color', 'w');
title(['$\alpha=',num2str(alphaSet(ii)),'$'],'Interpreter','latex','FontSize',14)
yticks(10.^(-16:-1))
ylim([1e-16,5e-5])
export_fig (['fig\clu',num2str(ii),'.eps'])
export_fig (['fig\clu',num2str(ii),'.pdf'])
end