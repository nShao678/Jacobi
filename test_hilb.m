rng(1);
close all
n = 100;
A = zeros(n,n);
for ii = 1:n
    for jj = 1:n
        A(ii,jj) = sqrt(2*jj-1)/ii*prod(1-jj./(ii+1:ii+jj-1));
    end
end

A0 = A;
digits(20);            
H = sym(A);           
H_vpa = vpa(H);  
sval = sort(double(svd(H_vpa)));


% [~,~,Vs] = svd(single(A));
% % [Us,~] = qr(double(Us),'econ');
% [Vs,~] = qr(double(Vs));
% A0 = A0*Vs;


figure
hold on
plot(sval,'k-','LineWidth',2,'DisplayName','ref')

set(gca,'yscale','log')
iterMax = 300000;
A = A0;
for iter = 1:iterMax

    % AA = A'*A;
    % offA = abs(AA-diag(diag(AA)));
    % offA = offA./diag(AA)+diag(AA).\offA;
    % [val, idx] = max(offA(:));   
    % 
    % [i, j] = ind2sub(size(AA), idx);
    % A = jacobi(A, [i,j],1);
    A = jacobi(A, randperm(n,2),1);
    if iter == 50000
        plot(sort(vecnorm(A)),'-','LineWidth',2,'DisplayName','Randomized-50000')
    elseif iter == 100000
        plot(sort(vecnorm(A)),'-','LineWidth',2,'DisplayName','Randomized-100000')
    elseif iter == 200000
        plot(sort(vecnorm(A)),'-','LineWidth',2,'DisplayName','Randomized-200000')
    end
end

iterMax = 100000;
A = A0;
for iter = 1:iterMax

    AA = A'*A;
    offA = abs(AA-diag(diag(AA)));
    den = abs(diag(AA)-diag(AA)')+eye(n);
    den = sqrt(den.^2+4*offA.^2)+den;
    offA = 2*offA.^2./den;
    offA = offA./diag(AA)+diag(AA).\offA;
%     offA = offA.^2./den;
    [val, idx] = max(offA(:));   
    [i, j] = ind2sub(size(AA), idx);
    A = jacobi(A, [i,j],1);
    if iter == 1000
        plot(sort(vecnorm(A)),'--','LineWidth',2,'DisplayName','New-1000')
    elseif iter == 10000
        plot(sort(vecnorm(A)),'--','LineWidth',2,'DisplayName','New-10000')
    elseif iter == 50000
        plot(sort(vecnorm(A)),'--','LineWidth',2,'DisplayName','New-100000')
    end
%     err(2,iter) = norm(sort(vecnorm(A))-sval);
end
axis([-inf,100,0,1])
set(gcf, 'Color', 'w');
legend('FontSize',14,'box','off','Location','southeast')

hold off

export_fig (['fig\hilbh.eps']);
export_fig (['fig\hilbh.pdf']);

figure 
hold on
plot(sval,'k-','LineWidth',2,'DisplayName','ref')

plot(sort(svd(A0)),'b-','LineWidth',2,'DisplayName','svd')

plot(sort(vecnorm(A)),'r-','LineWidth',2,'DisplayName','Jacobi')

hold off

axis([-inf,100,0,1])
set(gca,'yscale','log')
set(gcf, 'Color', 'w');
legend('FontSize',14,'box','off','Location','southeast')
export_fig (['fig\hilb.eps']);
export_fig (['fig\hilb.pdf']);
% axis([-inf,40,0,1e-30])
% export_fig (['fig\hilbz.eps']);
% export_fig (['fig\hilbz.pdf']);