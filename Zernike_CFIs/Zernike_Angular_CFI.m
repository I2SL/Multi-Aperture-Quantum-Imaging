x = linspace(0,10,10000);
N_MAX = [1,2,3,4,5,10];
QFI_angular = x.^2;

figure;
hold on
plot(x,QFI_angular,'black','LineWidth',2)

for n_max = N_MAX
    
    n = repmat((0:(n_max))',[1,length(x)]);
    xx = repmat(x,[n_max+1,1]);
    
    CFI_theta = 4/3 * sum( n.*(n+2).*(besselj(n,xx) + besselj(n+2,xx)).^2,1);
    plot(x,CFI_theta)
end
hold off

title({'Zernike Aperture Basis', 'Angular Parameter CFI $\mathcal{J}_{\theta \theta}$'},'interpreter','latex')
xlabel('$r/\sigma$','interpreter','latex')
ylabel('$\frac{4}{3} \sum_{n=0}^{n_{max}}n(n+2)[J_{n}(x) + J_{n+2}(x)]^2$','interpreter','latex')
leg = legend(['QFI',string(num2cell(N_MAX))]);
title(leg,'$n_{max}$','interpreter','latex')