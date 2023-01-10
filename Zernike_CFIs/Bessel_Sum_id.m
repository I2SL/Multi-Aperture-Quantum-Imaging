x = linspace(0,100,10000);
N_MAX = [10,50,100,500,1000];

figure;
hold on
for n_max = N_MAX
    
    n = repmat((0:(n_max))',[1,length(x)]);
    xx = repmat(x,[n_max+1,1]);
    
    F = sum( (besselj(n-1,xx) - besselj(n+3,xx)).^2,1);
    plot(x,F)
end
hold off

xlabel('$x$','interpreter','latex')
ylabel('$\sum_{n=0}^{n_{max}}[J_{n-1}(x) - J_{n+3}(x)]^2$','interpreter','latex')
leg = legend(string(num2cell(N_MAX)));
title(leg,'$n_{max}$','interpreter','latex')



figure
hold on
for n_max = N_MAX
    
    n = repmat((0:(n_max))',[1,length(x)]);
    xx = repmat(x,[n_max+1,1]);
    
    CFI_theta = 4/3 * sum( n.*(n+2).*(besselj(n,xx) + besselj(n+2,xx)).^2,1);
    plot(x,CFI_theta)
end
hold off

xlabel('$x$','interpreter','latex')
ylabel('$\frac{4}{3} \sum_{n=0}^{n_{max}}n(n+2)[J_{n}(x) + J_{n+2}(x)]^2$','interpreter','latex')
leg = legend(string(num2cell(N_MAX)));
title(leg,'$n_{max}$','interpreter','latex')