a_max = 1000;
alpha = linspace(0,2,a_max);
N_MAX = [1,2,3,4,5,10];
QFI = 1/4 * 4;

% plot QFI for the radial parameter
figure(1)
plot(alpha, ones(size(alpha)),'black','LineWidth',2)
hold on
for n_max = N_MAX

n = repmat((0:(n_max-1))',[1,a_max]);
a = repmat(alpha,[n_max,1]);


%% compute CFI for truncated Zernike basis up to j_max
CFI_a_0 = @(n,a)  ( abs( besselj(n-1,a) - besselj(n+3,a) ).^2 );
CFI_zern = CFI_a_0(n,a);
CFI_zern = sum(CFI_zern,1);


%% plot fractional CFI
plot(alpha, CFI_zern/QFI,'--','LineWidth',1.5)

end
hold off

title({'Zernike Aperture Basis', 'Radial Parameter CFI $\mathcal{J}_{rr}$'},'interpreter','latex')
xlabel('$r/\sigma$','interpreter','latex')
ylabel('$\mathcal{J}_{rr}/\mathcal{K}_{rr}$','interpreter','latex')
leg = legend(['QFI',string(num2cell(N_MAX))]);
title(leg,'$n_{max}$','interpreter','latex')
axis 'square'