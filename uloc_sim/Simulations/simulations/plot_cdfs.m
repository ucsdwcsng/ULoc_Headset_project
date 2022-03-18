%% Plot all errors
% err1 is Dloc error, act_fin_error is spotfi error
clearvars
close all
path = '/media/user1/easystore/datasets/aditya_sim/';
files = {'atk_cart_errors.mat', 'atk_vr_errors.mat', 'atk_spotfi_errors_w_cart.mat'};
cart_err = load([path, files{1}]);
vr_err = load([path, files{2}]);
spotfi_err_w_cart = load([path, files{3}]);
%%
figure(1)
h(1,1) = cdfplot(spotfi_err_w_cart.act_fin_error), hold on;
xlabel('Localization Error(m)'), ylabel('CDF'), hold on;
h(1,2) = cdfplot(cart_err.err1), hold on;

figure(1)
h(2, 1) = cdfplot(vr_err.act_fin_error), hold on;
h(2, 2) = cdfplot(vr_err.err1), hold on;
xlim([0, 4])

set( h(1,1), 'LineStyle', '-', 'Color', 'b');
set( h(2,1), 'LineStyle', '--', 'Color', 'b');
set( h(1,2), 'LineStyle', '-', 'Color', 'r');
set( h(2,2), 'LineStyle', '--', 'Color', 'r');
legend('SpotFi error compared against Cartographer', 'DLoc error trained with Cartographer', ...
        'SpotFi error compared against VR', 'DLoc error trained with VR')
title('Error comparisons in Atkinson')
%% plot for jacobs
files = {'jacobs_cart_errors.mat', 'jacobs_vr_errors.mat', 'jacobs_spotfi_errors_w_cart.mat'};
cart_err2 = load([path, files{1}]);
vr_err2 = load([path, files{2}]);
spotfi_err_w_cart2 = load([path, files{3}]);
%%
figure(2)
h(1,1) = cdfplot(spotfi_err_w_cart2.act_fin_error), hold on;
xlabel('Localization Error(m)'), ylabel('CDF'), hold on;
h(1,2) = cdfplot(cart_err2.err1), hold on;

figure(2)
h(2, 1) = cdfplot(vr_err2.act_fin_error), hold on;
h(2, 2) = cdfplot(vr_err2.err1), hold on;


set( h(1,1), 'LineStyle', '-', 'Color', 'b');
set( h(2,1), 'LineStyle', '--', 'Color', 'b');
set( h(1,2), 'LineStyle', '-', 'Color', 'r');
set( h(2,2), 'LineStyle', '--', 'Color', 'r');
legend('SpotFi error compared against Cartographer', 'DLoc error trained with Cartographer', ...
        'SpotFi error compared against VR', 'DLoc error trained with VR')
title('Error comparisons in Jacobs')
%%
figure(3)
all_spotfi_cart_err = [spotfi_err_w_cart.act_fin_error, spotfi_err_w_cart2.act_fin_error];
all_dloc_cart_err = [cart_err.err1, cart_err2.err1];
all_spotfi_vr_err = [vr_err.act_fin_error, vr_err2.act_fin_error];
all_dloc_vr_err = [vr_err.err1, vr_err2.err1];

h(1,1) = cdfplot(all_spotfi_cart_err), hold on;
xlabel('Localization Error(m)'), ylabel('CDF'), hold on;
h(1,2) = cdfplot(all_dloc_cart_err), hold on;

figure(3)
h(2, 1) = cdfplot(all_spotfi_vr_err), hold on;
h(2, 2) = cdfplot(all_dloc_vr_err), hold on;
xlim([0, 4])

set( h(1,1), 'LineStyle', '-', 'Color', 'b');
set( h(2,1), 'LineStyle', '--', 'Color', 'b');
set( h(1,2), 'LineStyle', '-', 'Color', 'r');
set( h(2,2), 'LineStyle', '--', 'Color', 'r');
legend('SpotFi error compared against Cartographer', 'DLoc error trained with Cartographer', ...
        'SpotFi error compared against VR', 'DLoc error trained with VR')
title('All error comparisons')
% /media/user1/easystore/DeepLoc/results/e301_real2dec_train/rw_test/untitled.m