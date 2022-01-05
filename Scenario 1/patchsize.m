function patchsize
%for creating the patch size figure in scenario 1
x = [1 2 3 4 5 6 7 8 9 10];
y = [90.74 68.96 62.22 45.53 34.74 28.99 23.66 19.96 19.71 17.21];

figure(5)
plot(x,y,'b-',x,y,'b*')
set(gca, 'FontSize', 18)
xlabel('Patch size')
ylabel('% of ticks remaining')
title(['Patch Size vs % of Ticks Remaining'],'FontSize', 20)
xlim([1,10])