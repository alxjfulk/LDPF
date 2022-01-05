%for creating the Time between burns figure
NB = [1:1:5];
NT = 100*[8.7573e+04/5.0887e+05 1.6995e+05/5.0887e+05 1.2289e+05/5.0887e+05 2.0083e+05/5.0887e+05 3.5741e+05/5.0887e+05];
plot(NB,NT,'-or')
set(gca, 'FontSize', 20)
ylim([0 100])
xlabel('Time Between Burns')
ylabel('% of Ticks Remaining')
