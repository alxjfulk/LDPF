%for creating the number of burns figure
NB = [1:1:9];
NT = 100*[4.4897e+05/5.0887e+05 4.1513e+05/5.0887e+05 3.7940e+05/5.0887e+05 3.4023e+05/5.0887e+05 2.9716e+05/5.0887e+05 2.5011e+05/5.0887e+05 1.9929e+05/5.0887e+05 1.4515e+05/5.0887e+05 8.7573e+04/5.0887e+05];
plot(NB,NT,'-*m')
set(gca, 'FontSize', 20)
ylim([0 100])
xlabel('Number of Burns')
ylabel('% of Ticks Remaining')
