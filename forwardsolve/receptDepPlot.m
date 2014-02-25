%% plot deposition in receptors 

figure(1)
bar(dep(1:nrec)*1e6)
set(gca,'XTick',[1:nrec]);
set(gca, 'XTickLabel', recept.label(1:nrec))
xlabel('Receptor'), ylabel('Amount deposited (mg)')
grid on
set(gca,'XGrid', 'off')
shg

print -depsc 'depbar.eps'
