import pandas as pd
from matplotlib import pyplot as plt

tcgaMAF = pd.read_csv('../../data_tables/maf_files/pre_processed/Staudt_39_TCGA-DLBC_phs000178.aggregated.maf', sep='\t', low_memory=False)
drivers = pd.read_csv('../../data_tables/mutsig2cv_gistic_qvalues/DLBCL_550_training_noPDE4DIP_noHISTartifacts.sig_genes.txt', sep='\t', index_col=0)
drivers = drivers.loc[drivers['q'] <= 0.10].index

ys = tcgaMAF['n_alt_count'] / (tcgaMAF['n_alt_count'] + tcgaMAF['n_ref_count'])
xs = tcgaMAF['t_alt_count'] / (tcgaMAF['t_alt_count'] + tcgaMAF['t_ref_count'])

drivers = [x for x in drivers if (x != 'OR8I2') and (x != 'MUC6')]

driverMAF = tcgaMAF[tcgaMAF['Hugo_Symbol'].isin(drivers)]
driverMAF.index = driverMAF['Hugo_Symbol']

ys_drivers = driverMAF['n_alt_count'] / (driverMAF['n_alt_count'] + driverMAF['n_ref_count'])
xs_drivers = driverMAF['t_alt_count'] / (driverMAF['t_alt_count'] + driverMAF['t_ref_count'])

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7))
ax1.scatter(xs, ys, alpha=0.45, s=5, c='#00498E')
ax2.scatter(xs_drivers, ys_drivers, alpha=0.45, s=5, c='#00498E')
for i, txt in enumerate(driverMAF.index):
    if ys_drivers[i] > 0.20:
        ax2.annotate(txt, (xs_drivers[i], ys_drivers[i]), fontsize=6, xytext=(xs_drivers[i]+0.01, ys_drivers[i]))
ax2.set_ylim([0, 1.05])
ax1.set_ylim([0, 1.05])
ax2.set_xlim([0, 1])
ax1.set_title('All Mutations', size=16)
ax2.set_title('Driver Mutations', size=16)
fig.subplots_adjust(bottom=0.15, left=0.15)
fig.text(0.5, 0.02, 'Tumor VAF', ha='center', size=16)
fig.text(0.05, 0.37, 'Normal VAF', ha='center', size=16, rotation=90)

plt.savefig('../../plots/paper_figures/vaf_plots/tcga_vaf.pdf', format='pdf')
plt.savefig('../../plots/tcga_vaf.png', format='png')