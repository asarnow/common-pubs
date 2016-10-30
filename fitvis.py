#!/usr/bin/env python2.7
from fitness import *
from collections import defaultdict


def main(args):
    aa = pd.read_hdf(args.input[0], key='AA')
    cd = pd.read_hdf(args.input[0], key='Codon')
    bc = pd.read_hdf(args.input[0], key='Barcode')
    fits = {'AA': aa, 'Codon': cd, 'Barcode': bc}

    meta = parse_library('data/allele_dic_with_WT.pkl')
    idx = meta.reset_index().set_index(['Pos', 'AA', 'Codon'])  # Metadata index.

    pos = list(set(meta['Pos']) - {None, np.nan, 0})  # Unique Ub positions.
    aa = list(set(meta['AA']) - {None, np.nan, ""})  # Unique amino acids.
    aa_hydro = ['*'] + sorted(hydrophobicity(), key=hydrophobicity().get, reverse=True)
    aa_iggy = ['*', 'W', 'F', 'Y', 'L', 'I', 'M', 'V', 'C', 'A', 'G', 'P', 'S', 'T', 'N', 'Q', 'H', 'R', 'K', 'D', 'E']
    aa_iggy.reverse()

    cod = list(set(meta['Codon']) - {None, np.nan, ""})  # Unique codons.
    # cod = sorted(cod, key=lambda x: gclut[x])  # Lexigraphic by AA.
    cod = sorted(cod, key=lambda x: aa_iggy.index(gclut[x]))  # To match aa_iggy.

    control = 'control'
    experiments = fits['AA'].items.tolist()
    experiments.remove(control + 'R1')
    experiments = [ex[:-2] for ex in experiments]
    experiments = list(set(experiments))
    if args.exp is not None:
        experiments = [args.exp]

    score = 'slope'

    # TODO Sorted list of positions by: mean, max, std, min, isnan of mean/med/etc. maps.
    # TODO Where isnan set in all experiments.
    # TODO Exclude non-surface residues and repeat.
    # TODO Group residues by type and repeat.

    topdir = args.output[0]

    resmap = defaultdict(dict)
    adjmap = defaultdict(dict)
    resmean = defaultdict(dict)
    adjmean = defaultdict(dict)
    r1 = defaultdict(dict)
    r2 = defaultdict(dict)
    ctrl = {k: fits[k][control + "R1"] for k in fits}
    adjusted = defaultdict(dict)
    expt = defaultdict(dict)

    resmap[control]['AA'] = compute_hmap(ctrl['AA'][score], pos, aa_iggy, 'Pos', 'AA', np.nanmedian)
    resmap[control]['Codon'] = compute_hmap(ctrl['Codon'][score], pos, cod, 'Pos', 'Codon', np.nanmedian)
    resmap[control]['Barcode'] = compute_hmap(
                ctrl['Barcode'][score].reset_index().groupby(['Pos', 'AA']).mean()[score],
                pos, aa_iggy, 'Pos', 'AA', np.nanmedian)
    resmean[control] = {k: np.nanmean(resmap[control][k], axis=0) for k in resmap[control]}

    # Compile data and calculate heat map matrices for plotting.
    for ex in experiments:
        for k in fits:
            if k == "AA":
                r1[ex][k] = fits[k][ex + "R1"][score]
                r2[ex][k] = fits[k][ex + "R2"][score]
                expt[ex][k] = (r1[ex][k] + r2[ex][k]) / 2
                adjusted[ex][k] = expt[ex][k] - ctrl[k][score]
                resmap[ex][k] = compute_hmap(expt[ex][k], pos, aa_iggy, 'Pos', 'AA', np.nanmedian)
                adjmap[ex][k] = compute_hmap(adjusted[ex][k], pos, aa_iggy, 'Pos', 'AA', np.nanmedian)
            elif k == "Codon":
                r1[ex][k] = fits[k][ex + "R1"][score]
                r2[ex][k] = fits[k][ex + "R2"][score]
                expt[ex][k] = (r1[ex][k] + r2[ex][k]) / 2
                adjusted[ex][k] = expt[ex][k] - ctrl[k][score]
                resmap[ex][k] = compute_hmap(expt[ex][k], pos, cod, 'Pos', 'Codon', np.nanmedian)
                adjmap[ex][k] = compute_hmap(adjusted[ex][k], pos, cod, 'Pos', 'Codon', np.nanmedian)
            elif k == "Barcode":
                r1[ex][k] = fits[k][ex + "R1"][score].reset_index().groupby(['Pos', 'AA']).mean()[score]
                r2[ex][k] = fits[k][ex + "R2"][score].reset_index().groupby(['Pos', 'AA']).mean()[score]
                expt[ex][k] = (r1[ex][k] + r2[ex][k]) / 2
                adjusted[ex][k] = expt[ex][k] - ctrl[k][score].reset_index().groupby(['Pos', 'AA']).mean()[score]
                resmap[ex][k] = compute_hmap(expt[ex][k], pos, aa_iggy, 'Pos', 'AA', np.nanmedian)
                adjmap[ex][k] = compute_hmap(adjusted[ex][k], pos, aa_iggy, 'Pos', 'AA', np.nanmedian)
            resmean[ex] = {k: np.nanmean(resmap[ex][k], axis=0) for k in resmap[ex]}
            adjmean[ex] = {k: np.nanmean(adjmap[ex][k], axis=0) for k in adjmap[ex]}

    # Plot data and create image files.
    # Compare control standard errors.
    draw_stderr_hist(ctrl['AA']['stderr'].values,
                     ctrl['Codon']['stderr'].values,
                     ctrl['Barcode']['stderr'].values,
                     None,
                     os.path.join(topdir, control + '_stderr' + '.png'))
    # Compare replicates.
    for ex in r1:
        aa = (r1[ex]['AA'], r2[ex]['AA'])
        cd = (r1[ex]['Codon'], r2[ex]['Codon'])
        bc = (r1[ex]['Barcode'], r2[ex]['Barcode'])
        draw_replicate_corr(aa, cd, bc, None, os.path.join(topdir, ex + '_daycorr' + '.png'))
    # Fitness histograms.
    for ex in expt:
        draw_fitness_hist(expt[ex]['AA'].values,
                          expt[ex]['Codon'].values,
                          expt[ex]['Barcode'].values,
                          None,
                          os.path.join(topdir, ex + '_fithist.png'))
    # Experiment heat maps (including control).
    for ex in resmap:
        for k in resmap[ex]:
            if k == "AA" or k == "Barcode":
                draw_hmap(resmap[ex][k], aa_iggy, os.path.join(topdir, k + '_resmap_' + ex + '.png'))
            if k == "Codon":
                draw_hmap(resmap[ex][k], cod, os.path.join(topdir, k + '_resmap_' + ex + '.png'))
    # Adjusted (control-subtracted) heat maps.
    for ex in adjmap:
        for k in adjmap[ex]:
            if k == "AA" or k == "Barcode":
                draw_hmap(adjmap[ex][k], aa_iggy, os.path.join(topdir, k + '_adjmap_' + ex + '.png'))
            if k == "Codon":
                draw_hmap(adjmap[ex][k], cod, os.path.join(topdir, k + '_adjmap_' + ex + '.png'))

    return 0


def draw_stderr_hist(aa, cd, bc, wd, fname=None):
    fig = plt.figure(figsize=(20, 10))
    plt.subplot(2, 2, 1)
    weights = np.ones_like(aa) / float(len(aa))
    plt.hist(aa, 50, weights=weights)
    plt.title('By Residue', size=28)
    plt.xlabel('Standard Error', size=28)
    plt.ylabel('Frequency', size=28)
    plt.xticks(size=20)
    plt.yticks(size=20)

    plt.subplot(2, 2, 2)
    weights = np.ones_like(cd) / float(len(cd))
    plt.hist(cd, 50, weights=weights)
    plt.title('By Codon', size=28)
    plt.xlabel('Standard Error', size=28)
    plt.ylabel('Frequency', size=28)
    plt.xticks(size=20)
    plt.yticks(size=20)

    plt.subplot(2, 2, 3)
    weights = np.ones_like(bc) / float(len(bc))
    plt.hist(bc, 50, weights=weights)
    plt.title('By Lineage', size=28)
    plt.xlabel('Standard Error', size=28)
    plt.ylabel('Frequency', size=28)
    plt.xticks(size=20)
    plt.yticks(size=20)
    plt.tight_layout(h_pad=1.5)

    if wd is not None:
        plt.subplot(2, 2, 4)
        weights = np.ones_like(wd[np.isfinite(wd)]) / float(len(wd[np.isfinite(wd)]))
        plt.hist(wd[np.isfinite(wd)], 50, weights=weights)
        plt.title('By Lineage, Weighted Average', size=28)
        plt.xlabel('Standard Deviation', size=28)
        plt.ylabel('Frequency', size=28)
        plt.xticks(size=20)
        plt.yticks(size=20)

    plt.tight_layout(h_pad=1.5)

    if fname is not None:
        plt.savefig(fname, bbox_inches='tight')
        plt.close(fig)


def draw_replicate_corr(aa, cd, bc, wd, fname=None):
    fig = plt.figure(figsize=(20, 15))
    plt.subplot(2, 2, 1)
    seaborn.regplot(aa[0], aa[1], robust=False, ci=68)
    plt.title('By Residue', size=28)
    plt.xlabel('Fitness (Replicate 1)', size=28)
    plt.ylabel('Fitness (Replicate 2)', size=28)
    plt.xticks(size=20)
    plt.yticks(size=20)

    plt.subplot(2, 2, 2)
    seaborn.regplot(cd[0], cd[1], robust=False, ci=68)
    plt.title('By Codon', size=28)
    plt.xlabel('Fitness (Replicate 1)', size=28)
    plt.ylabel('Fitness (Replicate 2)', size=28)
    plt.xticks(size=20)
    plt.yticks(size=20)

    plt.subplot(2, 2, 3)
    seaborn.regplot(bc[0], bc[1], robust=False, ci=68)
    plt.title('By Lineage', size=28)
    plt.xlabel('Fitness (Replicate 1)', size=28)
    plt.ylabel('Fitness (Replicate 2)', size=28)
    plt.xticks(size=20)
    plt.yticks(size=20)

    if wd is not None:
        plt.subplot(2, 2, 4)
        seaborn.regplot(wd[0], wd[1], robust=False, ci=68)
        plt.title('By Lineage, Weighted Average', size=28)
        plt.xlabel('Fitness (Replicate 1)', size=28)
        plt.ylabel('Fitness (Replicate 2)', size=28)
        plt.xticks(size=20)
        plt.yticks(size=20)

    plt.tight_layout(h_pad=1.5)

    if fname is not None:
        plt.savefig(fname, bbox_inches='tight')
        plt.close(fig)


def draw_fitness_hist(dat1, dat2, dat3, dat4, fname=None):
    fig = plt.figure(figsize=(20, 10))
    plt.subplot(2, 2, 1)
    hdat = dat1[np.isfinite(dat1)]
    weights = np.ones_like(hdat) / float(len(hdat))
    plt.hist(hdat, 50, weights=weights)
    plt.title('By Residue', size=28)
    plt.xlabel('Fitness', size=28)
    plt.ylabel('Frequency', size=28)
    plt.xticks(size=20)
    plt.yticks(size=20)

    plt.subplot(2, 2, 2)
    hdat = dat2[np.isfinite(dat2)]
    weights = np.ones_like(hdat) / float(len(hdat))
    plt.hist(hdat, 50, weights=weights)
    plt.title('By Codon', size=28)
    plt.xlabel('Fitness', size=28)
    plt.ylabel('Frequency', size=28)
    plt.xticks(size=20)
    plt.yticks(size=20)

    plt.subplot(2, 2, 3)
    hdat = dat3[np.isfinite(dat3)]
    weights = np.ones_like(hdat) / float(len(hdat))
    plt.hist(hdat, 50, weights=weights)
    plt.title('By Lineage', size=28)
    plt.xlabel('Fitness', size=28)
    plt.ylabel('Frequency', size=28)
    plt.xticks(size=20)
    plt.yticks(size=20)

    if dat4 is not None:
        plt.subplot(2, 2, 4)
        hdat = dat4[np.isfinite(dat4)]
        weights = np.ones_like(hdat) / float(len(hdat))
        plt.hist(hdat, 50, weights=weights)
        plt.title('By Lineage, Weighted Average', size=28)
        plt.xlabel('Fitness', size=28)
        plt.ylabel('Frequency', size=28)
        plt.xticks(size=20)
        plt.yticks(size=20)

    plt.tight_layout(h_pad=1.5)
    if fname is not None:
        plt.savefig(fname, bbox_inches='tight')
        plt.close(fig)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser("Map read sequences to known barcodes.")
    parser.add_argument("-p", "--processes", action="store", type=int, default=1,
                        dest="numproc", help="Number of parallel processes.")
    parser.add_argument("-b", "--bardefs", action="store", type=str, default=None,
                        dest="bardefs", help="Barcode definitions (allele dictionary).")
    parser.add_argument("-d", "--expdefs", action="store", type=str, default=None,
                        dest="expdefs", help="Experiment definitions.")
    parser.add_argument("-e", "--experiment", action="store", type=str, default=None,
                        dest="exp", metavar="EXPNAME", help="Experiment label.")
    parser.add_argument("-s", "--score", action="store", dest="score", type=str,
                        default='slope', metavar="SCORE", help="Fitness score label (e.g. 'slope').")
    parser.add_argument("--prefix", action="store", type=str, default="",
                        metavar="PREFIX", help="Prefix for output figure files.")
    parser.add_argument("input", nargs=1, metavar="COUNTS")
    parser.add_argument("output", nargs=1, metavar="FITNESS")
    sys.exit(main(parser.parse_args()))
