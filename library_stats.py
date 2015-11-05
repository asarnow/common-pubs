#!/usr/bin/env python
import cPickle as pic
%pylab


def hist2d_exact(x, y):
    xu = list(set(x))
    yu = list(set(y))
    valmap = lambda x : { x[i] : i for i in xrange(0,len(x)) }
    xmap = valmap(xu)
    ymap = valmap(yu)
    h = zeros([len(xu), len(yu)])
    for xi,yi in zip(x,y):
        h[xmap[xi], ymap[yi]] += 1
    return h


allele = pic.load(open('allele_dic.pkl','rb'))
translate = pic.load(open('translate.pkl','rb'))
aa2num = pic.load(open('aminotonumber.pkl','rb'))

codons = [x[0].split('_')[1] for x in allele.values()]
aa = [translate[x.replace('T','U')] for x in codons]
pos = [int(x[0].split('_')[0]) for x in allele.values()]
aa_num = [aa2num[x] for x in aa]
# hm = hist2d(pos, aa_num, [arange(2, 78), arange(0, 22)])
hm = hist2d_exact(aa_num, pos)
fig = figure()
pcolor(hm)
xlim(0, 76)
ylim(0, 21)
ax = gca()
fig.set_facecolor('white')
ax.set_yticks([x+0.6 for x in aa2num.values()])
ax.set_yticklabels(aa2num.keys())
ax.set_ylabel('Destination Codon')
ax.set_xlabel('Residue Number')
cb = colorbar()
cb.set_label('Number of Occurrences')

fig = figure()
fig.set_facecolor('white')
h = hist2d_exact(codons, pos)
pcolor(h)
xlim(0, 76)
ylim(0, 64)
xlabel('Residue Number')
ylabel('Destination Codon')
cb = colorbar()
cb.set_label('Number of Occurrences')
