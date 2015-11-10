import cPickle as pik

# group&replicate: [[T0index, gentime], [T1index, gentime], [T2index, gentime]]

groups = {'shmooR1':[['TTGACT', 0], ['GGAACT',1.87], ['TGACAT', 3.82]],
          'shmooR2':[['GGACGG', 0], ['CTCTAC', 1.9], ['GCGGAC', 3.6]],
          'etohR1':[['ATCGTG',0], ['TGAGTG', 3.14], ['CGCCTG', 5.14]],
          'etohR2':[['GCCATG',0], ['AAAATG', 1.76], ['TGTTGG', 4.02]],
          'onionR1':[['CGTGAT',0], ['ACATCG', 2.2], ['GCCTAA', 'nan']],
          'onionR2':[['TGGTCA',0], ['CACTGT', 2.6], ['ATTGGC', 3.8]],
          'whangeeR1':[['TTTCAC',0], ['GGCCAC',1.85], ['CGAAAC',4.10]],
          'whangeeR2':[['CGTACG',0], ['CCACTC',1.97], ['GCTACC',3.45]],
          'pyndR1':[['ATCAGT',0], ['GCTCAT', 2.43], ['AGGAAT', 4.11]],
          'pyndR2':[['CTTTTG',0], ['TAGTTG', 2.02], ['CCGGTG', 4.16]],
          'apcR1':[['GATCTG',0], ['TCAAGT',2.1], ['CTGATC',3.87]],
          'apcR2':[['AAGCTA',0], ['GTAGCC',2.1], ['TACAAG',3.87]],
          'controlR1':[['ATTCCG',0], ['AGCTAG',2.95], ['GTATAG',5.13]]}
output = open('groups.pkl', 'wb')
pik.dump(groups, output)
output.close()
