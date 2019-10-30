import sys
import csv

if __name__ == '__main__':
    if len(sys.argv) < 3 or len(sys.argv) > 5:
        print 'Usage: jaspar2gene <input> <output> [-col N]', \
        "\n\t<input> is a TSV file with a JASPAR identifier in column N", \
        "\n\t<output> is a copy of input with JASPAR identifier replaced by the gene name", \
        "\n\tIf column is unspecified any occurrence will be replaced."
        sys.exit(1)

    col = None
    J2TFmap = {
        'MA0002.2': 'RUNX1',
        'MA0003.3': 'TFAP2A',
        'MA0004.1': 'Arnt',
        'MA0006.1': 'Ahr::Arnt',
        'MA0007.3': 'Ar',
        'MA0009.2': 'T',
        'MA0014.2': 'PAX5',
        'MA0017.2': 'NR2F1',
        'MA0018.2': 'CREB1',
        'MA0019.1': 'Ddit3::Cebpa',
        'MA0024.3': 'E2F1',
        'MA0025.1': 'NFIL3',
        'MA0027.2': 'EN1',
        'MA0028.2': 'ELK1',
        'MA0029.1': 'Mecom',
        'MA0030.1': 'FOXF2',
        'MA0031.1': 'FOXD1',
        'MA0032.2': 'FOXC1',
        'MA0033.2': 'FOXL1',
        'MA0035.3': 'Gata1',
        'MA0036.2': 'GATA2',
        'MA0037.2': 'GATA3',
        'MA0038.1': 'Gfi1',
        'MA0039.2': 'Klf4',
        'MA0040.1': 'Foxq1',
        'MA0041.1': 'Foxd3',
        'MA0042.2': 'FOXI1',
        'MA0043.2': 'HLF',
        'MA0046.2': 'HNF1A',
        'MA0047.2': 'Foxa2',
        'MA0048.2': 'NHLH1',
        'MA0050.2': 'IRF1',
        'MA0051.1': 'IRF2',
        'MA0052.3': 'MEF2A',
        'MA0056.1': 'MZF1',
        'MA0057.1': 'MZF1(var.2)',
        'MA0058.3': 'MAX',
        'MA0059.1': 'MAX::MYC',
        'MA0060.2': 'NFYA',
        'MA0062.2': 'Gabpa',
        'MA0063.1': 'Nkx2-5',
        'MA0065.2': 'Pparg::Rxra',
        'MA0066.1': 'PPARG',
        'MA0067.1': 'Pax2',
        'MA0068.2': 'PAX4',
        'MA0069.1': 'Pax6',
        'MA0070.1': 'PBX1',
        'MA0071.1': 'RORA',
        'MA0072.1': 'RORA(var.2)',
        'MA0073.1': 'RREB1',
        'MA0074.1': 'RXRA::VDR',
        'MA0075.2': 'Prrx2',
        'MA0076.2': 'ELK4',
        'MA0077.1': 'SOX9',
        'MA0078.1': 'Sox17',
        'MA0079.3': 'SP1',
        'MA0080.4': 'SPI1',
        'MA0081.1': 'SPIB',
        'MA0083.3': 'SRF',
        'MA0084.1': 'SRY',
        'MA0087.1': 'Sox5',
        'MA0088.2': 'ZNF143',
        'MA0089.1': 'MAFG::NFE2L1',
        'MA0090.2': 'TEAD1',
        'MA0091.1': 'TAL1::TCF3',
        'MA0092.1': 'Hand1::Tcf3',
        'MA0093.2': 'USF1',
        'MA0095.2': 'YY1',
        'MA0098.3': 'ETS1',
        'MA0099.2': 'FOS::JUN',
        'MA0100.2': 'Myb',
        'MA0101.1': 'REL',
        'MA0102.3': 'CEBPA',
        'MA0103.2': 'ZEB1',
        'MA0104.3': 'Mycn',
        'MA0105.4': 'NFKB1',
        'MA0106.3': 'TP53',
        'MA0107.1': 'RELA',
        'MA0108.2': 'TBP',
        'MA0109.1': 'HLTF',
        'MA0111.1': 'Spz1',
        'MA0112.3': 'ESR1',
        'MA0113.3': 'NR3C1',
        'MA0114.3': 'Hnf4a',
        'MA0115.1': 'NR1H2::RXRA',
        'MA0116.1': 'Znf423',
        'MA0117.2': 'Mafb',
        'MA0119.1': 'NFIC::TLX1',
        'MA0122.2': 'NKX3-2',
        'MA0124.2': 'Nkx3-1',
        'MA0125.1': 'Nobox',
        'MA0130.1': 'ZNF354C',
        'MA0131.2': 'HINFP',
        'MA0132.2': 'PDX1',
        'MA0135.1': 'Lhx3',
        'MA0136.2': 'ELF5',
        'MA0137.3': 'STAT1',
        'MA0138.2': 'REST',
        'MA0139.1': 'CTCF',
        'MA0140.2': 'GATA1::TAL1',
        'MA0141.3': 'ESRRB',
        'MA0142.1': 'Pou5f1::Sox2',
        'MA0143.3': 'Sox2',
        'MA0144.2': 'STAT3',
        'MA0145.3': 'TFCP2',
        'MA0146.2': 'Zfx',
        'MA0147.2': 'Myc',
        'MA0148.3': 'FOXA1',
        'MA0149.1': 'EWSR1-FLI1',
        'MA0150.2': 'Nfe2l2',
        'MA0151.1': 'Arid3a',
        'MA0152.1': 'NFATC2',
        'MA0153.2': 'HNF1B',
        'MA0154.3': 'EBF1',
        'MA0155.1': 'INSM1',
        'MA0156.2': 'FEV',
        'MA0157.2': 'FOXO3',
        'MA0158.1': 'HOXA5',
        'MA0159.1': 'RARA::RXRA',
        'MA0160.1': 'NR4A2',
        'MA0161.1': 'NFIC',
        'MA0162.2': 'EGR1',
        'MA0163.1': 'PLAG1',
        'MA0164.1': 'Nr2e3',
        'MA0258.2': 'ESR2',
        'MA0259.1': 'ARNT::HIF1A',
        'MA0442.1': 'SOX10',
        'MA0461.2': 'Atoh1',
        'MA0462.1': 'BATF::JUN',
        'MA0463.1': 'Bcl6',
        'MA0464.2': 'BHLHE40',
        'MA0465.1': 'CDX2',
        'MA0466.2': 'CEBPB',
        'MA0467.1': 'Crx',
        'MA0468.1': 'DUX4',
        'MA0469.2': 'E2F3',
        'MA0470.1': 'E2F4',
        'MA0471.1': 'E2F6',
        'MA0472.2': 'EGR2',
        'MA0473.2': 'ELF1',
        'MA0474.2': 'ERG',
        'MA0475.2': 'FLI1',
        'MA0476.1': 'FOS',
        'MA0477.1': 'FOSL1',
        'MA0478.1': 'FOSL2',
        'MA0479.1': 'FOXH1',
        'MA0480.1': 'Foxo1',
        'MA0481.1': 'FOXP1',
        'MA0482.1': 'Gata4',
        'MA0483.1': 'Gfi1b',
        'MA0484.1': 'HNF4G',
        'MA0485.1': 'Hoxc9',
        'MA0486.2': 'HSF1',
        'MA0488.1': 'JUN',
        'MA0489.1': 'JUN(var.2)',
        'MA0490.1': 'JUNB',
        'MA0491.1': 'JUND',
        'MA0492.1': 'JUND(var.2)',
        'MA0493.1': 'Klf1',
        'MA0494.1': 'Nr1h3::Rxra',
        'MA0495.1': 'MAFF',
        'MA0496.1': 'MAFK',
        'MA0497.1': 'MEF2C',
        'MA0498.2': 'MEIS1',
        'MA0499.1': 'Myod1',
        'MA0500.1': 'Myog',
        'MA0501.1': 'MAF::NFE2',
        'MA0502.1': 'NFYB',
        'MA0503.1': 'Nkx2-5(var.2)',
        'MA0504.1': 'NR2C2',
        'MA0505.1': 'Nr5a2',
        'MA0506.1': 'NRF1',
        'MA0507.1': 'POU2F2',
        'MA0508.1': 'PRDM1',
        'MA0509.1': 'Rfx1',
        'MA0510.2': 'RFX5',
        'MA0511.2': 'RUNX2',
        'MA0512.2': 'Rxra',
        'MA0513.1': 'SMAD2::SMAD3::SMAD4',
        'MA0514.1': 'Sox3',
        'MA0515.1': 'Sox6',
        'MA0516.1': 'SP2',
        'MA0517.1': 'STAT1::STAT2',
        'MA0518.1': 'Stat4',
        'MA0519.1': 'Stat5a::Stat5b',
        'MA0520.1': 'Stat6',
        'MA0521.1': 'Tcf12',
        'MA0522.2': 'TCF3',
        'MA0523.1': 'TCF7L2',
        'MA0524.2': 'TFAP2C',
        'MA0525.2': 'TP63',
        'MA0526.1': 'USF2',
        'MA0527.1': 'ZBTB33',
        'MA0528.1': 'ZNF263',
        'MA0591.1': 'Bach1::Mafk',
        'MA0592.2': 'Esrra',
        'MA0593.1': 'FOXP2',
        'MA0594.1': 'Hoxa9',
        'MA0595.1': 'SREBF1',
        'MA0596.1': 'SREBF2',
        'MA0597.1': 'THAP1',
        'MA0598.2': 'EHF',
        'MA0599.1': 'KLF5',
        'MA0600.2': 'RFX2',
        'MA0601.1': 'Arid3b',
        'MA0602.1': 'Arid5a',
        'MA0603.1': 'Arntl',
        'MA0604.1': 'Atf1',
        'MA0605.1': 'Atf3',
        'MA0606.1': 'NFAT5',
        'MA0607.1': 'Bhlha15',
        'MA0608.1': 'Creb3l2',
        'MA0609.1': 'Crem',
        'MA0610.1': 'DMRT3',
        'MA0611.1': 'Dux',
        'MA0612.1': 'EMX1',
        'MA0613.1': 'FOXG1',
        'MA0614.1': 'Foxj2',
        'MA0615.1': 'Gmeb1',
        'MA0616.1': 'Hes2',
        'MA0617.1': 'Id2',
        'MA0618.1': 'LBX1',
        'MA0619.1': 'LIN54',
        'MA0620.1': 'Mitf',
        'MA0621.1': 'mix-a',
        'MA0622.1': 'Mlxip',
        'MA0623.1': 'Neurog1',
        'MA0624.1': 'NFATC1',
        'MA0625.1': 'NFATC3',
        'MA0626.1': 'Npas2',
        'MA0627.1': 'Pou2f3',
        'MA0628.1': 'POU6F1',
        'MA0629.1': 'Rhox11',
        'MA0630.1': 'SHOX',
        'MA0631.1': 'Six3',
        'MA0632.1': 'Tcfl5',
        'MA0633.1': 'Twist2',
        'MA0634.1': 'ALX3',
        'MA0635.1': 'BARHL2',
        'MA0636.1': 'BHLHE41',
        'MA0637.1': 'CENPB',
        'MA0638.1': 'CREB3',
        'MA0639.1': 'DBP',
        'MA0640.1': 'ELF3',
        'MA0641.1': 'ELF4',
        'MA0642.1': 'EN2',
        'MA0643.1': 'Esrrg',
        'MA0644.1': 'ESX1',
        'MA0645.1': 'ETV6',
        'MA0646.1': 'GCM1',
        'MA0647.1': 'GRHL1',
        'MA0648.1': 'GSC',
        'MA0649.1': 'HEY2',
        'MA0650.1': 'HOXA13',
        'MA0651.1': 'HOXC11',
        'MA0652.1': 'IRF8',
        'MA0653.1': 'IRF9',
        'MA0654.1': 'ISX',
        'MA0655.1': 'JDP2',
        'MA0656.1': 'JDP2(var.2)',
        'MA0657.1': 'KLF13',
        'MA0658.1': 'LHX6',
        'MA0659.1': 'MAFG',
        'MA0660.1': 'MEF2B',
        'MA0661.1': 'MEOX1',
        'MA0662.1': 'MIXL1',
        'MA0663.1': 'MLX',
        'MA0664.1': 'MLXIPL',
        'MA0665.1': 'MSC',
        'MA0666.1': 'MSX1',
        'MA0667.1': 'MYF6',
        'MA0668.1': 'NEUROD2',
        'MA0669.1': 'NEUROG2',
        'MA0670.1': 'NFIA',
        'MA0671.1': 'NFIX',
        'MA0672.1': 'NKX2-3',
        'MA0673.1': 'NKX2-8',
        'MA0674.1': 'NKX6-1',
        'MA0675.1': 'NKX6-2',
        'MA0676.1': 'Nr2e1',
        'MA0677.1': 'Nr2f6',
        'MA0678.1': 'OLIG2',
        'MA0679.1': 'ONECUT1',
        'MA0680.1': 'PAX7',
        'MA0681.1': 'Phox2b',
        'MA0682.1': 'Pitx1',
        'MA0683.1': 'POU4F2',
        'MA0684.1': 'RUNX3',
        'MA0685.1': 'SP4',
        'MA0686.1': 'SPDEF',
        'MA0687.1': 'SPIC',
        'MA0688.1': 'TBX2',
        'MA0689.1': 'TBX20',
        'MA0690.1': 'TBX21',
        'MA0691.1': 'TFAP4',
        'MA0692.1': 'TFEB',
        'MA0693.1': 'Vdr',
        'MA0694.1': 'ZBTB7B',
        'MA0695.1': 'ZBTB7C',
        'MA0696.1': 'ZIC1',
        'MA0697.1': 'ZIC3',
        'MA0698.1': 'ZBTB18',
        'MA0699.1': 'LBX2',
        'MA0700.1': 'LHX2',
        'MA0701.1': 'LHX9',
        'MA0702.1': 'LMX1A',
        'MA0703.1': 'LMX1B',
        'MA0704.1': 'Lhx4',
        'MA0705.1': 'Lhx8',
        'MA0706.1': 'MEOX2',
        'MA0707.1': 'MNX1',
        'MA0708.1': 'MSX2',
        'MA0709.1': 'Msx3',
        'MA0710.1': 'NOTO',
        'MA0711.1': 'OTX1',
        'MA0712.1': 'OTX2',
        'MA0713.1': 'PHOX2A',
        'MA0714.1': 'PITX3',
        'MA0715.1': 'PROP1',
        'MA0716.1': 'PRRX1',
        'MA0717.1': 'RAX2',
        'MA0718.1': 'RAX',
        'MA0719.1': 'RHOXF1',
        'MA0720.1': 'Shox2',
        'MA0721.1': 'UNCX',
        'MA0722.1': 'VAX1',
        'MA0723.1': 'VAX2',
        'MA0724.1': 'VENTX',
        'MA0725.1': 'VSX1',
        'MA0726.1': 'VSX2',
        'MA0727.1': 'NR3C2',
        'MA0728.1': 'Nr2f6(var.2)',
        'MA0729.1': 'RARA',
        'MA0730.1': 'RARA(var.2)',
        'MA0731.1': 'BCL6B',
        'MA0732.1': 'EGR3',
        'MA0733.1': 'EGR4',
        'MA0734.1': 'GLI2',
        'MA0735.1': 'GLIS1',
        'MA0736.1': 'GLIS2',
        'MA0737.1': 'GLIS3',
        'MA0738.1': 'HIC2',
        'MA0739.1': 'Hic1',
        'MA0740.1': 'KLF14',
        'MA0741.1': 'KLF16',
        'MA0742.1': 'Klf12',
        'MA0743.1': 'SCRT1',
        'MA0744.1': 'SCRT2',
        'MA0745.1': 'SNAI2',
        'MA0746.1': 'SP3',
        'MA0747.1': 'SP8',
        'MA0748.1': 'YY2',
        'MA0749.1': 'ZBED1',
        'MA0750.1': 'ZBTB7A',
        'MA0751.1': 'ZIC4',
        'MA0752.1': 'ZNF410',
        'MA0753.1': 'ZNF740',
        'MA0754.1': 'CUX1',
        'MA0755.1': 'CUX2',
        'MA0756.1': 'ONECUT2',
        'MA0757.1': 'ONECUT3',
        'MA0758.1': 'E2F7',
        'MA0759.1': 'ELK3',
        'MA0760.1': 'ERF',
        'MA0761.1': 'ETV1',
        'MA0762.1': 'ETV2',
        'MA0763.1': 'ETV3',
        'MA0764.1': 'ETV4',
        'MA0765.1': 'ETV5',
        'MA0766.1': 'GATA5',
        'MA0767.1': 'GCM2',
        'MA0768.1': 'LEF1',
        'MA0769.1': 'Tcf7',
        'MA0770.1': 'HSF2',
        'MA0771.1': 'HSF4',
        'MA0772.1': 'IRF7',
        'MA0773.1': 'MEF2D',
        'MA0774.1': 'MEIS2',
        'MA0775.1': 'MEIS3',
        'MA0776.1': 'MYBL1',
        'MA0777.1': 'MYBL2',
        'MA0778.1': 'NFKB2',
        'MA0779.1': 'PAX1',
        'MA0780.1': 'PAX3',
        'MA0781.1': 'PAX9',
        'MA0782.1': 'PKNOX1',
        'MA0783.1': 'PKNOX2',
        'MA0784.1': 'POU1F1',
        'MA0785.1': 'POU2F1',
        'MA0786.1': 'POU3F1',
        'MA0787.1': 'POU3F2',
        'MA0788.1': 'POU3F3',
        'MA0789.1': 'POU3F4',
        'MA0790.1': 'POU4F1',
        'MA0791.1': 'POU4F3',
        'MA0792.1': 'POU5F1B',
        'MA0793.1': 'POU6F2',
        'MA0794.1': 'PROX1',
        'MA0795.1': 'SMAD3',
        'MA0796.1': 'TGIF1',
        'MA0797.1': 'TGIF2',
        'MA0798.1': 'RFX3',
        'MA0799.1': 'RFX4',
        'MA0800.1': 'EOMES',
        'MA0801.1': 'MGA',
        'MA0802.1': 'TBR1',
        'MA0803.1': 'TBX15',
        'MA0804.1': 'TBX19',
        'MA0805.1': 'TBX1',
        'MA0806.1': 'TBX4',
        'MA0807.1': 'TBX5',
        'MA0808.1': 'TEAD3',
        'MA0809.1': 'TEAD4',
        'MA0810.1': 'TFAP2A(var.2)',
        'MA0811.1': 'TFAP2B',
        'MA0812.1': 'TFAP2B(var.2)',
        'MA0813.1': 'TFAP2B(var.3)',
        'MA0814.1': 'TFAP2C(var.2)',
        'MA0815.1': 'TFAP2C(var.3)',
        'MA0816.1': 'Ascl2',
        'MA0817.1': 'BHLHE23',
        'MA0818.1': 'BHLHE22',
        'MA0819.1': 'CLOCK',
        'MA0820.1': 'FIGLA',
        'MA0821.1': 'HES5',
        'MA0822.1': 'HES7',
        'MA0823.1': 'HEY1',
        'MA0824.1': 'ID4',
        'MA0825.1': 'MNT',
        'MA0826.1': 'OLIG1',
        'MA0827.1': 'OLIG3',
        'MA0828.1': 'SREBF2(var.2)',
        'MA0829.1': 'Srebf1(var.2)',
        'MA0830.1': 'TCF4',
        'MA0831.1': 'TFE3',
        'MA0832.1': 'Tcf21',
        'MA0833.1': 'ATF4',
        'MA0834.1': 'ATF7',
        'MA0835.1': 'BATF3',
        'MA0836.1': 'CEBPD',
        'MA0837.1': 'CEBPE',
        'MA0838.1': 'CEBPG',
        'MA0839.1': 'CREB3L1',
        'MA0840.1': 'Creb5',
        'MA0841.1': 'NFE2',
        'MA0842.1': 'NRL',
        'MA0843.1': 'TEF',
        'MA0844.1': 'XBP1',
        'MA0845.1': 'FOXB1',
        'MA0846.1': 'FOXC2',
        'MA0847.1': 'FOXD2',
        'MA0848.1': 'FOXO4',
        'MA0849.1': 'FOXO6',
        'MA0850.1': 'FOXP3',
        'MA0851.1': 'Foxj3',
        'MA0852.1': 'Foxk1',
        'MA0853.1': 'Alx4',
        'MA0854.1': 'Alx1',
        'MA0855.1': 'RXRB',
        'MA0856.1': 'RXRG',
        'MA0857.1': 'Rarb',
        'MA0858.1': 'Rarb(var.2)',
        'MA0859.1': 'Rarg',
        'MA0860.1': 'Rarg(var.2)',
        'MA0861.1': 'TP73',
        'MA0862.1': 'GMEB2',
        'MA0863.1': 'MTF1',
        'MA0864.1': 'E2F2',
        'MA0865.1': 'E2F8',
        'MA0866.1': 'SOX21',
        'MA0867.1': 'SOX4',
        'MA0868.1': 'SOX8',
        'MA0869.1': 'Sox11',
        'MA0870.1': 'Sox1',
        'MA0871.1': 'TFEC',
        'MA0872.1': 'TFAP2A(var.3)',
        'MA0873.1': 'HOXD12',
        'MA0874.1': 'Arx',
        'MA0875.1': 'BARX1',
        'MA0876.1': 'BSX',
        'MA0877.1': 'Barhl1',
        'MA0878.1': 'CDX1',
        'MA0879.1': 'Dlx1',
        'MA0880.1': 'Dlx3',
        'MA0881.1': 'Dlx4',
        'MA0882.1': 'DLX6',
        'MA0883.1': 'Dmbx1',
        'MA0884.1': 'DUXA',
        'MA0885.1': 'Dlx2',
        'MA0886.1': 'EMX2',
        'MA0887.1': 'EVX1',
        'MA0888.1': 'EVX2',
        'MA0889.1': 'GBX1',
        'MA0890.1': 'GBX2',
        'MA0891.1': 'GSC2',
        'MA0892.1': 'GSX1',
        'MA0893.1': 'GSX2',
        'MA0894.1': 'HESX1',
        'MA0895.1': 'HMBOX1',
        'MA0896.1': 'Hmx1',
        'MA0897.1': 'Hmx2',
        'MA0898.1': 'Hmx3',
        'MA0899.1': 'HOXA10',
        'MA0900.1': 'HOXA2',
        'MA0901.1': 'HOXB13',
        'MA0902.1': 'HOXB2',
        'MA0903.1': 'HOXB3',
        'MA0904.1': 'Hoxb5',
        'MA0905.1': 'HOXC10',
        'MA0906.1': 'HOXC12',
        'MA0907.1': 'HOXC13',
        'MA0908.1': 'HOXD11',
        'MA0909.1': 'HOXD13',
        'MA0910.1': 'Hoxd8',
        'MA0911.1': 'Hoxa11',
        'MA0912.1': 'Hoxd3',
        'MA0913.1': 'Hoxd9',
        'MA0914.1': 'ISL2',
        'MA1099.1': 'Hes1'}

    if len(sys.argv) > 3:
        if sys.argv[3] == '-col':
            col = int(sys.argv[4])
    with open(sys.argv[1]) as csvfile:
        r = csv.reader(csvfile, delimiter='\t')
        with open(sys.argv[2], 'wt') as tsvfile:
            w = csv.writer(tsvfile, delimiter='\t')
            for row in r:
                if col != None:
                    try:
                        newname = J2TFmap[row[col]]
                        row[col] = newname
                    except:
                        print 'Failed to replace', row[col]
                else:
                    for c in range(len(row)):
                        try:
                            newname = J2TFmap[row[c]]
                            row[c] = newname
                        except:
                            pass
                w.writerow(row)
        tsvfile.close()
    csvfile.close()