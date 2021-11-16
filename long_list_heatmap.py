from glbase3 import *
import os

# config.draw_mode = "pdf"
os.chdir('/Users/jplab/Desktop/2021-11/data/11-15_cell_surface_protein')

matched_surface_proteins = './Tony_cell_surface_protein_list_'+'A'+'_match.tsv'


expn = expression(
        filename= matched_surface_proteins,
        format={'force_tsv': True, 'ensg': 6, 'name':7, 'skiplines': 0}, 
        expn='column[0:6]',
        # cond_names=[
        # # 'B14_0.4_rp1', 'B14_0.4_rp2', 'B14_0.8_rp1', 'B14_0.8_rp2', 'B14_1.2_rp1', 'B14_1.2_rp2', 'B14_1.6_rp1', 'B14_1.6_rp2'
        # ]
        )      

expn.coerce(int)

# remove if 5 conditions lower than 100
expn = expn.filter_low_expressed(200,4)


        # first log2 transform
        # second zscore
        # expn.log(2,.1)
        # expn_3groups.log(2,.1)

        # expn.heatmap(filename="%s_level_nozscore.pdf" % sample_name, heat_wid=0.25, grid=True,
        #     col_cluster=False, heat_hei=0.018*len(expn), colbar_label="Z-score_"+sample_name, row_cluster=True)

        # trans to z-score
expn.row_Z()

expn.saveTSV('surface_protein_matched'+"_zscore.tsv")

# # save to tsv for plotting 
# expn.saveTSV("3d_zscore_5groups_nomergereps.tsv")
# expn_3groups.saveTSV("3d_zscore_3groups_nomergereps.tsv")


# handle long name    
# print(sample_name+': ',len(sample_name))
sample_name = 'tony_cell_surface_protein_'+'A'
title_sample_name = ''

if len(sample_name) > 50:
    title_sample_name = sample_name[0:25] + '-' + '\n' + sample_name[25:50] + '-' + '\n' + sample_name[50:]
    
# handle long list heatmap height issue
expn_length = len(expn)
print(expn_length)

figsize=(12, 22)
heat_hei = 0.01 * expn_length
grid=True
if expn_length <= 50:
    heat_hei = 0.03 * expn_length
    figsize=(3, 7)
if expn_length >= 120:
    heat_hei = 0.0045 * expn_length
if expn_length >= 150:
    heat_hei = 0.0045 * expn_length
    figsize=(12, 26)
if expn_length >= 250:
    heat_hei = 0.002 * expn_length
    figsize=(12, 35)
if expn_length >= 500:
    heat_hei = 0.0013 * expn_length
    figsize=(12, 48)
    grid=False
if expn_length >= 800: 
    heat_hei = 0.0009 * expn_length
    figsize=(12, 50)
    grid=False   

# use this to add name to heatmap
highlights = [
    'CA12', 'ABCC2', 'ATP1B1', 'ITGB1', 'ABCC3', 'ABCC1', 'CD44', 'CNTN1', 'NRCAM', 'S1PR3', 'APP', 'EGFR', 'SLC38A2', 'CDH17', 'FURIN', 'ERMP1', 'CPD', 'SLC12A2', 'ABCC4', 'ADAM15', 'SLC3A2', 'TFPI', 'TMEM106B', 'SLC23A2', 
    'TNFRSF1A', 'PTPRM', 'TMEM123', 'SLC5A3', 'MEGF8', 'DCBLD2', 'PTPRF', 'PRNP', 'MRC2', 'AMIGO2', 'ITGA3', 'CD46', 'NT5E', 'MYOF', 'APLP2', 'SLC7A2', 'PLXNB2', 'PLA2R1', 'TNFSF15', 'PCDH9', 'TGFBR2', 'TM4SF1', 'TMEM132A', 
    'LGR4', 'OSMR', 'IL1R1', 'FAT1', 'HLA-A', 'SLC39A8', 'IGF1R', 'SLC7A5', 'TPCN1', 'MYADM', 'CRIM1', 'CNTNAP3B', 'TM4SF4', 'NRP2', 'PARM1', 'GABRB3', 'CLDN2', 'CSF1', 'MEGF9', 'SEMA4B', 'SPINT2', 'ITGB4', 'NRG1', 'CDH11',
    'LAMP2', 'DDR1', 'MPZL1', 'MUC13', 'CLSTN1', 'ATP1A1', 'CD38', 'HLA-DMB', 'BCAM', 'ADAM9', 'ABCA2', 'CEACAM6', 'CELSR1', 'PAM', 'MR1', 'DAG1', 'SLC1A5', 'IL6ST', 'CD59', 'TSPAN3', 'TSPAN4', 'PTPRK', 'LY6E', 'ATP13A3', 
    'SLC38A1', 'UNC93B1', 'ALCAM', 'TM4SF18', 'SLC16A7', 'ABCA3', 'EVA1C', 'FRAS1', 'CD63', 'SLC40A1', 'SLC2A1', 'VLDLR', 'CD55', 'BSG', 'ATP1B3', 'KITLG', 'THSD7A', 'LGR6', 'FGFR1', 'ATRN', 'NIPAL3', 'TGOLN2', 'QSOX1', 'HEG1',
    'GPR107', 'ITM2B', 'CD151', 'SPPL2A', 'ABCG2', 'CD70', 'MANSC1', 'RNF13', 'SLC39A14', 'PTPRU', 'ITFG1', 'CD164', 'ITM2C', 'GRAMD1B', 'PLXND1', 'ANTXR1', 'TMX4', 'NPTN', 'BAMBI', 'SLC16A1', 'CD24', 'ITGA6', 'SLC4A7',
    'ADAM10', 'LMAN2', 'CD9', 'ADCY9', 'LRP10', 'TTYH3', 'TMEM9', 'LTBR', 'SERINC3', 'SCNN1A', 'TMEM8A', 'ITGAV', 'SLC22A5', 'GINM1', 'HLA-C', 'TMED7', 'TM9SF4', 'MFSD6', 'PLXNA3', 'NTRK3', 'HAVCR1', 'DCBLD1', 'CDH1', 'TENM3',
    'TM9SF2', 'EPHA2', 'PIEZO1', 'CEACAM1', 'CLSTN3', 'DSG2', 'NCSTN', 'PODXL', 'GJA1', 'SLC30A1', 'SLC6A6', 'GPRC5A', 'NLGN2', 'CD47', 'GIPR', 'SERINC1', 'ITGB5', 'SLC6A12', 'CD276', 'NPC1', 'PVR', 'SUSD2', 'ADIPOR2', 'ANKH', 
    'TM7SF3', 'LRP12', 'PCDH1', 'CLSTN2', 'TM9SF3', 'SLC29A4', 'EFNA1', 'LDLRAD4', 'SLC47A1', 'ANO5', 'ROBO1', 'TGFBR3', 'FZD7', 'LRP5', 'SCARB1', 'TPBG', 'SLC44A1', 'CNNM2', 'TCTN3', 'LNPEP', 'NRP1', 'ATP2B4', 'ABCC5', 'SDC1', 
    'SLCO4A1', 'NOTCH2', 'GABRA5', 'STS', 'FLRT3', 'TXNDC15', 'LRP1', 'GPR37', 'F2RL1', 'LRP6', 'SLC12A6', 'APLP1', 'PIGT', 'TSPAN6', 'ATRAID', 'CELSR2', 'ACVR1B', 'NEO1', 'NIPAL2', 'NPHS1', 'GLIPR1', 'IL13RA1', 'NETO2', 'LAMP1', 
    'TSPAN15', 'TSPAN14', 'SLC20A2', 'PTGFRN', 'TM4SF20', 'IL3RA', 'MDGA1', 'RAET1G', 'CNTNAP3', 'SERINC5', 'CELSR3', 'IGF2R', 'MUC1', 'GPR108', 'SLC2A13', 'SGCB', 'IFNAR1', 'TSPAN17', 'PMEPA1', 'CD109', 'THBD', 'SGCD', 'FGFR2', 'GGT7', 
    'FGFRL1', 'SLC7A1', 'SCN1B', 'MFAP3L', 'SLC22A23', 'TLR4', 'TSPAN31', 'GPC1', 'PTPRJ', 'SLC16A5', 'SLC36A4', 'SLC9A7', 'PLXNA2', 'CDH19', 'MET', 'OSTM1', 'PKD2', 'ABCA1', 'SUSD1', 'ABCA7', 'SCARB2', 'MFSD12', 'LRP3', 'ADRA1D', 'SIGIRR', 
    'HM13', 'SLC39A10', 'SCAP', 'TNFRSF10D', 'IL6R', 'KIAA0319', 'LMBRD1', 'SLC46A3', 'PTTG1IP', 'CD14', 'ABCC10', 'PLXNB1', 'LRP4', 'CLDN12', 'SEZ6L2', 'SLC51B', 'FGFR4', 'SLC4A8', 'SLC39A6', 'RPN1', 'VSIG10', 'SLC12A7', 'CLDND1', 
    'JAG1', 'TMEM245', 'RNF130', 'SLC10A3', 'FZD8', 'TSPAN13', 'PIGO', 'ATP6V0A2', 'NRG4', 'TPRA1', 'MICA', 'ADCY6', 'PLAUR', 'P2RX4', 'SLC23A1', 'SLC16A4', 'IL4R', 'NPR3', 'TMCO3', 'RNF43', 'UBAC2', 'BDKRB2', 'PLXNA1', 'IFNGR1', 'TMEM87B', 'BMPR2', 'ANO6', 'SPPL2B', 'MUC16', 'SDK2', 'TCIRG1', 'LMBR1', 'CHPT1', 'GABRE', 'LDLR', 'ATP13A1', 'TSPAN33', 'SLC22A4', 'SORL1', 'SLC1A4', 'ATG9A', 'DIRC2', 'SIRPA', 'GLP2R', 'SLC22A3', 'ERBB2', 'SLCO1B3', 'C14orf132', 'ITGA5', 'SLC9A6', 'LRRC37A3', 'ASIC1', 'TMEM150A', 'EPHB4', 'TMEM8B', 'HFE', 'UNC5D', 'RNF167', 'EBP', 'C11orf24', 'GP1BA', 'LIFR', 'SCARA5', 'CSF2RA', 'SLC5A6', 'SLC29A2', 'SLC29A1', 'ROS1', 'SLC35A5', 'BMPR1A', 'M6PR', 'NCR3LG1', 'CTNS', 'TMEM37', 'TMEM219', 'NTNG2', 'MFSD11', 'IGSF3', 'TNFRSF21', 'SLCO3A1', 'SLC11A2', 'RHBDF2', 'GPRC5C', 'SLC6A9', 'BST1', 'MC1R', 'MFAP3', 'FZD6', 'EPHB6', 'SGCE', 'SLC44A2', 'CDH16', 'ROR1', 'SLC41A3', 'SUSD4', 'PIK3IP1', 'LRFN4', 'TLR6', 'ATP13A2', 'NUP210', 'MFSD5', 'CDON', 'STT3B', 'LRP8', 'TGFBR1', 'ITGA11', 'SEMA4G', 'CACNA1G', 'GDPD5', 'INSR', 'DGCR2', 'SLC9A1', 'SLC31A1', 'PCDH7', 'CXADR', 'HLA-E', 'SORT1', 'ACVR2B', 'CEACAM5', 'TSPAN9', 'ADCY3', 'TMEM30A', 'EMP3', 'STIM1', 'SLC19A2', 'MMP17', 'F2R', 'TSPAN1', 'SLC9A2', 'CNTNAP1', 'SLC36A1', 'FGFR3', 'ERBB3', 'FAM171A1', 'SCN9A', 'TSPAN5', 'CD58', 'SLITRK5', 'NOTCH3', 'CDH6', 'SLC46A1', 'CNNM3', 'CD320', 'LRFN1', 'ACP2', 'TMX3', 'KCNK5', 'CASD1', 'GFRA1', 'BTN2A1', 'CADM1', 'EPHB2', 'ZDHHC11', 'SORCS2', 'PILRB', 'TMEM184A', 'DAGLA', 'CD82', 'FKRP', 'LRFN3', 'ITGA10', 'GPR137', 'SLC41A1', 'ITGA2', 'EPHA4', 'EPHB3', 'FCGRT', 'FAM171B', 'RNF149', 'TMEM104', 'ADAM8', 'RELL1', 'IL11RA', 'CLDN15', 'KLRG1', 'EDA2R', 'LTB4R', 'SLC12A4', 'TMEM161A', 'FRRS1', 'QSOX2', 'IL17RA', 'PCDHAC1', 'P2RY2', 'JAG2', 'GPR180', 'VSIG10L', 'TMEM116', 'TCTN2', 'ADORA2B', 'PGAP1', 'NAGPA', 'ATRNL1', 'DSC2', 'LYSMD3', 'SUCO', 'SLC22A17', 'RNFT1', 'BACE2', 'HRH1', 'IGSF1', 'PANX2', 'DAGLB', 'FZD2', 'NRXN3', 'SSR1', 'FZD5', 'ITGB8', 'ERMAP', 'FAM174B', 'LPAR1', 'TMEM87A', 'PROCR', 'ADAM12', 'SYPL1', 'ADAM17', 'LRIG3', 'ECE1', 'ITGB2', 'CDH2', 'TNFRSF11A', 'IL17RC', 'PTPRG', 'BOC', 'GPC6', 'SEMA4C', 'CD37', 'LMBRD2', 'SCARF1', 'DNER', 'SERINC2', 'PTPRA', 'HTR1D', 'CHL1', 'GLRB', 'SLC2A11', 'SLC52A2', 'CACNG8', 'SLC12A9', 'CD302', 'GPRC5B', 'PLD5', 'CDH5', 'TYRO3', 'SLC6A16', 'SLC6A8', 'GABBR1', 'P2RY6', 'SLC29A3', 'ASTN2', 'IFNGR2', 'CHRNA5', 'TMEM62', 'PRTG', 'SLC26A6', 'SLC39A9', 'AXL', 'ADAM23', 'LRIG1', 'FAM189B', 'NPY4R', 'ENPP1', 'FLVCR1', 'SLC37A3', 'CNNM4', 'ZDHHC11B', 'ABCA5', 'SLC26A9', 'ZDHHC5', 'VIPR1', 'SLC15A4', 'SCNN1D', 'SLCO1A2', 'AREG', 'MMP16', 'LMAN2L', 'NRG2', 'SLC45A4', 'EREG', 'UPK3B', 'IGSF8', 'EPOR', 'IL1RAP', 'GPR141', 'VASN', 'SEMA6C', 'IL27RA', 'LRIG2', 'TP53I13', 'RAET1E', 'SLCO2B1', 'ABCB9', 'LYNX1', 'ANTXR2', 'PLXNC1', 'CACNG6', 'MST1R', 'PKD1', 'GRIN2B', 'CADM4', 'SLC38A9', 'MAMDC4', 'ACVR1', 'SLC26A2', 'ABCA12', 'NLGN4Y', 'GRIN2D', 'GPR158', 'CHRM3', 'SLC2A6', 'BCAN', 'PTPRS', 'DSC3', 'SIDT2', 'CLDN1', 'SLC8B1', 'CD22', 'CXCL16', 'SUSD5', 'SLC26A1', 'LRP11', 'TMEM182', 'C3orf80', 'PTCH1', 'SLC2A3', 'SEMA4F', 'PIEZO2', 'TMEM9B', 'PSEN1', 'RYK', 'TTYH1', 'EVC2', 'CLMP', 'RTN4RL2', 'SLC17A5', 'MXRA8', 'EMP2', 'FAS', 'SLC43A2', 'CEACAM19', 'F11R', 'CDH24', 'SLC19A1', 'LPAR6', 'ADRB2', 'BTC', 'TNFRSF10A', 'TMEM179B', 'SLC43A1', 'CNR1', 'C1orf159', 'CPM', 'SLC33A1', 'UNC5B', 'GJB2', 'DLL3', 'FAT3', 'PCDHA5', 'SLC44A5', 'FAM174A', 'SLC1A2', 'CCR7', 'ENG', 'RECK', 'RGMB', 'MCOLN1', 'FAIM2', 'CLDN4', 'PODXL2', 'SLC7A6', 'HBEGF', 'CACHD1', 'PCDHAC2', 'ITGAM', 'SMO', 'TMEM67', 'LDLRAD3', 'TM9SF1', 'GRID1', 'LYPD1', 'SLC16A6', 'GPR160', 'SLC6A15', 'EPHA5', 'ENPP4', 'LEPR', 'BTN3A1', 'ELFN2', 'LRRN2', 'GPR137B', 'S1PR5', 'BDKRB1', 'TMEM140', 'GPR155', 'GPNMB', 'SLC15A2', 'KIRREL3', 'ADCY5', 'TGFA', 'AGTR1', 'SLAMF7', 'PRRT3', 'EFNB1', 'ZFYVE27', 'CD274', 'GJB1', 'PCDHGA2', 'ANO7', 'L1CAM', 'KREMEN1', 'SLC22A15', 'F3', 'SCN8A', 'SLC13A4', 'BTN3A2', 'SLC2A12', 'RELT', 'KLRC2', 'PTPRH', 'EFNA5', 'PHEX', 'PCSK5', 'BTN2A2', 'NIPAL1', 'IL17RB', 'LMBR1L', 'ADRA2C', 'HLA-DMA', 'ITGB6', 'LRRC37A2', 'TSPAN8', 'ITGA7', 'ABCA8', 'NOTCH1', 'EGF', 'SLCO1C1', 'SHISA4', 'NPR2', 'PQLC2', 'SLC5A11', 'GRIN3B', 'IGSF9B', 'AQP2', 'TLR3', 'ITGA1', 'TMEM145', 'HYAL2', 'RTN4R', 'DPEP1', 'HLA-B', 'DISP1', 'ACKR3', 'BVES', 'EMP1', 'PCDHA4', 'PLXNB3', 'EMB', 'FAT4', 'IL18R1', 'FZD3', 'TMEM63B', 'TNFRSF25', 'GPR153', 'SEMA6A', 'LHFPL5', 'SDK1', 'FZD4', 'MUC4', 'PCDHA10', 'SLC37A4', 'TRABD2A', 'GPR161', 'PCDHA12', 'ANPEP', 'SLC37A1', 'PTPRB', 'ULBP3', 'CRB3', 'TMEM63C', 'MMP14', 'MFSD8', 'SIDT1', 'FZD1', 'P2RX6', 'GPR39', 'PCDHA11', 'PCDHGB4', 'TMEM63A', 'SLC4A4', 'SLC13A2', 'SEMA4D', 'PCDH12', 'BTNL9', 'LRRC37B', 'NALCN', 'PANX1', 'OPN3', 'SSPN', 'NPY1R', 'EDNRA', 'GPR35', 'ITGAX', 'IGSF11', 'LRRC32', 'CSPG4', 'TIMD4', 'GRM1', 'LSAMP', 'ITGAE', 'SLC2A8', 'F2RL2', 'LRRC37A', 'KLRC3', 'SLC41A2', 'TMC7', 'SUCNR1', 'P2RY1', 'MFSD2A', 'EFNB3', 'IL17RD', 'SYP', 'CLDN18', 'ENPEP', 'NCAM1', 'TMEM158', 'GGT1', 'SLC6A17', 'ANO9', 'SSTR5', 'CDHR3', 'MICB', 'CHRNB1', 'LTB4R2', 'SLC1A1', 'TMPRSS6', 'SLC38A11', 'GPM6B', 'GPR157', 'CDH4', 'TLR5', 'TMEM231', 'EFNB2', 'SCARF2', 'NCMAP', 'ACVR2A', 'SLC34A3', 'CD177', 'GRIK5', 'PCDHB2', 'SPNS2', 'GPC2', 'TMEM255B', 'PCDHGB5', 'PTGER2', 'CACNG4', 'IGDCC4', 'BACE1', 'SEMA6B', 'NLGN3', 'NIPAL4', 'FZD9', 'SLC39A4', 'GPIHBP1', 'SSTR2', 'ROBO3', 'CALHM2', 'RET', 'MMP25', 'ESYT3', 'PCDHGA1', 'UNC5C', 'VNN1', 'FLRT2', 'CHRNA7', 'MCAM', 'PCDHB6', 'ABCA4', 'TMEM132B', 'CSPG5', 'DSCAML1', 'AMIGO1', 'MEGF11', 'PCDHB14', 'CHRNB2', 'ZP3', 'TMEM26', 'TNFSF8', 'TLR1', 'HLA-F', 'ADAM19', 'ULBP1', 'GPR176', 'LPAR2', 'FAT2', 'PCDHB13', 'NFASC', 'GPR27', 'GPR132', 'PCDHGB1', 'BEST4', 'DUOX2', 'NKAIN3', 'SPRN', 'BTN3A3', 'LAG3', 'RTN4RL1', 'SLC4A5', 'FLVCR2', 'PCDHGA5', 'EPHA10', 'RNF150', 'KCNMB4', 'PSEN2', 'CD40', 'CDHR2', 'ADAM28', 'HTR1F', 'TACR2', 'TAS2R4', 'GPR179', 'IL15RA', 'IL1R2', 'GPR162', 'PILRA', 'AGER', 'SLC6A13', 'PTPRN2', 'CACNG7', 'IFNLR1', 'PCDHA3', 'EDA', 'MPZL3', 'ITGB3', 'CLEC2D', 'SLC1A7', 'PCDHGA7', 'GPR156', 'ABCB4', 'PCDHGB3', 'S1PR2', 'MANSC4', 'PTPRD', 'C5AR1', 'MPZ', 'RHCG', 'CD160', 'DUOX1', 'ITGAD', 'CDHR5', 'MUC12', 'ADORA1', 'PCDHA2', 'SLC51A', 'SLC9A3', 'ULBP2', 'ESAM', 'HRH2', 'ERVW-1', 'P2RY11', 'CLN3', 'PCDHB11', 'SLC44A3', 'SLC12A5', 'SLITRK6', 'SLC18A2', 'KIAA1324', 'ACHE', 'TENM2', 'ITLN1', 'ICOSLG', 'TLR10', 'SLC12A8', 'CLDN11', 'VTCN1', 'TNFSF4', 'NRN1', 'TNFSF13B', 'PMEL', 'ADRA1B', 'EFNA4', 'FAM171A2', 'PSCA', 'TSPAN7', 'MUC3A', 'RHBDL2', 'MERTK', 'SEMA7A', 'ADRA2B', 'THSD1', 'CRLF2', 'ISLR2', 'PCDHGC3', 'CDH26', 'OXTR', 'SLC17A7', 'CDCP1', 'SLC44A4', 'TAS2R14', 'PTGER4', 'CD163L1', 'SDC2', 'CDH23', 'CHRNB4', 'RXFP4', 'GPR75', 'PIGR', 'FNDC4', 'NOTCH4', 'TECTA', 'FLT4', 'SLC19A3', 'GPR3', 'KCNMB3', 'GRIN3A', 'MILR1', 'TMEM106A', 'FFAR4', 'SCN3B', 'TMEM178A', 'PTK7', 'CD226', 'AQP1', 'CATSPERG', 'ADAM22', 'PCDHGA6', 'SLCO1B1', 'CHRNA3', 'GRM2', 'SLC10A4', 'OR1Q1', 'TAS1R3', 'TRABD2B', 'IFNAR2', 'ASIC5', 'ADAM32', 'DLK2', 'ENTPD1', 'TAS2R10', 'KIAA1549L', 'IL12RB2', 'PCDH10', 'PEAR1', 'PLXDC2', 'IL17RE', 'EFNA3', 'SLC16A8', 'DLL1', 'GABRQ', 'TSPAN2', 'NTRK2', 'TNFRSF19', 'IL21R', 'KREMEN2', 'RNF128', 'TENM4', 'IGSF9', 'ADCY7', 'SCN5A', 'LRRC4B', 'CCR10', 'IL10RB', 'PLXNA4', 'CSF1R', 'PDCD1LG2', 'DDR2', 'PCDHB8', 'RGMA', 'SLC2A5', 'AOC3', 'CHRNA9', 'ADAM11', 'CNTNAP2', 'BST2', 'TMEM171', 'GABBR2', 'EPCAM', 'C5AR2', 'PCDHB5', 'GRIK3', 'APLNR', 'DRD4', 'SSTR1', 'SLC52A3', 'CSMD2', 'SLAMF9', 'ADRA2A', 'EFNA2', 'SGCA', 'EPHB1', 'VSTM4', 'CD96', 'HHLA2', 'CALY', 'GPR63', 'PKD2L1', 'OLR1', 'FLT3', 'TSPAN18', 'OR2I1P', 'FAP', 'ZPLD1', 'SLCO2A1', 'KDR', 'IL7R', 'SV2C', 'TREM1', 'ELFN1', 'GPR143', 'GPR50', 'SLCO5A1', 'CCKBR', 'NOX4', 'TREH', 'ITGA8', 'SLC16A12', 'TECTB', 'LPAR5', 'SLC2A14', 'GJD2', 'TMEM132E', 'PTPRT', 'FPR1'
]
expn.heatmap(
    filename="%s.pdf" % sample_name, 
    # bracket = (-2.0, 2.0), 
    heat_wid=0.25, 
    col_cluster=False, 
    heat_hei=heat_hei, 
    colbar_label="Z-score: \n"+title_sample_name, 
    row_cluster=True,
    optimal_ordering = True,
    figsize=figsize,
    highlights=highlights,
    grid=grid
)