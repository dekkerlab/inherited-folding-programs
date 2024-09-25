# on limbo ...
_bw_path = "/abyss/sergpolly/data_ranger/bigwigs"
_clr_path = "/abyss/sergpolly/data_ranger/finalcoolers"
_mega_path = "/abyss/sergpolly/data_ranger/finalcoolers_mega"
_pubclr_path = "/abyss/sergpolly/data_ranger/publiccoolers"
# # on dell17
# _path = "/data/proj_sync"

bws = {
    # #
    # # our precious atac ! danke !
    # "mM.atac": f"{_bw_path}/mMp.reducedchroms.sorted.p4m5.sorted.genomecov.sorted.bw",
    # # "mT.atac": f"{_bw_path}/mT.mRp.clN.bigWig",
    # "mG.atac": f"{_bw_path}/mGp.reducedchroms.sorted.p4m5.sorted.genomecov.sorted.bw",
    # "pM.atac": f"{_bw_path}/pMp.reducedchroms.sorted.p4m5.sorted.genomecov.sorted.bw",
    # # "pT.atac": f"{_bw_path}/pT.mRp.clN.bigWig",
    # "pG.atac": f"{_bw_path}/pGp.reducedchroms.sorted.p4m5.sorted.genomecov.sorted.bw",
    # #
    # our precious atac ! danke !
    "mM.atac": f"{_bw_path}/Mp.genomcov.rawcountsscale.bw",  # mM and pM merged together ...
    "mG.atac": f"{_bw_path}/mGp.genomcov.rawcountsscale.bw",
    "pG.atac": f"{_bw_path}/pGp.genomcov.rawcountsscale.bw",
    #
    # relevant public tracks for DLD1 from Nezar/George
    "H3K27me3_new": f"{_bw_path}/GSE251932_DLD1.NT.H3K27me3_VS_Input.fc.signal.bw",
    "H3K9me3": f"{_bw_path}/GSE251932_DLD1.NT.H3K9me3_VS_Input.fc.signal.bw",
    # relevant public tracks for DLD1 (others)
    "H3K27me3": f"{_bw_path}/GSM6597768_CB005_H3K27me3_N_DLD-1_parental_500K_KMET_S5_uniq.b20s100dt.bw",
    "H3K4me3": f"{_bw_path}/GSM6597766_CB003_H3K4me3_N_DLD-1_parental_500K_KMET_S3_uniq.b20s100dt.bw",
    "H3K27ac": f"{_bw_path}/GSM6597770_CB007_H3K27ac_N_DLD-1_parental_500K_S7_uniq.b20s100dt.bw",
    "H3K36me3": f"{_bw_path}/H3K36me3_DRX013178.bw",
    "H3K4me1": f"{_bw_path}/H3K4me1_DRX013183.bw",
    "MED1": f"{_bw_path}/MED1_SRX11555394.bw",
    "ctcf": f"{_bw_path}/GSM4238557_CTCF_DLD1_19036958_19036956_hg38_noFr.bw",
    #
    # some rna-seq - Async rangap control and depletion ...
    "rg0r1fwd":f"{_bw_path}/rnaseq/rg0r1.forward.bigWig",
    "rg0r1rev":f"{_bw_path}/rnaseq/rg0r1.reverse.bigWig",
    "rg8r1fwd":f"{_bw_path}/rnaseq/rg8r1.forward.bigWig",
    "rg8r1rev":f"{_bw_path}/rnaseq/rg8r1.reverse.bigWig",
}


telo_reps_dict = {
    # mMito
    'mMito_R1': f'{_clr_path}/ranGAP1-0-Ms-R1.hg38.mapq_30.1000.mcool',
    'mMito_R2': f'{_clr_path}/ranGAP1-0-Ms-R2.hg38.mapq_30.1000.mcool',
    # Telo/Cyto 1 replicate each
    'mTelo': f'{_clr_path}/ranGAP1-0-Telo-R2.hg38.mapq_30.1000.mcool',
    'mCyto': f'{_clr_path}/ranGAP1-0-Cyto-R1.hg38.mapq_30.1000.mcool',
    # G1 m5hrs
    'm5h_R1': f'{_clr_path}/ranGAP1-0-G1s-R1.hg38.mapq_30.1000.mcool',
    'm5h_R2': f'{_clr_path}/ranGAP1-0-G1s-R2.hg38.mapq_30.1000.mcool',
    # G1 10 hrs
    'm10h_R1': f'{_clr_path}/ranGAP1-0-10hR-R1.hg38.mapq_30.1000.mcool',
    'm10h_R2': f'{_clr_path}/ranGAP1-0-10hR-R3.hg38.mapq_30.1000.mcool',
    # pMito
    'pMito_R1': f'{_clr_path}/ranGAP1-2-Ms-R1.hg38.mapq_30.1000.mcool',
    'pMito_R2': f'{_clr_path}/ranGAP1-2-Ms-R2.hg38.mapq_30.1000.mcool',
    # pTelo/Cyto 1 replicate each
    'pTelo': f'{_clr_path}/ranGAP1-4-Telo-R2.hg38.mapq_30.1000.mcool',
    'pCyto': f'{_clr_path}/ranGAP1-4-Cyto-R1.hg38.mapq_30.1000.mcool',
    # G1 p5hrs
    'p5h_R1': f'{_clr_path}/ranGAP1-7-G1s-R1.hg38.mapq_30.1000.mcool',
    'p5h_R2': f'{_clr_path}/ranGAP1-7-G1s-R2.hg38.mapq_30.1000.mcool',
    # G1 p10hrs
    'p10h_R1': f'{_clr_path}/ranGAP1-12-10hR-R1.hg38.mapq_30.1000.mcool',
    'p10h_R2': f'{_clr_path}/ranGAP1-12-10hR-R3.hg38.mapq_30.1000.mcool',
    # mp 10 hrs
    'mp10h_R1': f'{_clr_path}/ranGAP1-7-10hR-R1.hg38.mapq_30.1000.mcool',
    'mp10h_R2': f'{_clr_path}/ranGAP1-7-10hR-R3.hg38.mapq_30.1000.mcool',
    # Nup93 m5 hrs G1
    'N93m5_R1': f'{_clr_path}/N93-0-G1s-R1.hg38.mapq_30.1000.mcool',
    'N93m5_R2': f'{_clr_path}/N93-0-G1s-R2.hg38.mapq_30.1000.mcool',
    # Nup93 m10 hrs G1
    'N93m10_R1': f'{_clr_path}/N93-0-10hR-G1s-R3.hg38.mapq_30.1000.mcool',
    'N93m10_R2': f'{_clr_path}/N93-0-10hR-G1s-R4.hg38.mapq_30.1000.mcool',
    # Nup93 p5 hrs G1
    'N93p5_R1': f'{_clr_path}/N93-7-G1s-R1.hg38.mapq_30.1000.mcool',
    'N93p5_R2': f'{_clr_path}/N93-7-G1s-R2.hg38.mapq_30.1000.mcool',
    # Nup93 p10 hrs G1
    'N93p10_R1': f'{_clr_path}/N93-12-10hR-G1s-R3.hg38.mapq_30.1000.mcool',
    'N93p10_R2': f'{_clr_path}/N93-12-10hR-G1s-R4.hg38.mapq_30.1000.mcool',
    # Nup93 mp - we're not really using it
    'N93mp10_R1': f'{_clr_path}/N93-7-10hR-G1s-R3.hg38.mapq_30.1000.mcool',
    'N93mp10_R2': f'{_clr_path}/N93-7-10hR-G1s-R4.hg38.mapq_30.1000.mcool',
}


telo_dict = {
    "mMito": f"{_clr_path}/ranGAP1-0-Ms-R1R2.hg38.mapq_30.1000.mcool",
    "mTelo": f"{_clr_path}/ranGAP1-0-Telo-R2.hg38.mapq_30.1000.mcool",
    "mCyto": f"{_clr_path}/ranGAP1-0-Cyto-R1.hg38.mapq_30.1000.mcool",
    # "m5hR1": f"{_clr_path}/ASDLDranGAP1-0-5hHiCD3-R1__hg38.hg38.mapq_30.1000.mcool",
    # "m5hR2": f"{_clr_path}/ASDLDranGAP1-0-5hHiCD3-R2__hg38.hg38.mapq_30.1000.mcool",
    "m5hR1R2": f"{_clr_path}/ranGAP1-0-G1s-R1R2.hg38.mapq_30.1000.mcool",
    # "m10hR1": f"{_clr_path}/ASDLDranGAP1-0-10hHiCD3-R1__hg38.hg38.mapq_30.1000.mcool",
    # "m10hR2": f"{_clr_path}/ASDLDranGAP1-0-10hHiCD3-R2__hg38.hg38.mapq_30.1000.mcool",
    "m10hR1R2": f"{_clr_path}/ranGAP1-0-10hR-R1R3.hg38.mapq_30.1000.mcool",
    # ...
    "pMito": f"{_clr_path}/ranGAP1-2-Ms-R1R2.hg38.mapq_30.1000.mcool",  # 4X as deep
    "pTelo": f"{_clr_path}/ranGAP1-4-Telo-R2.hg38.mapq_30.1000.mcool",
    "pCyto": f"{_clr_path}/ranGAP1-4-Cyto-R1.hg38.mapq_30.1000.mcool",
    # "p5hR1": f"{_clr_path}/ASDLDranGAP1-7-5hHiCD3-R1__hg38.hg38.mapq_30.1000.mcool",
    # "p5hR2": f"{_clr_path}/ASDLDranGAP1-7-5hHiCD3-R2__hg38.hg38.mapq_30.1000.mcool",
    "p5hR1R2": f"{_clr_path}/ranGAP1-7-G1s-R1R2.hg38.mapq_30.1000.mcool",
    # "p10hR1": f"{_clr_path}/ASDLDranGAP1-7-10hHiCD3-R1__hg38.hg38.mapq_30.1000.mcool",
    # "p10hR2": f"{_clr_path}/ASDLDranGAP1-7-10hHiCD3-R2__hg38.hg38.mapq_30.1000.mcool",
    "p10hR1R2": f"{_clr_path}/ranGAP1-12-10hR-R1R3.hg38.mapq_30.1000.mcool",
    # special one! - MP
    # ...
    "mp10hR1R2": f"{_clr_path}/ranGAP1-7-10hR-R1R3.hg38.mapq_30.1000.mcool",
    #
    # Nup93 samples ...
    "N93m5": f"{_clr_path}/N93-0-G1s-R1R2.hg38.mapq_30.1000.mcool",
    "N93m10": f"{_clr_path}/N93-0-10hR-G1s-R3R4.hg38.mapq_30.1000.mcool",
    "N93p5": f"{_clr_path}/N93-7-G1s-R1R2.hg38.mapq_30.1000.mcool",
    "N93p10": f"{_clr_path}/N93-12-10hR-G1s-R3R4.hg38.mapq_30.1000.mcool",
    "N93mp10": f"{_clr_path}/N93-7-10hR-G1s-R3R4.hg38.mapq_30.1000.mcool",
}


pubclr_dict = {
    "dldmicroc": f"{_pubclr_path}/GSM5394172_control_mergeRep.mcool",
}


mega_telo_dict = {
    # "N93-0-10hR-G1s-R3R4" : f"{_mega_path}/N93-0-10hR-G1s-R3R4.hg38.mapq_30.1000.mcool",
    # "N93-0-G1s-R1R2" : f"{_mega_path}/N93-0-G1s-R1R2.hg38.mapq_30.1000.mcool",
    # "N93-12-10hR-G1s-R3R4" : f"{_mega_path}/N93-12-10hR-G1s-R3R4.hg38.mapq_30.1000.mcool",
    # "N93-7-10hR-G1s-R3R4" : f"{_mega_path}/N93-7-10hR-G1s-R3R4.hg38.mapq_30.1000.mcool",
    # "N93-7-G1s-R1R2" : f"{_mega_path}/N93-7-G1s-R1R2.hg38.mapq_30.1000.mcool",
    "N93pG1s_MEGA" : f"{_mega_path}/N93-aux-G1s-MEGA.hg38.mapq_30.1000.mcool",
    "N93mG1s_MEGA" : f"{_mega_path}/N93-noaux-G1s-MEGA.hg38.mapq_30.1000.mcool",
    # "ranGAP1-0-10hR-R1R3" : f"{_mega_path}/ranGAP1-0-10hR-R1R3.hg38.mapq_30.1000.mcool",
    # "ranGAP1-0-G1s-R1R2" : f"{_mega_path}/ranGAP1-0-G1s-R1R2.hg38.mapq_30.1000.mcool",
    # "ranGAP1-0-Ms-R1R2" : f"{_mega_path}/ranGAP1-0-Ms-R1R2.hg38.mapq_30.1000.mcool",
    # "ranGAP1-12-10hR-R1R3" : f"{_mega_path}/ranGAP1-12-10hR-R1R3.hg38.mapq_30.1000.mcool",
    # "ranGAP1-2-Ms-R1R2" : f"{_mega_path}/ranGAP1-2-Ms-R1R2.hg38.mapq_30.1000.mcool",
    # "ranGAP1-7-10hR-R1R3" : f"{_mega_path}/ranGAP1-7-10hR-R1R3.hg38.mapq_30.1000.mcool",
    # "ranGAP1-7-G1s-R1R2" : f"{_mega_path}/ranGAP1-7-G1s-R1R2.hg38.mapq_30.1000.mcool",
    "pG1s_MEGA" : f"{_mega_path}/ranGAP1-aux-G1s-MEGA.hg38.mapq_30.1000.mcool",
    "Ms_MEGA" : f"{_mega_path}/ranGAP1-Ms-MEGA.hg38.mapq_30.1000.mcool",
    "mG1s_MEGA" : f"{_mega_path}/ranGAP1-noaux-G1s-MEGA.hg38.mapq_30.1000.mcool",
}


# bws_vlim = {
#     "mM.atac": dict(vmin=0.2,vmax=4.),
#     # "mT.atac": dict(vmin=0.005,vmax=0.3),
#     "mG.atac": dict(vmin=0.2,vmax=4.),
#     # "pM.atac": dict(vmin=0.01,vmax=0.3),
#     # # "pT.atac": dict(vmin=0.005,vmax=0.3),
#     "pG.atac": dict(vmin=0.2,vmax=4.),
#     "H3K36me3": dict(vmin=0.03,vmax=0.08),
#     "H3K4me1": dict(vmin=0.04,vmax=0.1),
#     "MED1": dict(vmin=0.04,vmax=0.15),
#     "H3K27me3": dict(vmin=0.6,vmax=1),
#     "H3K4me3": dict(vmin=1,vmax=18),
#     "H3K27ac": dict(vmin=1,vmax=16),
#     "ctcf": dict(vmin=2,vmax=18),
#     # some rna-seq - Async rangap control and depletion ...
#     "rg0r1fwd": dict(vmin=1,vmax=13),
#     "rg0r1rev": dict(vmin=1,vmax=13),
#     "rg8r1fwd": dict(vmin=1,vmax=13),
#     "rg8r1rev": dict(vmin=1,vmax=13),

# }


bws_vlim = {
    "mM.atac": dict(vmin=0.005,vmax=0.3),
    # "mT.atac": dict(vmin=0.005,vmax=0.3),
    "mG.atac": dict(vmin=0.005,vmax=0.3),
    # "pM.atac": dict(vmin=0.005,vmax=0.3),
    # # "pT.atac": dict(vmin=0.005,vmax=0.3),
    "pG.atac": dict(vmin=0.005,vmax=0.3),
    "H3K36me3": dict(vmin=0.03,vmax=0.08),
    "H3K4me1": dict(vmin=0.04,vmax=0.1),
    "MED1": dict(vmin=0.04,vmax=0.15),
    "H3K27me3": dict(vmin=0.2,vmax=1),
    "H3K4me3": dict(vmin=0.5,vmax=12),
    "H3K27ac": dict(vmin=0.5,vmax=15),
    "ctcf": dict(vmin=2,vmax=5),
    # some rna-seq - Async rangap control and depletion ...
    "rg0r1fwd": dict(vmin=0.5,vmax=12),
    "rg0r1rev": dict(vmin=0.5,vmax=12),
    "rg8r1fwd": dict(vmin=0.5,vmax=12),
    "rg8r1rev": dict(vmin=0.5,vmax=12),

}


bws_atac_vlim = {
    "mM.atac": dict(vmin=0.2,vmax=4.),
    # "mT.atac": dict(vmin=0.005,vmax=0.3),
    "mG.atac": dict(vmin=0.2,vmax=4.),
    # "pM.atac": dict(vmin=0.01,vmax=0.3),
    # # "pT.atac": dict(vmin=0.005,vmax=0.3),
    "pG.atac": dict(vmin=0.2,vmax=4.),
    "H3K36me3": dict(vmin=0.03,vmax=0.08),
    "H3K4me1": dict(vmin=0.04,vmax=0.1),
    "MED1": dict(vmin=0.04,vmax=0.15),
    "H3K27me3": dict(vmin=0.2,vmax=1),
    "H3K4me3": dict(vmin=0.5,vmax=12),
    "H3K27ac": dict(vmin=0.5,vmax=15),
    "ctcf": dict(vmin=2,vmax=5),
    # some rna-seq - Async rangap control and depletion ...
    "rg0r1fwd": dict(vmin=0.5,vmax=12),
    "rg0r1rev": dict(vmin=0.5,vmax=12),
    "rg8r1fwd": dict(vmin=0.5,vmax=12),
    "rg8r1rev": dict(vmin=0.5,vmax=12),

}
