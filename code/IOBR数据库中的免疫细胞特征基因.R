rm(list = ls())
options(stringsAsFactors = F)
setwd('/home/datahup/syj/ESCC/ssGSEA-9.1/')
library(IOBR)
library(tidyverse)

# load(file = "/home/datahup/syj/ESCC/验证集15/TG.rdata")
# my_data <- TG_dat

#查看参数和算法------
tme_deconvolution_methods
signature_score_calculation_methods
#肿瘤微环境相关特征基因集----------
names(signature_tme)

#师兄的：Macrophage Plasmacytoid_dendritic_cel

signature_tme$MHC_Class_I
signature_tme$MHC_Class_II
signature_tme$B_cells_Rooney_et_al
signature_tme$B_cells_Danaher_et_al#有
signature_tme$B_cells_Bindea_et_al
signature_tme$DC_Danaher_et_al#有，比下一个明显一点
signature_tme$DC_Bindea_et_al#有
signature_tme$DC_01_CLEC9A
signature_tme$DC_02_CD1C
signature_tme$DC_03_LAMP3
signature_tme$DC_04_CD207
signature_tme$DC_05_LILRA4
signature_tme$DC_06_STMN1
signature_tme$Macrophages_Rooney_et_al
signature_tme$Macrophages_Danaher_et_al
signature_tme$Macrophages_Bindea_et_al#有，用师兄用这个都行

# Tfh (Follicular Helper T cells): 直接参与抗原呈递，促进 B 细胞的激活和抗体生成。在淋巴结的生发中心中发挥关键作用。
# Treg (Regulatory T cells): 虽然主要功能是调节免疫反应，但它们也可以通过调节其他 T 细胞对抗原的反应来间接影响抗原呈递。
# Tn (Naive T cells): 作为未成熟 T 细胞，Tn 细胞在首次接触抗原时会被抗原呈递细胞（如树突状细胞）激活。
# Tf17 (T helper 17 cells): 参与特定的免疫反应，虽然它们主要参与炎症反应，但也可能在抗原呈递过程中起到一定作用。
signature_tme$CD4_c0_Tcm
signature_tme$CD4_c1_Treg
signature_tme$CD4_c2_Tn
signature_tme$CD4_c3_Tfh
signature_tme$CD4_c4_Tstr
signature_tme$CD4_c5_CTL
signature_tme$CD4_c6_Tn_FHIT
signature_tme$CD4_c7_Tn_TCEA3
signature_tme$CD4_c8_Tf17#you
signature_tme$CD4_c9_Tn_TCF7_SLC40A1
signature_tme$CD4_c10_Tn_LEF1_ANKRD55
signature_tme$CD4_c11_Tisg

# Teff (Effector T cells): 直接识别和杀伤被感染的细胞，依赖于抗原呈递细胞（如树突状细胞）呈递抗原。
# Tcm (Central Memory T cells): 具备记忆功能，可以在再次接触抗原时迅速激活，抗原呈递是其功能发挥的关键。
# Tn (Naive T cells): 在首次接触抗原时，需通过抗原呈递细胞的激活才能分化为效应 T 细胞。
# Trm (Tissue Resident Memory T cells): 虽然主要存在于特定组织中，但它们的激活和功能也依赖于抗原呈递。
# Teff_KLRG1 (Effector T cells expressing KLRG1): 作为效应 T 细胞的一种亚型，与抗原呈递有关。
signature_tme$CD8_c0_t_Teff
signature_tme$CD8_c1_Tex
signature_tme$CD8_c2_Teff
signature_tme$CD8_c3_Tn
signature_tme$CD8_c4_Tstr
signature_tme$CD8_c5_Tisg
signature_tme$CD8_c6_Tcm
signature_tme$CD8_c7_p_Tex
signature_tme$CD8_c8_Teff_KLRG1
signature_tme$CD8_c9_Tsen
signature_tme$CD8_c10_Teff_CD244
signature_tme$CD8_c11_Teff_SEMA4A
signature_tme$CD8_c12_Trm
signature_tme$CD8_c13_Tn_TCF7

signature_tme$CD8_Rooney_et_al
signature_tme$CD8_T_cells_Danaher_et_al
signature_tme$CD8_T_cells_Bindea_et_al#YOU
signature_tme$CD_8_T_effector


#所有免疫细胞相关的基因集（√）------
names(signature_collection)

signature_collection$CD4


























