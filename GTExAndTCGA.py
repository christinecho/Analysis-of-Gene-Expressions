__author__ = 'christine'

GTEx = ['Pancreas', 'Esophagus', 'Heart', 'Colon', 'Bone Marrow', 'Liver', 'Skin', 'Kidney', 'Brain', 'Pituitary', 'Blood', 'Thyroid', 'Uterus', 'Adipose Tissue', 'Adrenal Gland', 'Breast', 'Stomach', 'Spleen', 'Fallopian Tube', 'Lung', 'Blood Vessel', 'Ovary', 'Testis', '<not provided>', 'Muscle', 'Bladder', 'Nerve', 'Cervix Uteri', 'Prostate', 'Vagina', 'Salivary Gland', 'Small Intestine']
TCGA = ['STAD', 'KIRP', 'READ', 'HNSC', 'THCA', 'PAAD', 'MESO', 'KICH', 'ESCA', 'CESC', 'OV', 'PCPG', 'BRCA', 'CHOL', 'KIRC', 'LUAD', 'SARC', 'LGG', 'THYM', 'BLCA', 'GBM', 'SKCM', 'ACC', 'LIHC', 'TGCT', 'COAD', 'LUSC', 'UCEC', 'UVM', 'PRAD', 'UCS', 'DLBC']
print(len(TCGA))
TCGAToGTEx = {'STAD': 'Stomach', 'KIRP': 'Kidney', 'THCA': 'Thyroid', 'PAAD': 'Pancreas', 'KICH': 'Kidney',
              'ESCA': 'Esophagus', 'CESC': 'Cervix Uteri', 'OV': 'Ovary', 'PCPG': 'Adrenal Gland',
              'BRCA': 'Breast', 'KIRC': 'Kidney', 'LUAD': 'Lung', 'LGG': 'Brain', 'BLCA': 'Bladder',
              'GBM': 'Brain', 'SKCM': 'Skin', 'ACC': 'Adrenal Gland', 'LIHC': 'Liver', 'TGCT': 'Testis',
              'COAD': 'Colon', 'LUSC': 'Lung', 'UCEC': 'Uterus', 'PRAD': 'Prostate', 'UCS': 'Uterus'}


#if value is tending to the fight direction, make that as good result