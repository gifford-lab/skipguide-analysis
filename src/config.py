import os

# Number of cores when multiprocessing
NUM_PROCESSES = 60

# Root Directory
PRJ_DIR = os.path.abspath(os.path.join(os.getcwd(), os.pardir))

# Data Paths
DATA_DIR = os.path.join(PRJ_DIR, 'data')
READS_DIR = os.path.join(DATA_DIR, 'reads')
READ_FN_PREFIX = '190510Gif_D19-212'

LIB_DIR = os.path.join(DATA_DIR, 'libsa')
EXP_DESIGN = os.path.join(LIB_DIR, '061318_exonskipping_library.csv')

# Output Figures/Tables Paths
OUTPUT_DIR = os.path.join(PRJ_DIR, 'output')
IMAGES_DIR = os.path.join(OUTPUT_DIR, 'images')
TABLES_DIR = os.path.join(OUTPUT_DIR, 'tables')

# Tools Paths
INDELPHI_DIR = os.path.join(PRJ_DIR, 'tools/inDelphi/inDelphi-model/')
INDELPHI_CELL_TYPE = 'mESC'
MAXENTSCAN_DIR = os.path.join(PRJ_DIR, 'tools/MaxEntScan/fordownload/')

# Cache Paths
CACHE_DIR = os.path.join(PRJ_DIR, 'cache')
READS_CACHE_DIR = os.path.join(CACHE_DIR, 'reads')
PICKLE_CACHE_DIR = os.path.join(CACHE_DIR, 'pickles')
TEMP_DIR = os.path.join(CACHE_DIR, 'temp')

# MESC
PRECAS_UNSPLICED_SUFFIX = '020'
POSTCAS_UNSPLICED1_SUFFIX = '021'
POSTCAS_UNSPLICED2_SUFFIX = '022'
PRECAS_SPLICED_SUFFIX = '023'
POSTCAS_SPLICED1_SUFFIX = '024'
POSTCAS_SPLICED2_SUFFIX = '025'
PRECAS_GDNA_SUFFIX = '026'
POSTCAS_GDNA1_SUFFIX = '027'
POSTCAS_GDNA2_SUFFIX = '028'

PRECAS_GDNA_BC_SEQ = os.path.join(READS_CACHE_DIR, READ_FN_PREFIX + PRECAS_GDNA_SUFFIX + '_BC_seq.csv')
POSTCAS_GDNA1_BC_SEQ = os.path.join(READS_CACHE_DIR, READ_FN_PREFIX + POSTCAS_GDNA1_SUFFIX + '_BC_seq.csv')
POSTCAS_GDNA2_BC_SEQ = os.path.join(READS_CACHE_DIR, READ_FN_PREFIX + POSTCAS_GDNA2_SUFFIX + '_BC_seq.csv')
PRECAS_UNSPLICED_BC_SEQ = os.path.join(READS_CACHE_DIR, READ_FN_PREFIX + PRECAS_UNSPLICED_SUFFIX + '_BC_seq.csv')
POSTCAS_UNSPLICED1_BC_SEQ = os.path.join(READS_CACHE_DIR, READ_FN_PREFIX + POSTCAS_UNSPLICED1_SUFFIX + '_BC_seq.csv')
POSTCAS_UNSPLICED2_BC_SEQ = os.path.join(READS_CACHE_DIR, READ_FN_PREFIX + POSTCAS_UNSPLICED2_SUFFIX + '_BC_seq.csv')
PRECAS_SPLICED_BC_SEQ = os.path.join(READS_CACHE_DIR, READ_FN_PREFIX + PRECAS_SPLICED_SUFFIX + '_BC_seq.csv')
POSTCAS_SPLICED1_BC_SEQ = os.path.join(READS_CACHE_DIR, READ_FN_PREFIX + POSTCAS_SPLICED1_SUFFIX + '_BC_seq.csv')
POSTCAS_SPLICED2_BC_SEQ = os.path.join(READS_CACHE_DIR, READ_FN_PREFIX + POSTCAS_SPLICED2_SUFFIX + '_BC_seq.csv')

# Regular expression config
MAX_GDNA_R1_SUBSTITUTIONS = 5
MAX_GDNA_R2_SUBSTITUTIONS = 5
MAX_UNSPLICED_R1_SUBSTITUTIONS = 3
MAX_UNSPLICED_R2_SUBSTITUTIONS = 5
MAX_SPLICED_R1_SUBSTITUTIONS = 5
MAX_SPLICED_SKP_R1_SUBSTITUTIONS = 3
MAX_SPLICED_BC_R1_SUBSTITUTIONS = 5
MAX_SPLICED_R2_SUBSTITUTIONS = 5

# Reads Quality Threshold
BC_QUAL_THRES = 11
SEQ_QUAL_THRES = 11

MAX_INDEL_LEN = 60
DELETION_SIGNATURES = ['DS']
INSERTION_SIGNATURES = ['IS']

# ============= Cached Variable Names ===============
BC_GRNA_PRECAS_MAP = 'bc_grna_precas_map'
BC_GRNA1_POSTCAS_MAP = 'bc_grna1_postcas_map'
BC_GRNA2_POSTCAS_MAP = 'bc_grna2_postcas_map'

BC_UNSPLICED_PRECAS_MAP = 'bc_unspliced_precas_map'
BC_UNSPLICED1_POSTCAS_MAP = 'bc_unspliced1_postcas_map'
BC_UNSPLICED2_POSTCAS_MAP = 'bc_unspliced2_postcas_map'

BC_SPLICED_PRECAS_MAP = 'bc_spliced_precas_map'
BC_SPLICED1_POSTCAS_MAP = 'bc_spliced1_postcas_map'
BC_SPLICED2_POSTCAS_MAP = 'bc_spliced2_postcas_map'



BC_GID_PRECAS_MAP = 'bc_gid_precas_map'
BC_GID1_POSTCAS_MAP = 'bc_gid1_postcas_map'
BC_GID2_POSTCAS_MAP = 'bc_gid2_postcas_map'

BC_TID_UNSPLICED_PRECAS_MAP = 'bc_tid_unspliced_precas_map'
BC_GID_UNSPLICED_PRECAS_MAP = 'bc_gid_unspliced_precas_map'

BC_TID_INDEL_POSTCAS_MAP = 'bc_tid_indel_postcas_map'
BC_GID_INDEL_POSTCAS_MAP = 'bc_gid_indel_postcas_map'

UNNORM_GT_INDEL_DIST_MAP = 'unnorm_gt_indel_dist_map'
GT_INDEL_DIST_MAP = 'gt_indel_dist_map'

PREDICTED_GT_INDEL_DIST_MAP = 'predicted_gt_indel_dist_map'


BC_TID_SPLICEIDX_PRECAS_MAP = 'bc_tid_spliceidx_precas_map'
BC_TID_SPLICEIDX_POSTCAS_MAP = 'bc_tid_spliceidx_postcas_map'

INDEL_SPLICE_PRECAS_COUNT_MAP = 'indel_splice_precas_count_map'
INDEL_SPLICE_POSTCAS_COUNT_MAP = 'indel_splice_postcas_count_map'
GT_SPLICE_COUNT_MAP = 'gt_splice_count_map'
GT_PRECAS_SPLICE_COUNT_MAP = 'gt_precas_splice_count_map'

DAT_A_INDELPHI_SEQUENCE_IDENTITY = 'datA_indelphi_sequence_identity'

SPLICEAI_OUTPUT_DIR = os.path.join(CACHE_DIR, 'spliceai_out')
SPLICEAI_PRECAS_PREDS_PATH = os.path.join(SPLICEAI_OUTPUT_DIR, 'spliceai_precas_preds.npy')
SPLICEAI_POSTCAS_PREDS_PATH = os.path.join(SPLICEAI_OUTPUT_DIR, 'spliceai_postcas_preds.npy')
SPLICEAI_GT_REPAIR_PREDS = os.path.join(SPLICEAI_OUTPUT_DIR, 'spliceai_gt_indels') 

MMSPLICE_OUTPUT_DIR = os.path.join(CACHE_DIR, 'mmsplice_out')
MMSPLICE_PRECAS_DF_PATH = os.path.join(MMSPLICE_OUTPUT_DIR, 'mmsplice_precas_df.csv')
MMSPLICE_POSTCAS_DF_PATH = os.path.join(MMSPLICE_OUTPUT_DIR, 'mmsplice_postcas_df.csv')
MMSPLICE_GT_DF_DIR = os.path.join(MMSPLICE_OUTPUT_DIR, 'mmsplice_gt_indels')
MMSPLICE_GT_REPAIR_PREDS = os.path.join(MMSPLICE_OUTPUT_DIR, 'mmsplice_gt_indels') 