from config import *
import os
import numpy as np
import pandas as pd
import collections
import pickle

def reverse_complement(seq):
    rev = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    seq = seq.strip().upper()
    return ''.join([rev[seq[i]] for i in range(len(seq) - 1, -1, -1)])


def gc_content(seq):
    counts = collections.Counter(seq)
    return (counts['G'] + counts['C']) / len(seq)


def quality_scores_to_string(quality_scores):
    return ''.join(chr(q+33) for q in quality_scores)


def quality_string_to_scores(quality_string):
    return [ord(q)-33 for q in quality_string]


def is_good_quality(quality_string, quality_thres):
    # Checks whether all quality is >= quality_thres
    quality_scores = quality_string_to_scores(quality_string)
    return np.sum(np.array(quality_scores) < quality_thres) == 0


def sanitize_postcas_unspliced(seq, match_seq=None):
    # Attempt to remove target seq suffix
    suffix = 'TGATTACAC'
    max_missing_suffix = len(suffix) - 4
    
    if match_seq:
        rev_seq = seq[::-1]
        if rev_seq.find(match_seq[-4:][::-1]) > 0:
            max_missing_suffix = len(suffix) - 1
    
    suffix_idx = seq.find(suffix)
    if suffix_idx != -1:
        return seq[0:suffix_idx]
    
    for i in range(1, max_missing_suffix + 1):
        if seq[-(len(suffix) - i):] == suffix[0:-i]:
            return seq[0:-(len(suffix) - i)]
        
    return seq

def sanitize(bc_unspliced_postcas_map):
    for bc in tqdm(bc_unspliced_postcas_map):
        unsanitized = bc_unspliced_postcas_map[bc]
        bc_unspliced_postcas_map[bc] = collections.defaultdict(int)
        for og in unsanitized:
            bc_unspliced_postcas_map[bc][sanitize_postcas_unspliced(og)] += unsanitized[og]


def hamming_distance(s1, s2):
    # assumes s1 and s2 are of equal length
    return sum(s1[i] != s2[i] for i in range(len(s1)))


def kmers(seq, k):
    return set([seq[i:i+k] for i in range(len(seq) - k)])
    
    
def dice_coefficient(seq1, seq2_kmers, k=5):    
    s1, s2 = kmers(seq1, k), seq2_kmers
    return (len(s1 & s2) * 2.0) / (len(s1) + len(s2))


def expected_dice_coefficient(target_len, product_len, k=5):
    # expected dice coefficient if target and product are truly related
    # assumption is contiguous blocks are deleted
    deletion_size = abs(target_len - product_len)
    expected_kmer_matches = target_len - k - deletion_size
    
    # multiply by 0.9 for some leeway
    return max(0, 0.9 * ((2 * expected_kmer_matches)/(target_len - k + expected_kmer_matches)))


def pickle_exists(name):
    return os.path.exists(os.path.join(PICKLE_CACHE_DIR, name + '.p'))


def save_var(var, name):
    path = os.path.join(PICKLE_CACHE_DIR, name + '.p')
    with open(path, 'wb') as handle:
        pickle.dump(var, handle)


def load_var(name):
    path = os.path.join(PICKLE_CACHE_DIR, name + '.p')
    with open(path, 'rb') as handle:
        return pickle.load(handle)


def save_bc_seq(bc_seq_map, name):
    pure_bc_seq_map = {}
    for bc in bc_seq_map:
        pure_bc_seq_map[bc] = dict(bc_seq_map[bc])
    save_var(pure_bc_seq_map, name)
    
    
def load_bc_seq(name):
    pure_bc_seq_map = load_var(name)
    bc_seq_map = collections.defaultdict(lambda: collections.defaultdict(int))
    for bc in pure_bc_seq_map:
        for seq in pure_bc_seq_map[bc]:
            bc_seq_map[bc][seq] = pure_bc_seq_map[bc][seq]
    return bc_seq_map


def save_bc_id_indel(bc_id_indel_map, name):
    pure_bc_id_indel_map = {}
    for bc in bc_id_indel_map:
        if bc not in pure_bc_id_indel_map:
            pure_bc_id_indel_map[bc] = {}
        for lid in bc_id_indel_map[bc]:
            pure_bc_id_indel_map[bc][lid] = dict(bc_id_indel_map[bc][lid])
    save_var(pure_bc_id_indel_map, name)
    
    
def load_bc_id_indel(name):
    pure_bc_id_indel_map = load_var(name)
    bc_id_indel_map = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(int)))
    for bc in pure_bc_id_indel_map:
        for lid in pure_bc_id_indel_map[bc]:
            for indel in pure_bc_id_indel_map[bc][lid]:
                bc_id_indel_map[bc][lid][indel] = pure_bc_id_indel_map[bc][lid][indel]
    return bc_id_indel_map


# ========== Indel Distribution Utils ==========
def generate_empty_distribution():
    distribution = {1: {'A': 0, 'G': 0, 'T': 0, 'C': 0}}
    for i in range(1, MAX_INDEL_LEN + 1):
        distribution[-i] = {j: 0 for j in range(i + 1)}
    return distribution

def is_distribution_empty(distribution):
    distribution_list, _ = distribution_to_list(distribution)
    return sum(distribution_list) == 0

# Returns normalized distribution that is modified inplace
def normalize_distribution(distribution):
    distribution_list, _ = distribution_to_list(distribution)
    total = float(sum(distribution_list))
    
    if total == 0:
        return distribution
    
    for indel in distribution:
        for pos in distribution[indel]:
            distribution[indel][pos] /= total
    return distribution

def distribution_to_list(distribution):
    distribution_list = [distribution[1][b] for b in 'AGTC']
    labels = ['1A', '1G', '1T', '1C']
    for i in range(1, MAX_INDEL_LEN + 1):
        for j in range(i + 1):
            distribution_list.append(distribution[-i][j])
            labels.append(str(-i) + ',' + str(j))
    return distribution_list, labels


def generate_product(seq, cutsite, deletion_size):
    prods = set()
    for genotype_pos in range(deletion_size + 1):
        a = cutsite+genotype_pos-deletion_size
        b = cutsite+genotype_pos
        if a < 0 or b > len(seq):
            continue
        simulated_product = seq[0:a] + seq[b:]
        if simulated_product not in prods:
            yield (genotype_pos, simulated_product)
        prods.add(simulated_product)
        
        
def get_cutsite(g_seq, t_seq, g_orientation):
    if g_orientation == '+':
        cutsite = t_seq.find(g_seq) + len(g_seq) - 3
    elif g_orientation == '-':
        cutsite = t_seq.find(reverse_complement(g_seq)) + 3
    return cutsite


def get_cutsites(tid):
    design_t = exp_tid_target_map[tid]
    gids = exp_tid_gids_map[tid]
    cutsites = set()
    for gid in gids:
        grna = exp_gid_grna_map[gid]
        g_orientation = exp_grna_gid_map[grna][1]
        cutsites.add(get_cutsite(grna, design_t, g_orientation))
    return cutsites


def get_simulated_product(indel, t_seq):
    if indel[1] in DELETION_SIGNATURES:
        _, _, deletion_size, genotype_pos, cutsite = indel
        a = cutsite+genotype_pos-deletion_size
        b = cutsite+genotype_pos
        return t_seq[0:a] + t_seq[b:]
    elif indel[1] in INSERTION_SIGNATURES:
        _, _, insertion_len, insertion_genotype, cutsite = indel
        return t_seq[0:cutsite] + insertion_genotype + t_seq[cutsite:]
    
    return None


def indel_iterator(indel_splice_count_map):
    for indel in sorted(indel_splice_count_map.keys()):
        yield indel, indel_splice_count_map[indel]     


# ========== Math Utils ==========
def clip(x):
    return np.clip(x, 0.00001, 0.99999)


def logit(x):
    x = clip(x)
    return np.log(x) - np.log(1 - x)


def expit(x):
    return 1. / (1. + np.exp(-x))


# ========== Load library sequences ==========
exp_design = pd.read_csv(EXP_DESIGN, index_col='Identifier number')
exp_design.sort_index(inplace=True)

exp_grna_gid_map = {}
exp_gid_grna_map = {}
exp_target_tid_map = {}
exp_tid_target_map = {}
exp_gid_tid_map = {}
exp_grna_target_map = {}
exp_gid_spliceidx_map = {}

_tid = None
for idx, row in exp_design.iterrows():
    grna_seq = row['Designed gRNA (NGG orientation, 19 and 20)']
    grna_orientation = row['gRNA orientation']
    target_seq = row['Designed 61-bp target site (37i-24e, AG)']
    sequence = row['Sequence']
    sequence_spliceidx = row['Sequence position of AG (0-index)']
    
    exp_grna_gid_map[grna_seq] = (idx, grna_orientation)
    exp_gid_grna_map[idx] = grna_seq
    exp_grna_target_map[grna_seq] = target_seq
    
    exp_gid_spliceidx_map[idx] = sequence_spliceidx - sequence.find(target_seq) + 2
    
    if target_seq not in exp_target_tid_map:
        exp_target_tid_map[target_seq] = idx
        exp_tid_target_map[idx] = target_seq
        _tid = idx
        
    exp_gid_tid_map[idx] = _tid
    
exp_tid_gids_map = collections.defaultdict(set)
for gid in exp_gid_tid_map:
    tid = exp_gid_tid_map[gid]
    exp_tid_gids_map[tid].add(gid)
    

def gid_to_gt(gid):
    tid = exp_gid_tid_map[gid]
    design_t = exp_tid_target_map[tid]
    grna = exp_gid_grna_map[gid]
    g_orientation = exp_grna_gid_map[grna][1]
    pair = (grna, design_t, g_orientation)
    return pair