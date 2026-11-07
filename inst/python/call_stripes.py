# Script : call_stripes.py
#   This script contains functions to identify stripe anchors and length ranges
#   from signal probability matrix using a Bayesian NDLM,
#   as well as a function to load Hi-C matrices and incorporate interaction counts.


import os
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
import subprocess
import sys 
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from scipy.stats import norm
import rpy2.robjects as robjects
from rpy2.robjects import numpy2ri
import importlib.resources
import csv 
import gc 

# hicstraw install
try:
    import hictkpy
except ImportError:
    #print("hicstraw installing...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "hictkpy"])
    import hictkpy

def add_hic_counts_to_features(feature_file, hic_file, chrom, resolution):

    # get hic interaction counts
    normalization = 'NONE'  
    chr_num = chrom[3:]
    print(f"[{chrom}] {resolution} resolution Raw Matrix loading...")
    f = hictkpy.File(hic_file, resolution)
    sel = f.fetch(chr_num, chr_num, join=True, normalization=normalization)
    df = sel.to_df()
    df['x_idx'] = df['start1'] // resolution
    df['y_idx'] = df['start2'] // resolution

    # get pair features
    pair_df = pd.read_csv(feature_file)

    count_df = df[['x_idx', 'y_idx', 'count']].copy()
    count_df.columns = ['idx1', 'idx2', 'count']

    merged_df = pd.merge(
        pair_df, 
        count_df, 
        on=['idx1', 'idx2'], 
        how='left'
    )
    merged_df['count'] = merged_df['count'].fillna(0)

    # save file 
    feature_path = os.path.dirname(feature_file)
    new_file_name = chrom + '_pair_features.csv'
    new_file_path = os.path.join(feature_path, new_file_name)
    merged_df.to_csv(new_file_path, index=False)
    #print(f"GC, Map, Count, Distance features saved in : {new_file_path}")

    if os.path.exists(feature_file):
        os.remove(feature_file)
    
    del f, sel, df, count_df, merged_df
    gc.collect()


def get_diag_matrix(feature_file, max_d):

    df = pd.read_csv(feature_file)

    df['idx1'] = df['idx1'].astype(int)
    df['idx2'] = df['idx2'].astype(int)

    # get matrix size 
    max_matrix_idx = int(max(df['idx1'].max(), df['idx2'].max()))

    # fill matrix with signal probability 
    df_flipped = df.rename(columns={'idx1': 'idx2', 'idx2': 'idx1'})
    df_symmetric = pd.concat([df, df_flipped]).drop_duplicates(subset=['idx1', 'idx2'])
    df_pivot = df_symmetric.pivot(index='idx1', columns='idx2', values='signal_p')
    full_range = np.arange(0, max_matrix_idx + 1)
    final_matrix_df = df_pivot.reindex(index=full_range, columns=full_range)

    signal_matrix = final_matrix_df.to_numpy()
    signal_matrix_filtered = signal_matrix.copy()

    # apply zero to signal_p for interactions which have signal probability less than 0.5
    signal_matrix_filtered[signal_matrix_filtered < 0.5] = 0.0

    n_rows = signal_matrix_filtered.shape[0]
    diag_to_forward_mat = np.full((n_rows, max_d + 1), np.nan)
    diag_to_reverse_mat = np.full((n_rows, max_d + 1), np.nan)

    for k in range(max_d + 1):
        diag_vals_fw = np.diagonal(signal_matrix_filtered, offset=k)
        diag_vals_rv = np.diagonal(signal_matrix_filtered, offset=-k)
        len_fw = len(diag_vals_fw)
        len_rv = len(diag_vals_rv)
        diag_to_forward_mat[:len_fw, k] = diag_vals_fw
        diag_to_reverse_mat[n_rows-len_rv:, k] = diag_vals_rv

    del signal_matrix, signal_matrix_filtered, final_matrix_df, df_pivot
    gc.collect()

    return diag_to_forward_mat, diag_to_reverse_mat


def get_stripe_direction_matrix(diag_to_forward_mat, diag_to_reverse_mat, resolution=1e4, window_size=51*1e4):
    
    target_shape = diag_to_forward_mat.shape
    n_row = target_shape[0]
    n_col = target_shape[1]
    max_d = n_col -1 

    # window size (left_right/up_down +-w bins)
    window_size_bin = int(window_size / resolution)
    w = window_size_bin // 2 if window_size_bin % 2 == 0 else (window_size_bin - 1) // 2

    # get target index for analyzing 
    row_idx, col_idx = np.indices(target_shape)

    # get stripe direction matrix by comparing means of h-direction window and v-direction window
    direction_matrix_forward = np.full(target_shape, np.nan, dtype=float)
    direction_matrix_reverse = np.full(target_shape, np.nan, dtype=float)
    rows, cols = np.where(~np.isnan(diag_to_forward_mat))
    for i, j in zip(rows, cols):

        # get window indices
        w_indices_fw_st = max(i+j -w, i) - i 
        w_indices_fw_end = min(i+j +w, i+max_d, n_row-1) - i
        w_indices_fw = np.arange(w_indices_fw_st, w_indices_fw_end+1)

        w_indices_rv_end = (i+j) - max(i - w, i+j-max_d, 0)
        w_indices_rv_st = (i+j) - min(i + w, i+j)
        w_indices_rv = np.arange(w_indices_rv_st, w_indices_rv_end+1)

        # get left-right +-w signal probabilities of (i,j)
        valid_h_values = np.array(diag_to_forward_mat[i, w_indices_fw])

        # get up-down +-w signal probabilities of (i,j)
        valid_v_values = np.array(diag_to_reverse_mat[i+j, w_indices_rv])

        # get mean signal probability of left-right window
        sum_h = np.nansum(valid_h_values)
        count_h = np.sum(~np.isnan(valid_h_values))
        h_mean = (sum_h / count_h) if count_h > 0 else 0.0

        # get mean signal probability of up-down window
        sum_v = np.nansum(valid_v_values)
        count_v = np.sum(~np.isnan(valid_v_values))
        v_mean = (sum_v / count_v) if count_v > 0 else 0.0
        
        # get stripe direction by comparing signal probability means
        if h_mean > v_mean:
            direction_matrix_forward[i, j] = 1.0  
            direction_matrix_reverse[i+j, j] = 1.0  
        elif h_mean < v_mean:
            direction_matrix_forward[i, j] = 2.0  
            direction_matrix_reverse[i+j, j] = 2.0  
        else:
            direction_matrix_forward[i, j] = 3.0  
            direction_matrix_reverse[i+j, j] = 3.0  

                
    #print("Stripe direction matrix successfully built.")
    #print(
    #    f"Classified interactions -> Direction 1: {np.sum(direction_matrix_forward == 1):,}, Direction 2: {np.sum(direction_matrix_forward == 2):,}, Direction 3: {np.sum(direction_matrix_forward == 3):,}"
    #)
    
    return direction_matrix_forward, direction_matrix_reverse


def calculate_row_sums(diag_matrix, direction_matrix, direction):

    # calculate rowsum to call stripes
    # direction=1 for 5' stripes
    # direction=2 for 3' stripes
    cond_dir = (direction_matrix == direction) | (direction_matrix == 3)
    filtered_diag = np.where(cond_dir, diag_matrix, np.nan)
    row_sum_vector = np.nansum(filtered_diag, axis=1)
    
    return row_sum_vector


def get_stripe_candidates(rowsum_smoothed_mean, high_signal_idx, peak_prominence=0.1):
    
    peaks_idx, properties = find_peaks(rowsum_smoothed_mean, prominence=peak_prominence)

    # idx satisfying both peaks and high-signal
    common_candidates = set(peaks_idx).intersection(set(high_signal_idx))
    candidates_idx = sorted(list(common_candidates))

    return candidates_idx

def get_stripe_width_range(stripe_candidates, rowsum, rowsum_smoothed_mean, rowsum_smoothed_var, p_value, threshold_mean_to_peak=0.1, resolution=1e4, max_width=3e5):

    n = len(rowsum_smoothed_mean)
    rowsum_smoothed_std = np.sqrt(rowsum_smoothed_var)
    results = []
    max_bin_width = int(max_width // resolution)
    max_half_width = max_bin_width // 2 
        
    for p in stripe_candidates:

        peak_val = rowsum[p]
        p_val = p_value[p]
        m = rowsum_smoothed_mean[p]
        std = rowsum_smoothed_std[p]

        # setting cut_value to decide width
        # probability mass from mean to cut_value = (probability mass from mean to peak) * threshold
        target_cdf = 0.5 + (0.5 - p_val) * threshold_mean_to_peak
        cut_value = norm.ppf(target_cdf, loc=m , scale=std)
        
        start = p
        while start > 0 and (p - start) < max_half_width:
            if rowsum[start] < cut_value:
                break
            start -= 1
            
        end = p
        while end < n - 1 and (end - p) < max_half_width:
            if rowsum[end] < cut_value:
                break
            end += 1
            
        width = (end - start + 1)
        
        results.append({
            'peak': p,
            'width': width,
            'stripe_width_range': list(range(start, end + 1))
        })
        
    return results


def get_score_vector(stripe_width_range, rowsum_smoothed_mean, rowsum_smoothed_var, diag_mat, direction_mat, direction, threshold_baseline_p=0.95):

    n=len(stripe_width_range)
    mean_per_pixel = np.full(n, np.nan, dtype=float)
    var_per_pixel = np.full(n, np.nan, dtype=float)
    stripe_score_block = []

    # get signal probabilies in stripe-width range as scores-block 
    scores_width_range = diag_mat[stripe_width_range, :]
    dir_width_range = direction_mat[stripe_width_range, :]
    cond = (dir_width_range == direction) | (dir_width_range == 3)
    stripe_score_block = np.where(cond, scores_width_range, np.nan)

    # get raw score vector by averaging scores-block column-wise 
    score_vec_raw = np.nanmean(stripe_score_block, axis=0) 

    # get mean and variance of score within stripe-width-range
    divisors_width_range = np.sum(~np.isnan(stripe_score_block), axis=1)
    means_width_range = rowsum_smoothed_mean[stripe_width_range]
    vars_width_range = rowsum_smoothed_var[stripe_width_range]
    mean_per_pixel = np.where(
        divisors_width_range > 0, means_width_range / divisors_width_range, 0.0
    )
    var_per_pixel = np.where(
        divisors_width_range > 0, vars_width_range / (divisors_width_range**2), 0.0
    )
    mean_per_width = np.mean(mean_per_pixel)
    var_per_width = np.sum(var_per_pixel)/(n**2)
    
    # get score vector by subtracting baseline threshold from raw scores
    # define the baseline threshold as upper quantile of the score distribution
    baseline_score = norm.ppf(threshold_baseline_p, loc=mean_per_width, scale=np.sqrt(var_per_width))
    score_vector = score_vec_raw - baseline_score
    score_vector = np.nan_to_num(score_vector, nan=0.0)
    score_stats = {'mean' : mean_per_width - baseline_score, 'var' : var_per_width}
    
    return score_vector, score_stats


# Reference for find_max_slice function
#
# Original Source / Repository :
#   - Quagga : https://github.com/dmcbffeng/StripeCaller/blob/master/Quagga/caller/functions.py
#
# Reference :
#   Feng, F. et al. (2026)
#   Statistically rigorous and computationally efficient chromatin stripe detection with Quagga. Genome Res 6, 159–168 727
def find_max_slice(arr):
    """
    Given the enrichment score of each pixel along the candidate stripe,
    find the slice with the largest sum as a called stripe.

    Parameters
    ----------
    arr: np.array (1D)
        Enrichment score

    Returns
    ----------
    head: int
        The start location of the slice

    tail: int
        The end location of the slice

    _max: float
        The sum of scores of the slice

    """
    _max, head, tail = 0, 0, 0
    _max_ending, h, t = 0, 0, 0
    i = 0
    while i < len(arr):
        _max_ending = _max_ending + arr[i]
        if _max_ending < 0:
            h, t = i + 1, i + 1
            _max_ending = 0
        else:
            t = i + 1
        if _max_ending > _max:
            head, tail, _max = h, t, _max_ending
        i += 1
    return head, tail, _max

def load_params(params_path):
    with open(params_path, 'r') as f:
        lines = f.readlines()
        params = {}
        for line in lines:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) == 2:
                key, val = parts[0].strip(), parts[1].strip()
                val = float(val) 
            params[key] = val
    return params

def get_stripe_probability(p_value):

    # calculate stripe probabilities 
    # via EM-based Gaussian mixture clustering on logit-transformed p-values
    calc_signal_p_from_pvalues = robjects.r['calc_signal_p_from_pvalues']
    stripe_p_result = calc_signal_p_from_pvalues(p_value.tolist())
    stripe_probability = np.asarray(stripe_p_result).flatten()

    return stripe_probability


def get_stripe_result(stripe_candidates_info, diag_mat, direction_mat, rowsum_smoothed_mean, rowsum_smoothed_var, direction, resolution=1e4,
                      threshold_baseline_p=0.95, threshold_p_value_stripe=1e-6):

    stripe_result = []
    n_row = diag_mat.shape[0]
    for stripe_info in stripe_candidates_info:
        stripe_width_range = stripe_info['stripe_width_range']
        stripe_pixels = []
        stripe_pass_flag = False

        # get score vector within width range
        scores, score_stats = get_score_vector(stripe_width_range, rowsum_smoothed_mean, rowsum_smoothed_var, diag_mat, direction_mat, 
                         direction=direction, threshold_baseline_p=threshold_baseline_p)

        # call stripe from score vector 
        h, t, max_score = find_max_slice(scores)
        stripe_length_tmp = t-h+1
        stripe_score = max_score / stripe_length_tmp if stripe_length_tmp > 0 else 0.0
        length_start_l = []
        length_end_l = []
        for i in stripe_width_range:
            if direction == 1:
                c_start = min(i+h, n_row-1)
                c_end = min(i+t, n_row-1)
                pixels = [(i, j) for j in range(c_start, c_end + 1)]
            elif direction == 2:
                c_start = max(i-t, 0)
                c_end = max(i-h, 0)
                pixels = [(j, i) for j in range(c_start, c_end + 1)]
            else:
                raise ValueError('Direction must be 1 or 2')
            length_start_l.append(c_start)
            length_end_l.append(c_end)
            stripe_pixels.extend(pixels)
        length_start = min(length_start_l)
        length_end = min(max(length_end_l)+1, n_row-1)
        stripe_length = length_end - length_start
        stripe_length_range = list(range(length_start, length_end+1))
        
        # calculate p-value of stripe
        stripe_p_value = norm.sf(stripe_score, loc=score_stats['mean'], scale=np.sqrt(score_stats['var']))

        # get stripe direction
        if direction==1:
            stripe_direction='5_prime'
        elif direction==2:
            stripe_direction='3_prime'
        else:
            raise ValueError('Direction must be 1 or 2')

        # get stripe result
        stripe_result.append({'width_range' : stripe_width_range, 'length_range' : stripe_length_range, 'pixels' : stripe_pixels, 'direction' : stripe_direction, 
                              'width' : stripe_info['width'], 'length' : stripe_length, 'p_value' : stripe_p_value})

    
    for res in stripe_result:
        if res['length']*resolution < 6e4:
            res['pass'] = (res['p_value'] < threshold_p_value_stripe**2) 
        else:
            res['pass'] = (res['p_value'] < threshold_p_value_stripe) 


    return stripe_result


def call_stripes_on_rowsum(rowsum, chrom, diag_mat, direction_mat, params, direction, output_path, resolution=1e4, max_distance=2e6):

    # Fit Bayesian inference of Normal Dynamic Linear Model and EM on rowsum
    fit_NDLM_and_EM_on_rowsum = robjects.r['fit_NDLM_and_EM_on_rowsum']
    params_r = [ params['smoothing_factor'], params['var_prior_ndlm'], params['threshold_rowsums_em_p'] ]
    mod_result = fit_NDLM_and_EM_on_rowsum(rowsum.tolist(), int(direction), params_r, chrom, resolution, output_path)
    mod = mod_result.rx2('mod')
    
    # NDLM result
    numpy2ri.activate()
    rowsum_smoothed_mean = np.asarray(mod.rx2('mean')).flatten()
    rowsum_smoothed_var = np.asarray(mod.rx2('var')).flatten()

    # high signal row-index and p_values
    high_signal_idx = np.array(mod_result.rx2('high_signal_idx')) - 1
    p_value_rowsum = np.array(mod_result.rx2('p_value'))

    # get genomic index for stripe candidates
    stripe_candidates = get_stripe_candidates(rowsum_smoothed_mean, high_signal_idx, peak_prominence=params['peak_prominence'])

    # get stripe width range
    width_result = get_stripe_width_range(stripe_candidates, rowsum, rowsum_smoothed_mean, rowsum_smoothed_var, p_value_rowsum, 
                                          threshold_mean_to_peak=params['threshold_mean_to_peak'], resolution=resolution, max_width=int(params['max_stripe_width']))

    # call stripes 
    stripe_result = get_stripe_result(width_result, diag_mat, direction_mat, rowsum_smoothed_mean, rowsum_smoothed_var, direction, resolution=resolution,
                                       threshold_baseline_p=params['threshold_baseline_p'], threshold_p_value_stripe=params['threshold_p_value_stripe'])

    return stripe_result


def save_stripe_result_all(stripe_result_list, signal_p_path, output_path, chrom, resolution=1e4):

    filename = f"{chrom}_stripes.bedpe"
    output_file = os.path.join(output_path, filename)

    df_features = pd.read_csv(signal_p_path)
    pixel_dict = {}
    for row in df_features.itertuples(index=False):
        key = (row.idx1, row.idx2)
        pixel_dict[key] = {
            'start1': row.start1,
            'end1': row.end1,
            'start2': row.start2,
            'end2': row.end2,
            'count': row.count,
            'mu': row.mu,
            'signal_p': row.signal_p,
            'p_value_NB_MIX': row.p_value, 
        }

    # get passed stripes
    stripes_all = [res for res in stripe_result_list ]
    passed_stripes = [res for res in stripe_result_list if res.get('pass') == True]
    
    # write stripe bedpe file
    with open(output_file, mode='w', encoding='utf-8') as f:
        
        headers = [
            'chr',
            'start1',
            'end1',
            'start2',
            'end2',
            'stripe_number',
            'stripe_width_range',
            'stripe_length_range',
            'stripe_width',
            'stripe_length',
            'p_value_stripe',
            'direction',
            'count',
            'mu',
            'signal_p',
            'p_value_NB_MIX',
            'pass',
        ]
        f.write('\t'.join(headers) + '\n')

        stripe_number = 1
        for stripe in stripes_all:
            width_range_tmp = np.asarray(stripe.get('width_range')) * resolution
            stripe_width_range = f'{width_range_tmp[0]}-{width_range_tmp[-1]+resolution}'
            stripe_width = stripe.get('width') * resolution 
            length_range_tmp = np.asarray(stripe.get('length_range')) * resolution
            stripe_length_range = f'{length_range_tmp[0]}-{length_range_tmp[-1]}'
            stripe_length = stripe.get('length') * resolution
            p_value_stripe = stripe.get('p_value')
            direction = stripe.get('direction')
            pass_flag = stripe.get('pass')
            pixels = stripe.get('pixels', [])  # [(idx1, idx2), (idx1, idx2), ...]
            
            for idx1, idx2 in pixels:
                pixel_info = pixel_dict.get((idx1, idx2))
                if pixel_info:
                    row_data = [
                        str(chrom),
                        str(pixel_info['start1']),
                        str(pixel_info['end1']),
                        str(pixel_info['start2']),
                        str(pixel_info['end2']),
                        str(stripe_number),
                        str(stripe_width_range),
                        str(stripe_length_range),
                        str(stripe_width),
                        str(stripe_length),
                        str(p_value_stripe),
                        str(direction),
                        str(pixel_info['count']),
                        str(pixel_info['mu']),
                        str(pixel_info['signal_p']),
                        str(pixel_info['p_value_NB_MIX']),
                        str(pass_flag),
                    ]
                    f.write('\t'.join(row_data) + '\n')

            stripe_number += 1

    print(
        f'Successfully saved total {stripe_number} stripes ({len(passed_stripes)}'
        f' passed stripes) into BEDPE format at {output_file}'
    )

        
def call_stripes(signal_p_path, params_path, chrom, resolution=1e4, max_distance=2e6):

    # get signal matrix
    max_d = int(max_distance // resolution)
    diag_to_forward_mat, diag_to_reverse_mat = get_diag_matrix(signal_p_path, max_d)

    # get parameters
    params = load_params(params_path)
    print(params)

    # get direction matrix 
    direction_mat_forward, direction_mat_reverse  = get_stripe_direction_matrix(diag_to_forward_mat, diag_to_reverse_mat, resolution=resolution, window_size=int(params['window_size']))

    # rowsum for calling 5' stripes 
    rowsum_35 = calculate_row_sums(diag_to_forward_mat, direction_mat_forward, direction=1)

    # rowsum for calling 3' stripes 
    rowsum_53 = calculate_row_sums(diag_to_reverse_mat, direction_mat_reverse, direction=2)

    # output path
    plot_path = os.path.dirname(signal_p_path)
    output_path = os.path.dirname(os.path.dirname(signal_p_path))

    # call 3' -> 5' stripes 
    print("Calling 5' Stripes via Bayesian NDLM ...")
    stripe_result_3_5 = call_stripes_on_rowsum(rowsum_35, chrom, diag_to_forward_mat, direction_mat_forward, params, output_path=plot_path, direction=1, resolution=resolution, max_distance=max_distance)

    # call 5' -> 3' stripes 
    print("Calling 3' Stripes via Bayesian NDLM ...")
    stripe_result_5_3 = call_stripes_on_rowsum(rowsum_53, chrom, diag_to_reverse_mat, direction_mat_reverse, params, output_path=plot_path, direction=2, resolution=resolution, max_distance=max_distance)

    # save stripe results
    stripe_result = stripe_result_3_5 + stripe_result_5_3
    save_stripe_result_all(stripe_result, signal_p_path, output_path, chrom, resolution=resolution)
    











    





