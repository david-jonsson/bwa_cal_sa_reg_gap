/*
 * wba_match_cap_coro.h
 *
 *  Created on: Apr 12, 2021
 *      Author: David Jonsson
 */

#ifndef SRC_MATCH_GAP_CORO_H_
#define SRC_MATCH_GAP_CORO_H_



#include "bwt.h"
#include "bwtaln.h"

#ifdef __cplusplus
extern "C" {
#endif

void coro_bwa_cal_sa_reg_gap(int tid, bwt_t *const bwt, int n_seqs, bwa_seq_t *seqs, const gap_opt_t *opt);


#ifdef __cplusplus
}
#endif


#endif /* SRC_MATCH_GAP_CORO_H_ */
