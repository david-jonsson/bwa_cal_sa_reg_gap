/*
 * wba_match_cap_coro.cpp
 *
 *  Created on: Apr 12, 2021
 *      Author: David Jonsson
 */
#include <vector>
#include <cppcoro/generator.hpp>
#include <stdlib.h>
#include <stdio.h>
#include <cstring>

#include <iostream>

#include "bwtgap.h"
#include "match_gap_coro.h"

//-stdlib=libc++  -fcoroutines-ts
//-fcoroutines
#define aln_score(m,o,e,p) ((m)*(p)->s_mm + (o)*(p)->s_gapo + (e)*(p)->s_gape)
#define STATE_M 0
#define STATE_I 1
#define STATE_D 2

#define FORCE_INLINE __attribute__((always_inline)) inline

cppcoro::generator<bwt_aln1_t *> match_gap(int job_nr, bwt_t *const bwt, int len, const ubyte_t *seq, bwt_width_t *width,
                          bwt_width_t *seed_width, const gap_opt_t *opt, int *_n_aln, gap_stack_t *stack);
cppcoro::generator<bool> reg_gap(int *next_job, int tid, bwt_t *const bwt, int n_seqs, bwa_seq_t *seqs, const gap_opt_t *opt);

static inline int __occ_aux(uint64_t y, int c)
{
    // reduce nucleotide counting to bits counting
    y = ((c&2)? y : ~y) >> 1 & ((c&1)? y : ~y) & 0x5555555555555555ull;
    // count the number of 1s in y
    y = (y & 0x3333333333333333ull) + (y >> 2 & 0x3333333333333333ull);
    return ((y + (y >> 4)) & 0xf0f0f0f0f0f0f0full) * 0x101010101010101ull >> 56;
}

inline bwtint_t bwt_occ(const bwt_t *bwt, bwtint_t k, ubyte_t c)
{
    bwtint_t n;
    uint32_t *p, *end;

    if (k == bwt->seq_len) return bwt->L2[c+1] - bwt->L2[c];
    if (k == (bwtint_t)(-1)) return 0;
    k -= (k >= bwt->primary); // because $ is not in bwt

    // retrieve Occ at k/OCC_INTERVAL
    n = ((bwtint_t*)(p = bwt_occ_intv(bwt, k)))[c];
    p += sizeof(bwtint_t); // jump to the start of the first BWT cell

    // calculate Occ up to the last k/32
    end = p + (((k>>5) - ((k&~OCC_INTV_MASK)>>5))<<1);
    for (; p < end; p += 2) n += __occ_aux((uint64_t)p[0]<<32 | p[1], c);

    // calculate Occ
    n += __occ_aux(((uint64_t)p[0]<<32 | p[1]) & ~((1ull<<((~k&31)<<1)) - 1), c);
    if (c == 0) n -= ~k&31; // corrected for the masked bits

    return n;
}

// an analogy to bwt_occ() but more efficient, requiring k <= l
inline void bwt_2occ(const bwt_t *bwt, bwtint_t k, bwtint_t l, ubyte_t c, bwtint_t *ok, bwtint_t *ol)
{
    bwtint_t _k, _l;
    _k = (k >= bwt->primary)? k-1 : k;
    _l = (l >= bwt->primary)? l-1 : l;
    if (_l/OCC_INTERVAL != _k/OCC_INTERVAL || k == (bwtint_t)(-1) || l == (bwtint_t)(-1)) {
        *ok = bwt_occ(bwt, k, c);
        *ol = bwt_occ(bwt, l, c);
    } else {
        bwtint_t m, n, i, j;
        uint32_t *p;
        if (k >= bwt->primary) --k;
        if (l >= bwt->primary) --l;
        n = ((bwtint_t*)(p = bwt_occ_intv(bwt, k)))[c];
        p += sizeof(bwtint_t);
        // calculate *ok
        j = k >> 5 << 5;
        for (i = k/OCC_INTERVAL*OCC_INTERVAL; i < j; i += 32, p += 2)
            n += __occ_aux((uint64_t)p[0]<<32 | p[1], c);
        m = n;
        n += __occ_aux(((uint64_t)p[0]<<32 | p[1]) & ~((1ull<<((~k&31)<<1)) - 1), c);
        if (c == 0) n -= ~k&31; // corrected for the masked bits
        *ok = n;
        // calculate *ol
        j = l >> 5 << 5;
        for (; i < j; i += 32, p += 2)
            m += __occ_aux((uint64_t)p[0]<<32 | p[1], c);
        m += __occ_aux(((uint64_t)p[0]<<32 | p[1]) & ~((1ull<<((~l&31)<<1)) - 1), c);
        if (c == 0) m -= ~l&31; // corrected for the masked bits
        *ol = m;
    }
}

#define __occ_aux4(bwt, b)                                          \
    ((bwt)->cnt_table[(b)&0xff] + (bwt)->cnt_table[(b)>>8&0xff]     \
     + (bwt)->cnt_table[(b)>>16&0xff] + (bwt)->cnt_table[(b)>>24])

inline void bwt_occ4(const bwt_t *bwt, bwtint_t k, bwtint_t cnt[4])
{
    bwtint_t x;
    uint32_t *p, tmp, *end;
    if (k == (bwtint_t)(-1)) {
        memset(cnt, 0, 4 * sizeof(bwtint_t));
        return;
    }
    k -= (k >= bwt->primary); // because $ is not in bwt
    p = bwt_occ_intv(bwt, k);
    // printf("occ_0: %p %p %lu\n", p, bwt, k);

    memcpy(cnt, p, 4 * sizeof(bwtint_t));
    p += sizeof(bwtint_t); // sizeof(bwtint_t) = 4*(sizeof(bwtint_t)/sizeof(uint32_t))
    end = p + ((k>>4) - ((k&~OCC_INTV_MASK)>>4)); // this is the end point of the following loop
    for (x = 0; p < end; ++p) x += __occ_aux4(bwt, *p);
    tmp = *p & ~((1U<<((~k&15)<<1)) - 1);
    x += __occ_aux4(bwt, tmp) - (~k&15);
    cnt[0] += x&0xff; cnt[1] += x>>8&0xff; cnt[2] += x>>16&0xff; cnt[3] += x>>24;
}

// an analogy to bwt_occ4() but more efficient, requiring k <= l
inline void bwt_2occ4(const bwt_t *bwt, bwtint_t k, bwtint_t l, bwtint_t cntk[4], bwtint_t cntl[4])
{
    bwtint_t _k, _l;
    _k = k - (k >= bwt->primary);
    _l = l - (l >= bwt->primary);
    if (_l>>OCC_INTV_SHIFT != _k>>OCC_INTV_SHIFT || k == (bwtint_t)(-1) || l == (bwtint_t)(-1)) {
        bwt_occ4(bwt, k, cntk);
        bwt_occ4(bwt, l, cntl);
    } else {
        bwtint_t x, y;
        uint32_t *p, tmp, *endk, *endl;
        k -= (k >= bwt->primary); // because $ is not in bwt
        l -= (l >= bwt->primary);
        p = bwt_occ_intv(bwt, k);
        // printf("occ_4: %p %p %lu\n", p, bwt, k);

        memcpy(cntk, p, 4 * sizeof(bwtint_t));
        p += sizeof(bwtint_t); // sizeof(bwtint_t) = 4*(sizeof(bwtint_t)/sizeof(uint32_t))
        // prepare cntk[]
        endk = p + ((k>>4) - ((k&~OCC_INTV_MASK)>>4));
        endl = p + ((l>>4) - ((l&~OCC_INTV_MASK)>>4));
        for (x = 0; p < endk; ++p) x += __occ_aux4(bwt, *p);
        y = x;
        tmp = *p & ~((1U<<((~k&15)<<1)) - 1);
        x += __occ_aux4(bwt, tmp) - (~k&15);
        // calculate cntl[] and finalize cntk[]
        for (; p < endl; ++p) y += __occ_aux4(bwt, *p);
        tmp = *p & ~((1U<<((~l&15)<<1)) - 1);
        y += __occ_aux4(bwt, tmp) - (~l&15);
        memcpy(cntl, cntk, 4 * sizeof(bwtint_t));
        cntk[0] += x&0xff; cntk[1] += x>>8&0xff; cntk[2] += x>>16&0xff; cntk[3] += x>>24;
        cntl[0] += y&0xff; cntl[1] += y>>8&0xff; cntl[2] += y>>16&0xff; cntl[3] += y>>24;
    }
}


static void gap_reset_stack(gap_stack_t *stack)
{
    int i;
    for (i = 0; i != stack->n_stacks; ++i)
        stack->stacks[i].n_entries = 0;
    stack->best = stack->n_stacks;
    stack->n_entries = 0;
}

static inline void gap_push(gap_stack_t *stack, int i, bwtint_t k, bwtint_t l, int n_mm, int n_gapo, int n_gape, int n_ins, int n_del,
                            int state, int is_diff, const gap_opt_t *opt)
{
    int score;
    gap_entry_t *p;
    gap_stack1_t *q;
    score = aln_score(n_mm, n_gapo, n_gape, opt);
    q = stack->stacks + score;
    if (q->n_entries == q->m_entries) {
        q->m_entries = q->m_entries? q->m_entries<<1 : 4;
        q->stack = (gap_entry_t*)realloc(q->stack, sizeof(gap_entry_t) * q->m_entries);
    }
    p = q->stack + q->n_entries;
    p->info = (uint32_t)score<<21 | i; p->k = k; p->l = l;
    p->n_mm = n_mm; p->n_gapo = n_gapo; p->n_gape = n_gape;
    p->n_ins = n_ins; p->n_del = n_del;
    p->state = state;
    p->last_diff_pos = is_diff? i : 0;
    ++(q->n_entries);
    ++(stack->n_entries);
    if (stack->best > score) stack->best = score;
}

static inline void gap_pop(gap_stack_t *stack, gap_entry_t *e)
{
    gap_stack1_t *q;
    q = stack->stacks + stack->best;
    *e = q->stack[q->n_entries - 1];
    --(q->n_entries);
    --(stack->n_entries);
    if (q->n_entries == 0 && stack->n_entries) { // reset best
        int i;
        for (i = stack->best + 1; i < stack->n_stacks; ++i)
            if (stack->stacks[i].n_entries != 0) break;
        stack->best = i;
    } else if (stack->n_entries == 0) stack->best = stack->n_stacks;
}

static inline void gap_shadow(int x, int len, bwtint_t max, int last_diff_pos, bwt_width_t *w)
{
    int i, j;
    for (i = j = 0; i < last_diff_pos; ++i) {
        if (w[i].w > x) w[i].w -= x;
        else if (w[i].w == x) {
            w[i].bid = 1;
            w[i].w = max - (++j);
        } // else should not happen
    }
}

static inline int int_log2(uint32_t v)
{
    int c = 0;
    if (v & 0xffff0000u) { v >>= 16; c |= 16; }
    if (v & 0xff00) { v >>= 8; c |= 8; }
    if (v & 0xf0) { v >>= 4; c |= 4; }
    if (v & 0xc) { v >>= 2; c |= 2; }
    if (v & 0x2) c |= 1;
    return c;
}

inline cppcoro::generator<int> match_exact_alt(const bwt_t *bwt, int len, const ubyte_t *str, bwtint_t *k0, bwtint_t *l0)
{
    int i;
    bwtint_t k, l, ok, ol;
    k = *k0; l = *l0;
    bool
        return_zero = false;

    for (i = len - 1; i >= 0; --i) {
        ubyte_t c = str[i];
        if (c > 3)
        {
            return_zero = true;
            break;
        }
//            return 0; // there is an N here. no match
        bwtint_t
            temp_k = k - 1,
            temp_l = l;
        temp_k -= (temp_k >= bwt->primary);
        temp_l -= (temp_l >= bwt->primary);
        uint32_t
            *pre_p = bwt_occ_intv(bwt, temp_k),
            *pre_l = bwt_occ_intv(bwt, temp_l);
        __builtin_prefetch(pre_p);
        __builtin_prefetch(pre_p + 16);
        __builtin_prefetch(pre_l);
        __builtin_prefetch(pre_l + 16);

        co_yield 0;

        bwt_2occ(bwt, k - 1, l, c, &ok, &ol);
        k = bwt->L2[c] + ok + 1;
        l = bwt->L2[c] + ol;
        if (k > l)
        {
            return_zero = true;
            break;
        }
//            return 0; // no match
    }
    bwtint_t
        ret;
    if(return_zero)
        ret = 0;
    else
    {
        *k0 = k; *l0 = l;
        ret = l - k + 1;
    }
    co_yield ret;
}

inline cppcoro::generator<bwt_aln1_t *> match_gap(int job_nr, bwt_t *const bwt, int len, const ubyte_t *seq, bwt_width_t *width,
                          bwt_width_t *seed_width, const gap_opt_t *opt, int *_n_aln, gap_stack_t *stack)
{ // $seq is the reverse complement of the input read
    int best_score = aln_score(opt->max_diff+1, opt->max_gapo+1, opt->max_gape+1, opt);
    int best_diff = opt->max_diff + 1, max_diff = opt->max_diff;
    int best_cnt = 0;
    int max_entries = 0, j, _j, n_aln, m_aln;
    bwt_aln1_t *aln;

    m_aln = 4; n_aln = 0;
    aln = (bwt_aln1_t*)calloc(m_aln, sizeof(bwt_aln1_t));

    // check whether there are too many N
    for (j = _j = 0; j < len; ++j)
        if (seq[j] > 3) ++_j;
//    if (_j > max_diff) {
//        *_n_aln = n_aln;
//        co_yield aln;
//    }
//    else
    if (_j <= max_diff)
    {

        //for (j = 0; j != len; ++j) printf("#0 %d: [%d,%u]\t[%d,%u]\n", j, w[0][j].bid, w[0][j].w, w[1][j].bid, w[1][j].w);
        gap_reset_stack(stack); // reset stack
        gap_push(stack, len, 0, bwt->seq_len, 0, 0, 0, 0, 0, 0, 0, opt);

        while (stack->n_entries) {
            gap_entry_t e;
            int i, m, m_seed = 0, hit_found, allow_diff, allow_M, tmp;
            bwtint_t k, l, cnt_k[4], cnt_l[4], occ;

            if (max_entries < stack->n_entries) max_entries = stack->n_entries;
            if (stack->n_entries > opt->max_entries) break;

//            gap_stack1_t *q;
//            q = stack->stacks + stack->best;
//            __builtin_prefetch(&q->stack[q->n_entries - 1]);
//            co_yield 0;

            gap_pop(stack, &e); // get the best entry
            k = e.k; l = e.l; // SA interval
            i = e.info&0xffff; // length
            if (!(opt->mode & BWA_MODE_NONSTOP) && e.info>>21 > best_score + opt->s_mm) break; // no need to proceed

            m = max_diff - (e.n_mm + e.n_gapo);
            if (opt->mode & BWA_MODE_GAPE) m -= e.n_gape;
            if (m < 0) continue;
            if (seed_width) { // apply seeding
                m_seed = opt->max_seed_diff - (e.n_mm + e.n_gapo);
                if (opt->mode & BWA_MODE_GAPE) m_seed -= e.n_gape;
            }
            //printf("#1\t[%d,%d,%d,%c]\t[%d,%d,%d]\t[%u,%u]\t[%u,%u]\t%d\n", stack->n_entries, a, i, "MID"[e.state], e.n_mm, e.n_gapo, e.n_gape, width[i-1].bid, width[i-1].w, k, l, e.last_diff_pos);
            if (i > 0 && m < width[i-1].bid) continue;

            // check whether a hit is found
            hit_found = 0;
            if (i == 0) hit_found = 1;
            else if (m == 0 && (e.state == STATE_M || (opt->mode&BWA_MODE_GAPE) || e.n_gape == opt->max_gape))
            { // no diff allowed
//                if (bwt_match_exact_alt(bwt, i, seq, &k, &l)) hit_found = 1;

                cppcoro::generator<int>
                    me_alt = match_exact_alt(bwt, i, seq, &k, &l);
                cppcoro::generator<int>::iterator
                    me_alt_it = me_alt.begin();
                while(me_alt_it != me_alt.end())
                {
                    ++me_alt_it;
                    co_yield 0;
                }
                if(*me_alt_it) hit_found = 1;
                else continue; // no hit, skip
            }

            if (hit_found) { // action for found hits
                int score = aln_score(e.n_mm, e.n_gapo, e.n_gape, opt);
                int do_add = 1;
                //printf("#2 hits found: %d:(%u,%u)\n", e.n_mm+e.n_gapo, k, l);
                if (n_aln == 0) {
                    best_score = score;
                    best_diff = e.n_mm + e.n_gapo;
                    if (opt->mode & BWA_MODE_GAPE) best_diff += e.n_gape;
                    if (!(opt->mode & BWA_MODE_NONSTOP))
                        max_diff = (best_diff + 1 > opt->max_diff)? opt->max_diff : best_diff + 1; // top2 behaviour
                }
                if (score == best_score) best_cnt += l - k + 1;
                else if (best_cnt > opt->max_top2) break; // top2b behaviour
                if (e.n_gapo) { // check whether the hit has been found. this may happen when a gap occurs in a tandem repeat
                    for (j = 0; j != n_aln; ++j)
                        if (aln[j].k == k && aln[j].l == l) break;
                    if (j < n_aln) do_add = 0;
                }
                if (do_add) { // append
                    bwt_aln1_t *p;
                    gap_shadow(l - k + 1, len, bwt->seq_len, e.last_diff_pos, width);
                    if (n_aln == m_aln) {
                        m_aln <<= 1;
                        aln = (bwt_aln1_t*)realloc(aln, m_aln * sizeof(bwt_aln1_t));
                        memset(aln + m_aln/2, 0, m_aln/2*sizeof(bwt_aln1_t));
                    }
                    p = aln + n_aln;
                    p->n_mm = e.n_mm; p->n_gapo = e.n_gapo; p->n_gape = e.n_gape;
                    p->n_ins = e.n_ins; p->n_del = e.n_del;
                    p->k = k; p->l = l;
                    p->score = score;
                    //fprintf(stderr, "*** n_mm=%d,n_gapo=%d,n_gape=%d,n_ins=%d,n_del=%d\n", e.n_mm, e.n_gapo, e.n_gape, e.n_ins, e.n_del);
                    ++n_aln;
                }
                continue;
            }

            --i;

            bwtint_t
                temp_k = k - 1,
                temp_l = l;
            temp_k -= (temp_k >= bwt->primary);
            temp_l -= (temp_l >= bwt->primary);
            uint32_t
                *pre_p = bwt_occ_intv(bwt, temp_k),
                *pre_l = bwt_occ_intv(bwt, temp_l);
            __builtin_prefetch(pre_p);
            __builtin_prefetch(pre_p + 16);
            __builtin_prefetch(pre_l);
            __builtin_prefetch(pre_l + 16);

            co_yield 0;

            bwt_2occ4(bwt, k - 1, l, cnt_k, cnt_l); // retrieve Occ values
            occ = l - k + 1;
            // test whether diff is allowed
            allow_diff = allow_M = 1;
            if (i > 0) {
                int ii = i - (len - opt->seed_len);
                if (width[i-1].bid > m-1) allow_diff = 0;
                else if (width[i-1].bid == m-1 && width[i].bid == m-1 && width[i-1].w == width[i].w) allow_M = 0;
                if (seed_width && ii > 0) {
                    if (seed_width[ii-1].bid > m_seed-1) allow_diff = 0;
                    else if (seed_width[ii-1].bid == m_seed-1 && seed_width[ii].bid == m_seed-1
                             && seed_width[ii-1].w == seed_width[ii].w) allow_M = 0;
                }
            }
            // indels
            tmp = (opt->mode & BWA_MODE_LOGGAP)? int_log2(e.n_gape + e.n_gapo)/2+1 : e.n_gapo + e.n_gape;
            if (allow_diff && i >= opt->indel_end_skip + tmp && len - i >= opt->indel_end_skip + tmp) {
                if (e.state == STATE_M) { // gap open
                    if (e.n_gapo < opt->max_gapo) { // gap open is allowed
                        // insertion
                        gap_push(stack, i, k, l, e.n_mm, e.n_gapo + 1, e.n_gape, e.n_ins + 1, e.n_del, STATE_I, 1, opt);
                        // deletion
                        for (j = 0; j != 4; ++j) {
                            k = bwt->L2[j] + cnt_k[j] + 1;
                            l = bwt->L2[j] + cnt_l[j];
                            if (k <= l) gap_push(stack, i + 1, k, l, e.n_mm, e.n_gapo + 1, e.n_gape, e.n_ins, e.n_del + 1, STATE_D, 1, opt);
                        }
                    }
                } else if (e.state == STATE_I) { // extention of an insertion
                    if (e.n_gape < opt->max_gape) // gap extention is allowed
                        gap_push(stack, i, k, l, e.n_mm, e.n_gapo, e.n_gape + 1, e.n_ins + 1, e.n_del, STATE_I, 1, opt);
                } else if (e.state == STATE_D) { // extention of a deletion
                    if (e.n_gape < opt->max_gape) { // gap extention is allowed
                        if (e.n_gape + e.n_gapo < max_diff || occ < opt->max_del_occ) {
                            for (j = 0; j != 4; ++j) {
                                k = bwt->L2[j] + cnt_k[j] + 1;
                                l = bwt->L2[j] + cnt_l[j];
                                if (k <= l) gap_push(stack, i + 1, k, l, e.n_mm, e.n_gapo, e.n_gape + 1, e.n_ins, e.n_del + 1, STATE_D, 1, opt);
                            }
                        }
                    }
                }
            }
            // mismatches
            if (allow_diff && allow_M) { // mismatch is allowed
                for (j = 1; j <= 4; ++j) {
                    int c = (seq[i] + j) & 3;
                    int is_mm = (j != 4 || seq[i] > 3);
                    k = bwt->L2[c] + cnt_k[c] + 1;
                    l = bwt->L2[c] + cnt_l[c];
                    if (k <= l) gap_push(stack, i, k, l, e.n_mm + is_mm, e.n_gapo, e.n_gape, e.n_ins, e.n_del, STATE_M, is_mm, opt);
                }
            } else if (seq[i] < 4) { // try exact match only
                int c = seq[i] & 3;
                k = bwt->L2[c] + cnt_k[c] + 1;
                l = bwt->L2[c] + cnt_l[c];
                if (k <= l) gap_push(stack, i, k, l, e.n_mm, e.n_gapo, e.n_gape, e.n_ins, e.n_del, STATE_M, 0, opt);
            }
        }
    }


    *_n_aln = n_aln;
    //fprintf(stderr, "max_entries = %d\n", max_entries);

    co_yield aln;
}

inline cppcoro::generator<int> cal_width(const bwt_t *bwt, int len, const ubyte_t *str, bwt_width_t *width)
{
    bwtint_t k, l, ok, ol;
    int i, bid;
    bid = 0;
    k = 0; l = bwt->seq_len;
    for (i = 0; i < len; ++i) {
        ubyte_t c = str[i];
        if (c < 4) {
            bwtint_t
                temp_k = k - 1,
                temp_l = l;
            temp_k -= (temp_k >= bwt->primary);
            temp_l -= (temp_l >= bwt->primary);
            uint32_t
                *pre_p = bwt_occ_intv(bwt, temp_k),
                *pre_l = bwt_occ_intv(bwt, temp_l);
            __builtin_prefetch(pre_p);
            __builtin_prefetch(pre_p + 16);
            __builtin_prefetch(pre_l);
            __builtin_prefetch(pre_l + 16);
            co_yield 0;
            bwt_2occ(bwt, k - 1, l, c, &ok, &ol);
            k = bwt->L2[c] + ok + 1;
            l = bwt->L2[c] + ol;
        }
        if (k > l || c > 3) { // then restart
            k = 0;
            l = bwt->seq_len;
            ++bid;
        }
        width[i].w = l - k + 1;
        width[i].bid = bid;
    }
    width[len].w = 0;
    width[len].bid = ++bid;
    co_yield bid;
}

inline cppcoro::generator<bool> reg_gap(int *next_job, int tid, bwt_t *const bwt, int n_seqs, bwa_seq_t *seqs, const gap_opt_t *opt)
{

    int i, j, max_l = 0, max_len;
    gap_stack_t *stack;
    bwt_width_t *w, *seed_w;
    gap_opt_t local_opt = *opt;

    // initiate priority stack
    for (i = max_len = 0; i != n_seqs; ++i)
        if (seqs[i].len > max_len) max_len = seqs[i].len;
    if (opt->fnr > 0.0) local_opt.max_diff = bwa_cal_maxdiff(max_len, BWA_AVG_ERR, opt->fnr);
    if (local_opt.max_diff < local_opt.max_gapo) local_opt.max_gapo = local_opt.max_diff;
    stack = gap_init_stack(local_opt.max_diff, local_opt.max_gapo, local_opt.max_gape, &local_opt);

    seed_w = (bwt_width_t*)calloc(opt->seed_len+1, sizeof(bwt_width_t));
    w = 0;
//    for (i = 0; i != n_seqs; ++i) {
    while(*next_job < n_seqs)
    {
        bwa_seq_t
            *p = seqs + *next_job;
        ++(*next_job);
//        printf("%i %i\n", *next_job, n_seqs);
//        bwa_seq_t *p = seqs + i;
//#ifdef HAVE_PTHREAD
//        if (i % opt->n_threads != tid) continue;
//#endif
        p->sa = 0; p->type = BWA_TYPE_NO_MATCH; p->c1 = p->c2 = 0; p->n_aln = 0; p->aln = 0;
        if (max_l < p->len) {
            max_l = p->len;
            w = (bwt_width_t*)realloc(w, (max_l + 1) * sizeof(bwt_width_t));
            memset(w, 0, (max_l + 1) * sizeof(bwt_width_t));
        }
        //bwt_cal_width(bwt, p->len, p->seq, w);
        {
            cppcoro::generator<int>
                cal_w = cal_width(bwt, p->len, p->seq, w);
            cppcoro::generator<int>::iterator
                cal_w_it = cal_w.begin();
            while(cal_w_it != cal_w.end())
            {
                ++cal_w_it;
                co_yield 0;
            }
        }
        if (opt->fnr > 0.0) local_opt.max_diff = bwa_cal_maxdiff(p->len, BWA_AVG_ERR, opt->fnr);
        local_opt.seed_len = opt->seed_len < p->len? opt->seed_len : 0x7fffffff;
        if (p->len > opt->seed_len)
        {
            cppcoro::generator<int>
                cal_w = cal_width(bwt, opt->seed_len, p->seq + (p->len - opt->seed_len), seed_w);
            cppcoro::generator<int>::iterator
                cal_w_it = cal_w.begin();
            while(cal_w_it != cal_w.end())
            {
                ++cal_w_it;
                co_yield 0;
            }
        }
            //bwt_cal_width(bwt, opt->seed_len, p->seq + (p->len - opt->seed_len), seed_w);
        // core function
        for (j = 0; j < p->len; ++j) // we need to complement
            p->seq[j] = p->seq[j] > 3? 4 : 3 - p->seq[j];
        {
            cppcoro::generator<bwt_aln1_t *>
                m_gap = match_gap(*next_job - 1,bwt, p->len, p->seq, w, p->len <= opt->seed_len? 0 : seed_w, &local_opt, &p->n_aln, stack);
            cppcoro::generator<bwt_aln1_t *>::iterator
                m_gap_it = m_gap.begin();

            while(m_gap_it != m_gap.end())
            {
                ++m_gap_it;
                co_yield 0;
            }

            p->aln = *m_gap_it;
        }

        //fprintf(stderr, "mm=%lld,ins=%lld,del=%lld,gapo=%lld\n", p->aln->n_mm, p->aln->n_ins, p->aln->n_del, p->aln->n_gapo);
        // clean up the unused data in the record
        free(p->name); free(p->seq); free(p->rseq); free(p->qual);
        p->name = 0; p->seq = p->rseq = p->qual = 0;
    }
    free(seed_w); free(w);
    gap_destroy_stack(stack);
    co_yield 1;
}


#define NR_COROS 3
void coro_bwa_cal_sa_reg_gap(int tid, bwt_t *const bwt, int n_seqs, bwa_seq_t *seqs, const gap_opt_t *opt)
{

    bool
        done = false;
    int
        next_job = 0;

    cppcoro::generator<bool>
        reg_gaps[NR_COROS] =
        {
             reg_gap(&next_job, tid, bwt, n_seqs, seqs, opt),
             reg_gap(&next_job, tid, bwt, n_seqs, seqs, opt),
             reg_gap(&next_job, tid, bwt, n_seqs, seqs, opt),
//             reg_gap(&next_job, tid, bwt, n_seqs, seqs, opt),
        };
    cppcoro::generator<bool>::iterator
        its[NR_COROS];

    for(int i = 0; i < NR_COROS; ++i)
        its[i] = reg_gaps[i].begin();

    while(!done)
    {
        done = true;
        for(int i = 0; i < NR_COROS; ++i)
        {
            ++its[i];
            done &= *its[i];
        }
    }
}
