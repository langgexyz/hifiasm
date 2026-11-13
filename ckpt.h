#ifndef __HIFIASM_CKPT__
#define __HIFIASM_CKPT__

#include <stdint.h>
#include <stdio.h>
#include <pthread.h>
#include "CommandLines.h"
#include "Process_Read.h"

/*
 * Write-ahead-log style checkpoints for the error-correction / overlap
 * stages, which have no on-disk state today (the *.bin files are only
 * written once, right before graph construction).
 *
 * Two artifacts per stage:
 *  - <prefix>.r<round>.snap    full snapshot taken when a whole round
 *                             of error correction committed
 *  - <prefix>.r<round>.ecstr   append-only stream of per-read results
 *                             produced by cal_ec_multiple workers
 *  - <prefix>.ff.ovstr         same, for the final overlap pass
 *
 * Only "committed" per-read worker output is logged; every pass after
 * cal_ec_multiple (sl_ec_r, update, post_rev...) is replayed wholesale
 * on resume, so non-idempotent steps are never applied twice.
 */

#define HA_CKPT_HDR_MAGIC   "HICKPT01"
#define HA_CKPT_FIN_MAGIC   0x314E494650544b43ULL
/* v3: flt_tab section is length-prefixed and CRC-covered */
#define HA_CKPT_FMT_VER     3u
#define HA_CKPT_REC_MAGIC   0x5243u

#define HA_CKPT_STAGE_EC    0u
#define HA_CKPT_STAGE_OV    1u
#define HA_CKPT_STAGE_SNAP  2u

typedef struct {
	FILE *fp;
	char *path;
	uint8_t stage;
	int32_t round;
	uint64_t n_reads;
	uint64_t params_hash;
	uint8_t *done;          /* one byte per read */
	uint64_t n_done;
	double last_fsync;
	int32_t fsync_sec;
	pthread_mutex_t mut;
	long resume_off;   /* end offset of valid records after rescan */
} ha_ckpt_stream_t;

/* return a fingerprint of every option/inputs that changes worker output */
uint64_t ha_ckpt_params_hash(void);

/* snapshot: write state committed at end of round `round` (call from main thread) */
int ha_ckpt_write_round(int64_t round, void *flt_tab);
/* find the newest valid snapshot; returns its round or -1. fills *flt_tab. */
int64_t ha_ckpt_load_latest(void **flt_tab);

/* open/create stream for `stage`/`round`; rescans and fills done[] on resume */
ha_ckpt_stream_t *ha_ckpt_stream_open(int stage, int64_t round, uint64_t n_reads);
/* worker-side: append this read's committed results (EC pass) */
void ha_ckpt_dump_ec(ha_ckpt_stream_t *s, uint64_t rid);
/* worker-side: append this read's committed results (final overlap pass) */
void ha_ckpt_dump_ov(ha_ckpt_stream_t *s, uint64_t rid);
void ha_ckpt_stream_close(ha_ckpt_stream_t *s);

/* delete all ckpt artifacts once the definitive *.bin files exist */
void ha_ckpt_purge(void);

/* flt-tab (de)serialization wrappers implemented in htab.cpp */
int ckpt_ft_save(void *flt_tab, FILE *fp);
int ckpt_ft_load(void **flt_tab, FILE *fp);

/* FILE*-level (de)serialization of All_reads, implemented in Process_Read.cpp */
void write_All_reads_fp(All_reads *r, FILE *fp, uint32_t *crc);
int load_All_reads_fp(All_reads *r, FILE *fp, uint32_t *crc);

#endif
