#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <sys/stat.h>
#include <zlib.h>
#include "ckpt.h"
#include "ecovlp.h"
#include "htab.h"
#include "Overlaps.h"

extern All_reads R_INF;
extern cc_v scc, scb, sca;
extern uint32_t *het_cnt;
extern hifiasm_opt_t asm_opt;

/*******************************************************************************************
 * on-disk formats
 *
 * snapshot  <prefix>.r<round>.snap
 *   hdr : magic(8) | fmt_ver | stage=SNAP | round | n_reads | params_hash
 *   body: has_flt | [yak_ft] | All_reads | paf[] | revpaf[] |
 *         has_scc | [cc_v]   | has_scb   | [cc_v] | has_sca | [cc_v]
 *   fin : body_crc32 | FIN_MAGIC          (written via <file>.tmp + rename)
 *
 * stream    <prefix>.r<round>.ecstr / <prefix>.ff.ovstr
 *   hdr : same, stage=EC/OV
 *   rec : REC_MAGIC | stage | rid | plen | payload | crc32(payload)
 *   EC payload  : trio_flag | het | scc{n,f,a} | scb{n,a} | sca{n,a} |
 *                 paf{fc,ab,n,buffer} | rpaf{fc,ab,n,buffer}
 *   OV payload  : paf{n,buffer} | rpaf{n,buffer}
 *******************************************************************************************/

/* crc32() takes uInt (32-bit); chunk so >4GB sections are fully covered */
static uint32_t ckpt_crc32(uint32_t crc, const void *p, uint64_t n)
{
	const uint8_t *b = (const uint8_t*)p;
	while (n) {
		uInt c = (n > 0x7fffffffULL) ? 0x7fffffffU : (uInt)n;
		crc = crc32(crc, b, c);
		b += c; n -= c;
	}
	return crc;
}

#define CKPT_W(fp, p, n, c) do { \
	fwrite((p), 1, (n), (fp)); \
	if (c) *(c) = ckpt_crc32(*(c), (const Bytef*)(p), (n)); \
} while (0)

#define CKPT_RD(fp, p, n, c) do { \
	if (fread((p), 1, (n), (fp)) != (size_t)(n)) return 0; \
	if (c) *(c) = ckpt_crc32(*(c), (const Bytef*)(p), (n)); \
} while (0)

static uint64_t fnv1a(uint64_t h, const void *p, uint64_t n)
{
	const uint8_t *b = (const uint8_t*)p;
	uint64_t i;
	for (i = 0; i < n; ++i) { h ^= b[i]; h *= 1099511628211ULL; }
	return h;
}

uint64_t ha_ckpt_params_hash(void)
{
	uint64_t h = 14695981039346656037ULL;
	const hifiasm_opt_t *o = &asm_opt;
	int i;
	struct stat st;
	char buf[4096];
	FILE *fp;

	h = fnv1a(h, HA_VERSION, strlen(HA_VERSION));
	h = fnv1a(h, &o->k_mer_length, sizeof(o->k_mer_length));
	h = fnv1a(h, &o->mz_win, sizeof(o->mz_win));
	h = fnv1a(h, &o->bf_shift, sizeof(o->bf_shift));
	h = fnv1a(h, &o->high_factor, sizeof(o->high_factor));
	h = fnv1a(h, &o->max_ov_diff_ec, sizeof(o->max_ov_diff_ec));
	h = fnv1a(h, &o->max_ov_diff_final, sizeof(o->max_ov_diff_final));
	h = fnv1a(h, &o->number_of_round, sizeof(o->number_of_round));
	h = fnv1a(h, &o->number_of_pround, sizeof(o->number_of_pround));
	h = fnv1a(h, &o->adapterLen, sizeof(o->adapterLen));
	h = fnv1a(h, &o->min_overlap_Len, sizeof(o->min_overlap_Len));
	h = fnv1a(h, &o->max_hang_Len, sizeof(o->max_hang_Len));
	h = fnv1a(h, &o->gap_fuzz, sizeof(o->gap_fuzz));
	h = fnv1a(h, &o->min_drop_rate, sizeof(o->min_drop_rate));
	h = fnv1a(h, &o->max_drop_rate, sizeof(o->max_drop_rate));
	h = fnv1a(h, &o->flag, sizeof(o->flag));
	h = fnv1a(h, &o->polyploidy, sizeof(o->polyploidy));
	h = fnv1a(h, &o->dp_min_len, sizeof(o->dp_min_len));
	h = fnv1a(h, &o->dp_e, sizeof(o->dp_e));
	h = fnv1a(h, &o->is_ont, sizeof(o->is_ont));
	h = fnv1a(h, &o->is_sc, sizeof(o->is_sc));
	h = fnv1a(h, &o->rl_cut, sizeof(o->rl_cut));
	h = fnv1a(h, &o->sc_cut, sizeof(o->sc_cut));
	h = fnv1a(h, &o->mz_rewin, sizeof(o->mz_rewin));
	h = fnv1a(h, &o->mz_sample_dist, sizeof(o->mz_sample_dist));
	h = fnv1a(h, &o->max_kmer_cnt, sizeof(o->max_kmer_cnt));
	h = fnv1a(h, &o->min_hist_kmer_cnt, sizeof(o->min_hist_kmer_cnt));
	h = fnv1a(h, &o->is_dbg_het_cnt, sizeof(o->is_dbg_het_cnt));
	h = fnv1a(h, &o->hg_size, sizeof(o->hg_size));
	h = fnv1a(h, &o->num_reads, sizeof(o->num_reads));

	for (i = 0; i < o->num_reads; ++i) {
		if (!o->read_file_names[i]) continue;
		h = fnv1a(h, o->read_file_names[i], strlen(o->read_file_names[i]));
		if (stat(o->read_file_names[i], &st) == 0) {
			h = fnv1a(h, &st.st_size, sizeof(st.st_size));
			h = fnv1a(h, &st.st_mtime, sizeof(st.st_mtime));
		}
		fp = fopen(o->read_file_names[i], "r");
		if (fp) {
			size_t got = fread(buf, 1, sizeof(buf), fp);
			h = fnv1a(h, buf, got);
			fclose(fp);
		}
	}
	return h;
}

static int ha_ckpt_on(void)
{
	return asm_opt.load_index_from_disk && asm_opt.ckpt_int != 0 &&
		!(asm_opt.flag & HA_F_VERBOSE_GFA) && !asm_opt.sec_in && !asm_opt.dbg_ovec_cal;
}

/****************************** small helpers ******************************/

static char *ckpt_name(const char *prefix, const char *fmt, int64_t round)
{
	char *n = (char*)malloc(strlen(prefix) + 64);
	if (round >= 0) sprintf(n, fmt, prefix, (long long)round);
	else sprintf(n, fmt, prefix);
	return n;
}

static FILE *ckpt_fopen_atomic(const char *path, char **tmp_path)
{
	*tmp_path = (char*)malloc(strlen(path) + 8);
	sprintf(*tmp_path, "%s.tmp", path);
	return fopen(*tmp_path, "w");
}

static int ckpt_rename_sync(FILE *fp, const char *tmp, const char *path)
{
	int bad = ferror(fp);   /* catches any earlier failed fwrite */
	if (fflush(fp) != 0) bad = 1;
	if (fileno(fp) >= 0 && fsync(fileno(fp)) != 0) bad = 1;
	if (fclose(fp) != 0) bad = 1;
	if (bad || rename(tmp, path) != 0) { unlink(tmp); free((void*)tmp); return 0; }
	free((void*)tmp);
	return 1;
}

static int wr_hdr(FILE *fp, uint8_t stage, int64_t round, uint64_t n_reads)
{
	uint32_t v32;
	fwrite(HA_CKPT_HDR_MAGIC, 1, 8, fp);
	v32 = HA_CKPT_FMT_VER;   fwrite(&v32, sizeof(v32), 1, fp);
	v32 = stage;             fwrite(&v32, sizeof(v32), 1, fp);
	fwrite(&round, sizeof(round), 1, fp);
	fwrite(&n_reads, sizeof(n_reads), 1, fp);
	uint64_t ph = ha_ckpt_params_hash();
	fwrite(&ph, sizeof(ph), 1, fp);
	return 1;
}

/* header check shared by snap + stream scan; consumes fp position past hdr */
static int rd_hdr(FILE *fp, uint8_t want_stage, int64_t *round, uint64_t *n_reads)
{
	char magic[8];
	uint32_t v32, stage;
	uint64_t ph;
	if (fread(magic, 1, 8, fp) != 8 || memcmp(magic, HA_CKPT_HDR_MAGIC, 8)) return 0;
	if (fread(&v32, sizeof(v32), 1, fp) != 1 || v32 != HA_CKPT_FMT_VER) return 0;
	if (fread(&stage, sizeof(stage), 1, fp) != 1 || stage != want_stage) return 0;
	if (fread(round, sizeof(*round), 1, fp) != 1) return 0;
	if (fread(n_reads, sizeof(*n_reads), 1, fp) != 1) return 0;
	if (fread(&ph, sizeof(ph), 1, fp) != 1) return 0;
	if (ph != ha_ckpt_params_hash()) return -1;   /* stale: params differ */
	return 1;
}

/****************************** snapshot (L1) ******************************/

static void wr_paf_arr(FILE *fp, ma_hit_t_alloc *paf, uint64_t n, uint32_t *crc)
{
	uint64_t i;
	for (i = 0; i < n; ++i) {
		CKPT_W(fp, &paf[i].is_fully_corrected, sizeof(uint8_t), crc);
		CKPT_W(fp, &paf[i].is_abnormal, sizeof(uint8_t), crc);
		CKPT_W(fp, &paf[i].length, sizeof(paf[i].length), crc);
		if (paf[i].length) CKPT_W(fp, paf[i].buffer, sizeof(ma_hit_t)*paf[i].length, crc);
	}
}

static int rd_paf_arr(FILE *fp, ma_hit_t_alloc *paf, uint64_t n, uint32_t *crc, uint64_t fsz)
{
	uint64_t i;
	for (i = 0; i < n; ++i) {
		CKPT_RD(fp, &paf[i].is_fully_corrected, sizeof(uint8_t), crc);
		CKPT_RD(fp, &paf[i].is_abnormal, sizeof(uint8_t), crc);
		CKPT_RD(fp, &paf[i].length, sizeof(paf[i].length), crc);
		paf[i].size = paf[i].length;
		paf[i].buffer = NULL;
		if (paf[i].length) {
			/* a stored buffer can't exceed the file it lives in */
			if (paf[i].length > fsz/sizeof(ma_hit_t)) return 0;
			paf[i].buffer = (ma_hit_t*)malloc(sizeof(ma_hit_t)*paf[i].length);
			if (!paf[i].buffer) return 0;
			CKPT_RD(fp, paf[i].buffer, sizeof(ma_hit_t)*paf[i].length, crc);
		}
	}
	return 1;
}

static void wr_cc(FILE *fp, const cc_v *cc, uint64_t n, uint32_t *crc)
{
	uint64_t i;
	for (i = 0; i < n; ++i) {
		uint64_t m = cc->a[i].n;
		uint8_t f = cc->f ? cc->f[i] : 0;
		CKPT_W(fp, &m, sizeof(m), crc);
		CKPT_W(fp, &f, sizeof(f), crc);
		if (m) CKPT_W(fp, cc->a[i].a, sizeof(uint16_t)*m, crc);
	}
}

static int rd_cc(FILE *fp, cc_v *cc, uint64_t n, uint32_t *crc, uint64_t fsz)
{
	uint64_t i;
	for (i = 0; i < n; ++i) {
		uint64_t m = 0; uint8_t f = 0;
		CKPT_RD(fp, &m, sizeof(m), crc);
		CKPT_RD(fp, &f, sizeof(f), crc);
		if (m > fsz/sizeof(uint16_t)) return 0;
		if (m > cc->a[i].m) { REALLOC(cc->a[i].a, m); if (!cc->a[i].a) return 0; cc->a[i].m = m; }
		cc->a[i].n = m;
		if (m) CKPT_RD(fp, cc->a[i].a, sizeof(uint16_t)*m, crc);
		if (cc->f) cc->f[i] = f;
	}
	return 1;
}

static long sec_mrk(FILE *fp, long *prev)
{
	long now = ftell(fp), sz = now - *prev;
	*prev = now;
	return sz;
}

static int ckpt_write_one(const char *path, int64_t round, void *flt_tab,
		long *sz_flt, long *sz_reads, long *sz_paf, long *sz_tr)
{
	char *tmp = NULL;
	FILE *fp = ckpt_fopen_atomic(path, &tmp);
	uint32_t crc = 0;
	uint64_t u64;
	uint8_t u8;
	long prev;
	if (!fp) { free(tmp); return 0; }

	wr_hdr(fp, HA_CKPT_STAGE_SNAP, round, R_INF.total_reads);
	prev = ftell(fp);

	u8 = (flt_tab != NULL); CKPT_W(fp, &u8, 1, &crc);
	if (u8) {
		/* buffer the opaque flt_tab bytes; carry a section crc so the
		   reader verifies it BEFORE parsing (whole-file crc is only
		   checked at the footer, too late) */
		char *mb = NULL; size_t mn = 0;
		FILE *mf = open_memstream(&mb, &mn);
		int ok = mf && ckpt_ft_save(flt_tab, mf);
		if (mf) fclose(mf);
		if (!ok) { free(mb); fclose(fp); unlink(tmp); free(tmp); return 0; }
		u64 = mn; CKPT_W(fp, &u64, sizeof(u64), &crc);
		{ uint32_t scrc = ckpt_crc32(0, mb, mn);
		  CKPT_W(fp, &scrc, sizeof(scrc), &crc); }
		CKPT_W(fp, mb, mn, &crc);
		free(mb);
	}
	if (sz_flt) *sz_flt = sec_mrk(fp, &prev);

	write_All_reads_fp(&R_INF, fp, &crc);
	if (sz_reads) *sz_reads = sec_mrk(fp, &prev);

	wr_paf_arr(fp, R_INF.paf, R_INF.total_reads, &crc);
	wr_paf_arr(fp, R_INF.reverse_paf, R_INF.total_reads, &crc);
	if (sz_paf) *sz_paf = sec_mrk(fp, &prev);

	u8 = (scc.a != NULL); CKPT_W(fp, &u8, 1, &crc);
	if (u8) wr_cc(fp, &scc, R_INF.total_reads, &crc);
	u8 = (scb.a != NULL); CKPT_W(fp, &u8, 1, &crc);
	if (u8) wr_cc(fp, &scb, R_INF.total_reads, &crc);
	u8 = (sca.a != NULL); CKPT_W(fp, &u8, 1, &crc);
	if (u8) wr_cc(fp, &sca, R_INF.total_reads, &crc);
	if (sz_tr) *sz_tr = sec_mrk(fp, &prev);

	/* opt fields that workers read and that get mutated mid-run by
	   ha_opt_update_cov (not part of params_hash: the snapped value is the
	   already-bumped one, restoring it keeps resume identical) */
	{ int64_t v = asm_opt.max_n_chain;    CKPT_W(fp, &v, sizeof(v), &crc); }

	u64 = crc;   fwrite(&u64, sizeof(u64), 1, fp);
	u64 = HA_CKPT_FIN_MAGIC; fwrite(&u64, sizeof(u64), 1, fp);

	return ckpt_rename_sync(fp, tmp, path);
}

int ha_ckpt_write_round(int64_t round, void *flt_tab)
{
	char *path, *old;
	int ok;
	long sf = -1, sr = -1, sp = -1, st = -1;
	if (!ha_ckpt_on()) return 0;
	path = ckpt_name(asm_opt.output_file_name, "%s.r%lld.snap", round);
	ok = ckpt_write_one(path, round, flt_tab, &sf, &sr, &sp, &st);
	fprintf(stderr, "[M::%s::%.3f] %s checkpoint for round %lld"
		" (flt=%.1fMB reads=%.1fMB paf=%.1fMB tr=%.1fMB)\n", __func__,
		yak_realtime(), ok?"==> wrote":"!! failed to write", (long long)round,
		sf/1048576.0, sr/1048576.0, sp/1048576.0, st/1048576.0);
	if (ok) {
		old = ckpt_name(asm_opt.output_file_name, "%s.r%lld.ecstr", round);
		unlink(old); free(old); /* committed round: stream no longer needed */
	}
	if (ok && round > 0) {      /* keep only the newest snapshot */
		old = ckpt_name(asm_opt.output_file_name, "%s.r%lld.snap", round-1);
		unlink(old); free(old);
	}
	free(path);
	return ok;
}

/* free whatever load_All_reads_fp managed to allocate before failing;
   pointer arrays are calloc'd so unfilled slots are NULL and safe */
static void free_partial_reads(All_reads *r)
{
	uint64_t i;
	/* total_reads may be a garbage value from a truncated header; only
	   iterate when at least one per-read array was actually allocated */
	if (r->N_site || r->read_sperate || r->cigars || r->second_round_cigar || r->rsc)
		for (i = 0; i < r->total_reads; ++i) {
			if (r->N_site) free(r->N_site[i]);
			if (r->read_sperate) free(r->read_sperate[i]);
			if (r->cigars) { free(r->cigars[i].record); free(r->cigars[i].lost_base); }
			if (r->second_round_cigar) { free(r->second_round_cigar[i].record); free(r->second_round_cigar[i].lost_base); }
			if (r->rsc) free(r->rsc[i]);
		}
	free(r->N_site); free(r->read_length); free(r->read_size);
	free(r->read_sperate); free(r->name); free(r->name_index);
	free(r->trio_flag); free(r->cigars); free(r->second_round_cigar); free(r->rsc);
}

static void free_paf_arr(ma_hit_t_alloc *paf, uint64_t n)
{
	uint64_t i;
	if (!paf) return;
	for (i = 0; i < n; ++i) free(paf[i].buffer);
	free(paf);
}

static void free_cc_v(cc_v *c)
{
	uint64_t i;
	if (!c->a) return;
	for (i = 0; i < c->m; ++i) free(c->a[i].a);
	free(c->a); free(c->f);
	c->a = NULL; c->f = NULL; c->n = c->m = 0;
}

static int ckpt_read_one(const char *path, int64_t exp_round, void **flt_tab, All_reads *tmp)
{
	FILE *fp = fopen(path, "r");
	int64_t round = -1; uint64_t n_reads = 0, u64;
	uint32_t crc = 0, file_crc;
	uint8_t u8;
	int rh;
	struct stat st;
	uint64_t fsz = 0;
	void *lft = NULL;
	cc_v lscc, lscb, lsca;
	int64_t lchain = 0;
	if (!fp) return 0;
	memset(&lscc, 0, sizeof(lscc));
	memset(&lscb, 0, sizeof(lscb));
	memset(&lsca, 0, sizeof(lsca));
	if (fstat(fileno(fp), &st) == 0) fsz = (uint64_t)st.st_size;
	rh = rd_hdr(fp, HA_CKPT_STAGE_SNAP, &round, &n_reads);
	if (rh <= 0 || round != exp_round || n_reads == 0 || n_reads > fsz/16)
		{ fprintf(stderr, "[ckpt-dbg] %s hdr fail rh=%d round=%lld n=%lu\n", path, rh, (long long)round, (unsigned long)n_reads); goto fail; }

#define CKPT_RD_GF(p, n) do { \
	if (fread((p), 1, (n), fp) != (size_t)(n)) goto fail; \
	crc = ckpt_crc32(crc, (const Bytef*)(p), (n)); \
} while (0)

	/* pass 1: verify the whole body before ANY parsing. Read the footer
	   first, then crc everything between hdr and footer. Corrupt or
	   truncated files are rejected before a single field is trusted. */
	{	long fpos = ftell(fp);
		if (fsz < 16 || fseek(fp, (long)fsz - 16, SEEK_SET) != 0) goto fail;
		if (fread(&u64, sizeof(u64), 1, fp) != 1) goto fail;
		file_crc = (uint32_t)u64;
		if (fread(&u64, sizeof(u64), 1, fp) != 1 || u64 != HA_CKPT_FIN_MAGIC)
			{ fprintf(stderr, "[ckpt-dbg] %s fin fail\n", path); goto fail; }
		if (fseek(fp, fpos, SEEK_SET) != 0) goto fail;
		{	uint8_t buf[1<<20];
			uint64_t left = fsz - 16 - (uint64_t)fpos;
			while (left) {
				size_t c = left > sizeof(buf) ? sizeof(buf) : (size_t)left;
				if (fread(buf, 1, c, fp) != c) goto fail;
				crc = ckpt_crc32(crc, buf, c);
				left -= c;
			}
		}
		if (file_crc != crc) { fprintf(stderr, "[ckpt-dbg] %s crc mismatch %x vs %x\n", path, file_crc, crc); goto fail; }
		if (fseek(fp, fpos, SEEK_SET) != 0) goto fail;
	}
	crc = 0;   /* pass 2 re-accumulates to keep readers side-effect-free */

	CKPT_RD_GF(&u8, 1);
	if (u8) {
		uint64_t mn; uint8_t *mb = NULL; FILE *mf; int ok; uint32_t scrc;
		CKPT_RD_GF(&mn, sizeof(mn));
		CKPT_RD_GF(&scrc, sizeof(scrc));
		if (mn > fsz) goto fail;
		if (mn) {
			mb = (uint8_t*)malloc(mn);
			if (!mb) goto fail;
			if (fread(mb, 1, mn, fp) != mn) { free(mb); goto fail; }
			crc = ckpt_crc32(crc, mb, mn);
			/* verify the section before parsing it */
			if (ckpt_crc32(0, mb, mn) != scrc) { free(mb); fprintf(stderr, "[ckpt-dbg] %s ft crc\n", path); goto fail; }
			mf = fmemopen(mb, mn, "r");
			ok = mf && ckpt_ft_load(&lft, mf);
			if (mf) fclose(mf);
			free(mb);
			if (!ok) { fprintf(stderr, "[ckpt-dbg] %s ft fail\n", path); goto fail; }
		}
	}

	if (!load_All_reads_fp(tmp, fp, &crc)) { fprintf(stderr, "[ckpt-dbg] %s reads fail\n", path); goto fail; }
	if (tmp->total_reads != n_reads) { fprintf(stderr, "[ckpt-dbg] %s n_reads mismatch\n", path); goto fail; }

	tmp->paf = (ma_hit_t_alloc*)calloc(n_reads, sizeof(ma_hit_t_alloc));
	tmp->reverse_paf = (ma_hit_t_alloc*)calloc(n_reads, sizeof(ma_hit_t_alloc));
	if (!tmp->paf || !tmp->reverse_paf) goto fail;
	if (!rd_paf_arr(fp, tmp->paf, n_reads, &crc, fsz)) { fprintf(stderr, "[ckpt-dbg] %s paf fail\n", path); goto fail; }
	if (!rd_paf_arr(fp, tmp->reverse_paf, n_reads, &crc, fsz)) { fprintf(stderr, "[ckpt-dbg] %s rpaf fail\n", path); goto fail; }

	CKPT_RD_GF(&u8, 1);
	if (u8) {
		lscc.n = lscc.m = n_reads; CALLOC(lscc.a, n_reads); CALLOC(lscc.f, n_reads);
		if (!lscc.a || !lscc.f) goto fail;
		if (!rd_cc(fp, &lscc, n_reads, &crc, fsz)) { fprintf(stderr, "[ckpt-dbg] %s scc fail\n", path); goto fail; }
	}
	CKPT_RD_GF(&u8, 1);
	if (u8) {
		lscb.n = lscb.m = n_reads; CALLOC(lscb.a, n_reads);
		if (!lscb.a) goto fail;
		if (!rd_cc(fp, &lscb, n_reads, &crc, fsz)) { fprintf(stderr, "[ckpt-dbg] %s scb fail\n", path); goto fail; }
	}
	CKPT_RD_GF(&u8, 1);
	if (u8) {
		lsca.n = lsca.m = n_reads; CALLOC(lsca.a, n_reads);
		if (!lsca.a) goto fail;
		if (!rd_cc(fp, &lsca, n_reads, &crc, fsz)) { fprintf(stderr, "[ckpt-dbg] %s sca fail\n", path); goto fail; }
	}

	CKPT_RD_GF(&lchain, sizeof(lchain));

	if (fread(&u64, sizeof(u64), 1, fp) != 1) { fprintf(stderr, "[ckpt-dbg] %s crc missing\n", path); goto fail; }
	file_crc = (uint32_t)u64;
	if (fread(&u64, sizeof(u64), 1, fp) != 1 || u64 != HA_CKPT_FIN_MAGIC) { fprintf(stderr, "[ckpt-dbg] %s fin fail\n", path); goto fail; }
	if (file_crc != crc) { fprintf(stderr, "[ckpt-dbg] %s crc mismatch %x vs %x\n", path, file_crc, crc); goto fail; }
	fclose(fp);

	/* commit: file fully verified, now update globals */
	if (lft) { if (*flt_tab) ha_ft_destroy(*flt_tab); *flt_tab = lft; }
	free_cc_v(&scc); scc = lscc;
	free_cc_v(&scb); scb = lscb;
	free_cc_v(&sca); sca = lsca;
	asm_opt.max_n_chain = lchain;
	return 1;

fail:
	fclose(fp);
	if (lft) ha_ft_destroy(lft);
	free_cc_v(&lscc); free_cc_v(&lscb); free_cc_v(&lsca);
	free_paf_arr(tmp->paf, n_reads);
	free_paf_arr(tmp->reverse_paf, n_reads);
	free_partial_reads(tmp);
	memset(tmp, 0, sizeof(*tmp));
	return 0;
}

int64_t ha_ckpt_load_latest(void **flt_tab)
{
	int64_t r;
	char *path;
	All_reads tmp;
	uint64_t i;
	if (!ha_ckpt_on()) return -1;
	for (r = asm_opt.number_of_round - 1; r >= 0; --r) {
		path = ckpt_name(asm_opt.output_file_name, "%s.r%lld.snap", r);
		memset(&tmp, 0, sizeof(tmp));
		if (ckpt_read_one(path, r, flt_tab, &tmp)) {
			/* drop any partially-loaded R_INF (e.g. a lone *.ec.bin) */
			if (R_INF.total_reads && R_INF.N_site && R_INF.read_sperate)
				destory_read_bin(&R_INF);
			if (R_INF.paf) {
				for (i = 0; i < R_INF.total_reads; ++i) free(R_INF.paf[i].buffer);
				free(R_INF.paf);
			}
			if (R_INF.reverse_paf) {
				for (i = 0; i < R_INF.total_reads; ++i) free(R_INF.reverse_paf[i].buffer);
				free(R_INF.reverse_paf);
			}
			if (R_INF.rsc) {
				for (i = 0; i < R_INF.total_reads; ++i) free(R_INF.rsc[i]);
				free(R_INF.rsc);
			}
			memset(&R_INF, 0, sizeof(R_INF));
			R_INF = tmp;   /* struct copy: takes over buffers */
			fprintf(stderr, "[M::%s::%.3f] ==> resumed checkpoint of round %lld\n",
				__func__, yak_realtime(), (long long)r);
			free(path);
			return r;
		}
		free(path);
	}
	return -1;
}

/****************************** stream (L2) ******************************/

static char *stream_name(int stage, int64_t round)
{
	if (stage == HA_CKPT_STAGE_EC)
		return ckpt_name(asm_opt.output_file_name, "%s.r%lld.ecstr", round);
	return ckpt_name(asm_opt.output_file_name, "%s.ff.ovstr", -1);
}

/* serialize one read's committed ec results into stream buffer */
static uint64_t ec_payload(const ha_ckpt_stream_t *s, uint64_t rid, uint8_t *buf)
{
	uint8_t *w = buf;
	uint64_t n;
	ma_hit_t_alloc *p;
	(void)s;

	memcpy(w, &R_INF.trio_flag[rid], 1); w += 1;
	{ uint32_t h = het_cnt ? het_cnt[rid] : 0; memcpy(w, &h, 4); w += 4; }

	n = scc.a ? scc.a[rid].n : 0;
	memcpy(w, &n, 8); w += 8;
	{ uint8_t f = (scc.a && scc.f) ? scc.f[rid] : 0; memcpy(w, &f, 1); w += 1; }
	if (n) { memcpy(w, scc.a[rid].a, 2*n); w += 2*n; }

	n = scb.a ? scb.a[rid].n : 0;
	memcpy(w, &n, 8); w += 8;
	if (n) { memcpy(w, scb.a[rid].a, 2*n); w += 2*n; }

	n = sca.a ? sca.a[rid].n : 0;
	memcpy(w, &n, 8); w += 8;
	if (n) { memcpy(w, sca.a[rid].a, 2*n); w += 2*n; }

	p = &R_INF.paf[rid];
	memcpy(w, &p->is_fully_corrected, 1); w += 1;
	memcpy(w, &p->is_abnormal, 1); w += 1;
	n = p->length; memcpy(w, &n, 8); w += 8;
	if (p->length) { memcpy(w, p->buffer, sizeof(ma_hit_t)*p->length); w += sizeof(ma_hit_t)*p->length; }

	p = &R_INF.reverse_paf[rid];
	memcpy(w, &p->is_fully_corrected, 1); w += 1;
	memcpy(w, &p->is_abnormal, 1); w += 1;
	n = p->length; memcpy(w, &n, 8); w += 8;
	if (p->length) { memcpy(w, p->buffer, sizeof(ma_hit_t)*p->length); w += sizeof(ma_hit_t)*p->length; }

	return (uint64_t)(w - buf);
}

static uint64_t ov_payload(const ha_ckpt_stream_t *s, uint64_t rid, uint8_t *buf)
{
	uint8_t *w = buf;
	uint64_t n;
	ma_hit_t_alloc *p = &R_INF.paf[rid];
	(void)s;
	n = p->length; memcpy(w, &n, 8); w += 8;
	if (p->length) { memcpy(w, p->buffer, sizeof(ma_hit_t)*p->length); w += sizeof(ma_hit_t)*p->length; }
	p = &R_INF.reverse_paf[rid];
	n = p->length; memcpy(w, &n, 8); w += 8;
	if (p->length) { memcpy(w, p->buffer, sizeof(ma_hit_t)*p->length); w += sizeof(ma_hit_t)*p->length; }
	return (uint64_t)(w - buf);
}

/* replay one EC record into global slots */
static int ec_apply(ha_ckpt_stream_t *s, const uint8_t *r, uint64_t plen, uint64_t rid)
{
	const uint8_t *p = r, *e = r + plen;
	uint64_t n; uint8_t f;
	ma_hit_t_alloc *pa;

#define NEED(k) do { if ((uint64_t)(e-p) < (k)) return 0; } while (0)
	NEED(1); R_INF.trio_flag[rid] = *p++;
	NEED(4); { uint32_t v; memcpy(&v, p, 4); p += 4; if (het_cnt) het_cnt[rid] = v; }

	NEED(8); memcpy(&n, p, 8); p += 8;
	NEED(1); f = *p++;
	NEED(2*n);
	if (scc.a) {
		if (n > scc.a[rid].m) { REALLOC(scc.a[rid].a, n); scc.a[rid].m = n; }
		scc.a[rid].n = n;
		if (n) memcpy(scc.a[rid].a, p, 2*n);
		if (scc.f) scc.f[rid] = f;
	}
	p += 2*n;

	NEED(8); memcpy(&n, p, 8); p += 8;
	NEED(2*n);
	if (scb.a) {
		if (n > scb.a[rid].m) { REALLOC(scb.a[rid].a, n); scb.a[rid].m = n; }
		scb.a[rid].n = n;
		if (n) memcpy(scb.a[rid].a, p, 2*n);
	}
	p += 2*n;

	NEED(8); memcpy(&n, p, 8); p += 8;
	NEED(2*n);
	if (sca.a) {
		if (n > sca.a[rid].m) { REALLOC(sca.a[rid].a, n); sca.a[rid].m = n; }
		sca.a[rid].n = n;
		if (n) memcpy(sca.a[rid].a, p, 2*n);
	}
	p += 2*n;

	pa = &R_INF.paf[rid];
	NEED(2); pa->is_fully_corrected = p[0]; pa->is_abnormal = p[1]; p += 2;
	NEED(8); memcpy(&n, p, 8); p += 8;
	NEED(n*sizeof(ma_hit_t));
	if (pa->buffer) free(pa->buffer);
	pa->length = pa->size = n;
	pa->buffer = NULL;
	if (n) {
		pa->buffer = (ma_hit_t*)malloc(sizeof(ma_hit_t)*n);
		if (!pa->buffer) return 0;
		memcpy(pa->buffer, p, sizeof(ma_hit_t)*n);
	}
	p += sizeof(ma_hit_t)*n;

	pa = &R_INF.reverse_paf[rid];
	NEED(2); pa->is_fully_corrected = p[0]; pa->is_abnormal = p[1]; p += 2;
	NEED(8); memcpy(&n, p, 8); p += 8;
	NEED(n*sizeof(ma_hit_t));
	if (pa->buffer) free(pa->buffer);
	pa->length = pa->size = n;
	pa->buffer = NULL;
	if (n) {
		pa->buffer = (ma_hit_t*)malloc(sizeof(ma_hit_t)*n);
		if (!pa->buffer) return 0;
		memcpy(pa->buffer, p, sizeof(ma_hit_t)*n);
	}
	p += sizeof(ma_hit_t)*n;
#undef NEED
	return (p == e);
}

static int ov_apply(ha_ckpt_stream_t *s, const uint8_t *r, uint64_t plen, uint64_t rid)
{
	const uint8_t *p = r, *e = r + plen;
	ma_hit_t_alloc *pa;
	uint64_t n;
	(void)s;
#define NEED(k) do { if ((uint64_t)(e-p) < (k)) return 0; } while (0)
	pa = &R_INF.paf[rid];
	NEED(8); memcpy(&n, p, 8); p += 8;
	NEED(n*sizeof(ma_hit_t));
	if (pa->buffer) free(pa->buffer);
	pa->length = pa->size = n;
	pa->buffer = NULL;
	if (n) { pa->buffer = (ma_hit_t*)malloc(sizeof(ma_hit_t)*n); if (!pa->buffer) return 0; memcpy(pa->buffer, p, sizeof(ma_hit_t)*n); }
	p += sizeof(ma_hit_t)*n;

	pa = &R_INF.reverse_paf[rid];
	NEED(8); memcpy(&n, p, 8); p += 8;
	NEED(n*sizeof(ma_hit_t));
	if (pa->buffer) free(pa->buffer);
	pa->length = pa->size = n;
	pa->buffer = NULL;
	if (n) { pa->buffer = (ma_hit_t*)malloc(sizeof(ma_hit_t)*n); if (!pa->buffer) return 0; memcpy(pa->buffer, p, sizeof(ma_hit_t)*n); }
	p += sizeof(ma_hit_t)*n;
#undef NEED
	return (p == e);
}

/* scan an existing stream; apply good records; sets s->resume_off to the
   offset where the valid log ends (start of first bad/partial record) */
static int stream_rescan(ha_ckpt_stream_t *s, FILE *fp)
{
	int64_t round = -1; uint64_t n_reads = 0;
	uint8_t *payload = NULL;
	uint64_t payload_cap = 0;
	long good;
	int rh;

	rh = rd_hdr(fp, s->stage, &round, &n_reads);
	if (rh <= 0 || round != s->round || n_reads != s->n_reads) return 0;
	good = ftell(fp);
	uint64_t fsz = 0;
	{ struct stat st; if (fstat(fileno(fp), &st) == 0) fsz = (uint64_t)st.st_size; }

	for (;;) {
		uint16_t magic; uint32_t stage; uint64_t rid, plen; uint32_t crc;
		if (fread(&magic, sizeof(magic), 1, fp) != 1 || magic != HA_CKPT_REC_MAGIC) break;
		if (fread(&stage, sizeof(stage), 1, fp) != 1 || stage != s->stage) break;
		if (fread(&rid, sizeof(rid), 1, fp) != 1 || rid >= s->n_reads) break;
		if (fread(&plen, sizeof(plen), 1, fp) != 1) break;
		if (fsz && plen > fsz) break;
		if (plen > payload_cap) {
			uint8_t *nb = (uint8_t*)realloc(payload, plen);
			if (!nb) break;
			payload = nb; payload_cap = plen;
		}
		if (fread(payload, 1, plen, fp) != (size_t)plen) break;
		if (fread(&crc, sizeof(crc), 1, fp) != 1) break;
		if (crc != (uint32_t)crc32(0, payload, (uInt)plen)) break;
		if ((s->stage == HA_CKPT_STAGE_EC ? ec_apply(s, payload, plen, rid)
		                                : ov_apply(s, payload, plen, rid)) == 0) break;
		s->done[rid] = 1;
		good = ftell(fp);   /* record fully verified & applied */
	}
	s->n_done = 0;
	{ uint64_t i; for (i = 0; i < s->n_reads; ++i) s->n_done += s->done[i]; }
	free(payload);
	s->resume_off = good;
	return 1;
}

ha_ckpt_stream_t *ha_ckpt_stream_open(int stage, int64_t round, uint64_t n_reads)
{
	ha_ckpt_stream_t *s;
	char *path;
	FILE *old;

	if (!ha_ckpt_on()) return NULL;

	s = (ha_ckpt_stream_t*)calloc(1, sizeof(*s));
	s->stage = (uint8_t)stage;
	s->round = (int32_t)round;
	s->n_reads = n_reads;
	s->fsync_sec = asm_opt.ckpt_int;
	s->last_fsync = yak_realtime();
	pthread_mutex_init(&s->mut, NULL);
	s->done = (uint8_t*)calloc(n_reads, 1);
	s->params_hash = ha_ckpt_params_hash();
	path = stream_name(stage, round);
	s->path = path;

	old = fopen(path, "r");
	if (old && stream_rescan(s, old)) {
		fclose(old);
		truncate(path, s->resume_off);   /* drop torn tail */
		fprintf(stderr, "[M::%s::%.3f] ==> resumed %lu/%lu reads from %s\n",
			__func__, yak_realtime(), (unsigned long)s->n_done,
			(unsigned long)n_reads, path);
	} else {
		if (old) fclose(old);
		unlink(path);   /* incompatible header -> start over */
		s->n_done = 0;
		memset(s->done, 0, n_reads);
	}

	s->fp = fopen(path, "a");
	if (!s->fp) { free(s->done); free(path); pthread_mutex_destroy(&s->mut); free(s); return NULL; }
	if (s->n_done == 0) wr_hdr(s->fp, (uint8_t)stage, round, n_reads);
	return s;
}

static void stream_append(ha_ckpt_stream_t *s, uint64_t rid, const uint8_t *payload, uint64_t plen)
{
	uint16_t magic = HA_CKPT_REC_MAGIC;
	uint32_t stage = s->stage;
	uint32_t crc = (uint32_t)crc32(0, payload, (uInt)plen);

	int dosync;
	pthread_mutex_lock(&s->mut);
	fwrite(&magic, sizeof(magic), 1, s->fp);
	fwrite(&stage, sizeof(stage), 1, s->fp);
	fwrite(&rid, sizeof(rid), 1, s->fp);
	fwrite(&plen, sizeof(plen), 1, s->fp);
	fwrite(payload, 1, plen, s->fp);
	fwrite(&crc, sizeof(crc), 1, s->fp);
	fflush(s->fp);                        /* kill-safe: pushes to page cache */
	dosync = (s->fsync_sec > 0 && yak_realtime() - s->last_fsync >= s->fsync_sec);
	if (dosync) s->last_fsync = yak_realtime();
	pthread_mutex_unlock(&s->mut);
	/* fsync outside the lock: it covers whatever is dirty at call time,
	   so a concurrent append only adds to what gets synced */
	if (dosync) fsync(fileno(s->fp));
}

void ha_ckpt_dump_ec(ha_ckpt_stream_t *s, uint64_t rid)
{
	uint64_t need; uint8_t *buf;
	if (!s || s->done[rid]) return;
	need = 64 + 4*(scc.a?scc.a[rid].n:0) + 4*(scb.a?scb.a[rid].n:0)
		+ 4*(sca.a?sca.a[rid].n:0)
		+ (uint64_t)sizeof(ma_hit_t)*(R_INF.paf[rid].length + R_INF.reverse_paf[rid].length);
	buf = (uint8_t*)malloc(need);
	if (!buf) return;
	stream_append(s, rid, buf, ec_payload(s, rid, buf));
	free(buf);
	s->done[rid] = 1;
}

void ha_ckpt_dump_ov(ha_ckpt_stream_t *s, uint64_t rid)
{
	uint64_t need; uint8_t *buf;
	if (!s || s->done[rid]) return;
	need = 64 + (uint64_t)sizeof(ma_hit_t)*(R_INF.paf[rid].length + R_INF.reverse_paf[rid].length);
	buf = (uint8_t*)malloc(need);
	if (!buf) return;
	stream_append(s, rid, buf, ov_payload(s, rid, buf));
	free(buf);
	s->done[rid] = 1;
}

void ha_ckpt_stream_close(ha_ckpt_stream_t *s)
{
	if (!s) return;
	fflush(s->fp);
	if (s->fsync_sec >= 0) fsync(fileno(s->fp));
	fclose(s->fp);
	pthread_mutex_destroy(&s->mut);
	free(s->done);
	free(s->path);
	free(s);
}

void ha_ckpt_purge(void)
{
	char *n;
	int64_t r;
	/* sweep a fixed window: a previous run may have used a different -r */
	for (r = 0; r < 64; ++r) {
		n = ckpt_name(asm_opt.output_file_name, "%s.r%lld.snap", r); unlink(n); free(n);
		n = ckpt_name(asm_opt.output_file_name, "%s.r%lld.snap.tmp", r); unlink(n); free(n);
		n = ckpt_name(asm_opt.output_file_name, "%s.r%lld.ecstr", r); unlink(n); free(n);
	}
	n = ckpt_name(asm_opt.output_file_name, "%s.ff.ovstr", -1); unlink(n); free(n);
	/* stale .tmp litter from runs killed mid-atomic-write */
	{	char b[4096];
		snprintf(b, sizeof(b), "%s.ec.bin.tmp", asm_opt.output_file_name); unlink(b);
		snprintf(b, sizeof(b), "%s.ovlp.source.bin.tmp", asm_opt.output_file_name); unlink(b);
		snprintf(b, sizeof(b), "%s.ovlp.reverse.bin.tmp", asm_opt.output_file_name); unlink(b); }
}
