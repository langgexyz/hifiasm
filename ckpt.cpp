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

#define CKPT_W(fp, p, n, c) do { \
	fwrite((p), 1, (n), (fp)); \
	if (c) *(c) = crc32(*(c), (const Bytef*)(p), (uInt)(n)); \
} while (0)

#define CKPT_RD(fp, p, n, c) do { \
	if (fread((p), 1, (n), (fp)) != (size_t)(n)) return 0; \
	if (c) *(c) = crc32(*(c), (const Bytef*)(p), (uInt)(n)); \
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
	h = fnv1a(h, &o->mz_rewin, sizeof(o->mz_rewin));
	h = fnv1a(h, &o->mz_sample_dist, sizeof(o->mz_sample_dist));
	h = fnv1a(h, &o->max_kmer_cnt, sizeof(o->max_kmer_cnt));
	h = fnv1a(h, &o->min_hist_kmer_cnt, sizeof(o->min_hist_kmer_cnt));
	h = fnv1a(h, &o->is_dbg_het_cnt, sizeof(o->is_dbg_het_cnt));
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
	int fd;
	fflush(fp);
	fd = fileno(fp);
	if (fd >= 0) fsync(fd);
	if (fclose(fp) != 0) { unlink(tmp); free((void*)tmp); return 0; }
	if (rename(tmp, path) != 0) { unlink(tmp); free((void*)tmp); return 0; }
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

static int rd_paf_arr(FILE *fp, ma_hit_t_alloc *paf, uint64_t n, uint32_t *crc)
{
	uint64_t i;
	for (i = 0; i < n; ++i) {
		CKPT_RD(fp, &paf[i].is_fully_corrected, sizeof(uint8_t), crc);
		CKPT_RD(fp, &paf[i].is_abnormal, sizeof(uint8_t), crc);
		CKPT_RD(fp, &paf[i].length, sizeof(paf[i].length), crc);
		paf[i].size = paf[i].length;
		paf[i].buffer = NULL;
		if (paf[i].length) {
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

static int rd_cc(FILE *fp, cc_v *cc, uint64_t n, uint32_t *crc)
{
	uint64_t i;
	for (i = 0; i < n; ++i) {
		uint64_t m = 0; uint8_t f = 0;
		CKPT_RD(fp, &m, sizeof(m), crc);
		CKPT_RD(fp, &f, sizeof(f), crc);
		if (m > cc->a[i].m) { REALLOC(cc->a[i].a, m); cc->a[i].m = m; }
		cc->a[i].n = m;
		if (m) CKPT_RD(fp, cc->a[i].a, sizeof(uint16_t)*m, crc);
		if (cc->f) cc->f[i] = f;
	}
	return 1;
}

static int ckpt_write_one(const char *path, int64_t round, void *flt_tab)
{
	char *tmp = NULL;
	FILE *fp = ckpt_fopen_atomic(path, &tmp);
	uint32_t crc = 0;
	uint64_t u64;
	uint8_t u8;
	if (!fp) { free(tmp); return 0; }

	wr_hdr(fp, HA_CKPT_STAGE_SNAP, round, R_INF.total_reads);

	u8 = (flt_tab != NULL); CKPT_W(fp, &u8, 1, &crc);
	if (u8) { if (!ckpt_ft_save(flt_tab, fp)) { fclose(fp); unlink(tmp); free(tmp); return 0; } }
	/* flt_tab body is opaque; fold it into crc by hashing file bytes later is
	   not trivial -- crc covers the sections we serialize ourselves, the
	   FIN_MAGIC footer still guarantees the file is complete */

	write_All_reads_fp(&R_INF, fp, &crc);

	wr_paf_arr(fp, R_INF.paf, R_INF.total_reads, &crc);
	wr_paf_arr(fp, R_INF.reverse_paf, R_INF.total_reads, &crc);

	u8 = (scc.a != NULL); CKPT_W(fp, &u8, 1, &crc);
	if (u8) wr_cc(fp, &scc, R_INF.total_reads, &crc);
	u8 = (scb.a != NULL); CKPT_W(fp, &u8, 1, &crc);
	if (u8) wr_cc(fp, &scb, R_INF.total_reads, &crc);
	u8 = (sca.a != NULL); CKPT_W(fp, &u8, 1, &crc);
	if (u8) wr_cc(fp, &sca, R_INF.total_reads, &crc);

	u64 = crc;   fwrite(&u64, sizeof(u64), 1, fp);
	u64 = HA_CKPT_FIN_MAGIC; fwrite(&u64, sizeof(u64), 1, fp);

	return ckpt_rename_sync(fp, tmp, path);
}

int ha_ckpt_write_round(int64_t round, void *flt_tab)
{
	char *path, *old;
	int ok;
	if (!ha_ckpt_on()) return 0;
	path = ckpt_name(asm_opt.output_file_name, "%s.r%lld.snap", round);
	ok = ckpt_write_one(path, round, flt_tab);
	fprintf(stderr, "[M::%s::%.3f] %s checkpoint for round %lld\n", __func__,
		yak_realtime(), ok?"==> wrote":"!! failed to write", (long long)round);
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

static int ckpt_read_one(const char *path, int64_t exp_round, void **flt_tab, All_reads *tmp)
{
	FILE *fp = fopen(path, "r");
	int64_t round = -1; uint64_t n_reads = 0, u64;
	uint32_t crc = 0, file_crc;
	uint8_t u8;
	int rh;
	if (!fp) return 0;
	rh = rd_hdr(fp, HA_CKPT_STAGE_SNAP, &round, &n_reads);
	if (rh <= 0 || round != exp_round || n_reads == 0) { fclose(fp); fprintf(stderr, "[ckpt-dbg] %s hdr fail rh=%d round=%lld n=%lu\n", path, rh, (long long)round, (unsigned long)n_reads); return 0; }

	CKPT_RD(fp, &u8, 1, &crc);
	if (u8) { if (!ckpt_ft_load(flt_tab, fp)) { fclose(fp); fprintf(stderr, "[ckpt-dbg] %s ft fail\n", path); return 0; } }

	if (!load_All_reads_fp(tmp, fp, &crc)) { fclose(fp); fprintf(stderr, "[ckpt-dbg] %s reads fail\n", path); return 0; }

	tmp->paf = (ma_hit_t_alloc*)calloc(n_reads, sizeof(ma_hit_t_alloc));
	tmp->reverse_paf = (ma_hit_t_alloc*)calloc(n_reads, sizeof(ma_hit_t_alloc));
	if (!tmp->paf || !tmp->reverse_paf) { fclose(fp); return 0; }
	if (!rd_paf_arr(fp, tmp->paf, n_reads, &crc)) { fclose(fp); fprintf(stderr, "[ckpt-dbg] %s paf fail\n", path); return 0; }
	if (!rd_paf_arr(fp, tmp->reverse_paf, n_reads, &crc)) { fclose(fp); fprintf(stderr, "[ckpt-dbg] %s rpaf fail\n", path); return 0; }

	CKPT_RD(fp, &u8, 1, &crc);
	if (u8) {
		scc.n = scc.m = n_reads; CALLOC(scc.a, n_reads); CALLOC(scc.f, n_reads);
		if (!rd_cc(fp, &scc, n_reads, &crc)) { fclose(fp); fprintf(stderr, "[ckpt-dbg] %s scc fail\n", path); return 0; }
	}
	CKPT_RD(fp, &u8, 1, &crc);
	if (u8) {
		scb.n = scb.m = n_reads; CALLOC(scb.a, n_reads);
		if (!rd_cc(fp, &scb, n_reads, &crc)) { fclose(fp); fprintf(stderr, "[ckpt-dbg] %s scb fail\n", path); return 0; }
	}
	CKPT_RD(fp, &u8, 1, &crc);
	if (u8) {
		sca.n = sca.m = n_reads; CALLOC(sca.a, n_reads);
		if (!rd_cc(fp, &sca, n_reads, &crc)) { fclose(fp); fprintf(stderr, "[ckpt-dbg] %s sca fail\n", path); return 0; }
	}

	if (fread(&u64, sizeof(u64), 1, fp) != 1) { fclose(fp); fprintf(stderr, "[ckpt-dbg] %s crc missing\n", path); return 0; }
	file_crc = (uint32_t)u64;
	if (fread(&u64, sizeof(u64), 1, fp) != 1 || u64 != HA_CKPT_FIN_MAGIC) { fclose(fp); fprintf(stderr, "[ckpt-dbg] %s fin fail\n", path); return 0; }
	if (file_crc != crc) { fclose(fp); fprintf(stderr, "[ckpt-dbg] %s crc mismatch %x vs %x\n", path, file_crc, crc); return 0; }
	fclose(fp);
	if (tmp->total_reads != n_reads) { fprintf(stderr, "[ckpt-dbg] %s n_reads mismatch\n", path); return 0; }
	return 1;
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
	if (n) { pa->buffer = (ma_hit_t*)malloc(sizeof(ma_hit_t)*n); memcpy(pa->buffer, p, sizeof(ma_hit_t)*n); }
	p += sizeof(ma_hit_t)*n;

	pa = &R_INF.reverse_paf[rid];
	NEED(8); memcpy(&n, p, 8); p += 8;
	NEED(n*sizeof(ma_hit_t));
	if (pa->buffer) free(pa->buffer);
	pa->length = pa->size = n;
	pa->buffer = NULL;
	if (n) { pa->buffer = (ma_hit_t*)malloc(sizeof(ma_hit_t)*n); memcpy(pa->buffer, p, sizeof(ma_hit_t)*n); }
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

	for (;;) {
		uint16_t magic; uint32_t stage; uint64_t rid, plen; uint32_t crc;
		if (fread(&magic, sizeof(magic), 1, fp) != 1 || magic != HA_CKPT_REC_MAGIC) break;
		if (fread(&stage, sizeof(stage), 1, fp) != 1 || stage != s->stage) break;
		if (fread(&rid, sizeof(rid), 1, fp) != 1 || rid >= s->n_reads) break;
		if (fread(&plen, sizeof(plen), 1, fp) != 1) break;
		if (plen > payload_cap) {
			payload = (uint8_t*)realloc(payload, plen);
			if (!payload) break;
			payload_cap = plen;
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
		if (s->n_done) fprintf(stderr, "[M::%s] !! stale %s ignored\n", __func__, path);
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

	pthread_mutex_lock(&s->mut);
	fwrite(&magic, sizeof(magic), 1, s->fp);
	fwrite(&stage, sizeof(stage), 1, s->fp);
	fwrite(&rid, sizeof(rid), 1, s->fp);
	fwrite(&plen, sizeof(plen), 1, s->fp);
	fwrite(payload, 1, plen, s->fp);
	fwrite(&crc, sizeof(crc), 1, s->fp);
	fflush(s->fp);                        /* kill-safe: pushes to page cache */
	if (s->fsync_sec > 0 && yak_realtime() - s->last_fsync >= s->fsync_sec) {
		fsync(fileno(s->fp));
		s->last_fsync = yak_realtime();
	}
	pthread_mutex_unlock(&s->mut);
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
	free(s->wbuf);
	free(s->path);
	free(s);
}

void ha_ckpt_purge(void)
{
	char *n;
	int64_t r;
	for (r = -1; r < (int64_t)asm_opt.number_of_round + 1; ++r) {
		if (r >= 0) {
			n = ckpt_name(asm_opt.output_file_name, "%s.r%lld.snap", r); unlink(n); free(n);
			n = ckpt_name(asm_opt.output_file_name, "%s.r%lld.ecstr", r); unlink(n); free(n);
		}
	}
	n = ckpt_name(asm_opt.output_file_name, "%s.ff.ovstr", -1); unlink(n); free(n);
}
