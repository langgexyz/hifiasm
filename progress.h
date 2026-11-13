#ifndef HA_PROGRESS_H
#define HA_PROGRESS_H
#include <stdint.h>
#include <pthread.h>

/* lightweight progress/ETA reporter for kt_for-style loops.
 *
 * ha_prog_begin_run()  -- call once at program start (marks t=0)
 * ha_prog_stage(name,w)-- label a non-counted stage; w = its rough
 *                         share of total runtime (see weights below)
 * ha_prog_init(p,nm,tot,w) -- start a counted stage (one per kt_for)
 * ha_prog_add(p,n)     -- workers bump the counter; prints every ~10%
 * ha_prog_close(p)     -- stage done
 *
 * overall eta = elapsed_total * (1 - W) / W where
 * W = finished stage weights + current stage weight * stage fraction.
 * Weights are approximate (measured on HiFi 1Gb -t16) so the total eta
 * is an estimate, while the per-stage eta is exact-rate.
 */

typedef struct {
	char name[80];
	uint64_t total, next;
	volatile uint64_t done;
	volatile int init;
	double t0;
	double weight;
	pthread_mutex_t mut;
} ha_prog_t;

#ifdef __cplusplus
extern "C" {
#endif

void ha_prog_begin_run(void);
void ha_prog_stage(const char *name, double weight);
/* label a sub-step inside a stage; does not consume weight */
void ha_prog_sub(const char *name);
/* account weight of stages already done (e.g. checkpoint resume) */
void ha_prog_skip(double weight);
/* weight < 0: stage was already declared via ha_prog_stage -- reuse it */
void ha_prog_init(ha_prog_t *p, const char *name, uint64_t total, double weight);
void ha_prog_add(ha_prog_t *p, uint64_t n);
void ha_prog_close(ha_prog_t *p);
/* rough weight of error-correction round r of nr (geometric ~0.75 decay) */
double ha_prog_round_weight(int r, int nr);

#ifdef __cplusplus
}
#endif
#endif
