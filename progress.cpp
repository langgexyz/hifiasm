#include <stdio.h>
#include <string.h>
#include <math.h>
#include "progress.h"

extern double yak_realtime(void);

/* whole-run state */
static double g_trun = 0;     /* run start time */
static double g_wdone = 0;    /* weight of finished stages */
static double g_wcur = 0;     /* weight of the current stage */

void ha_prog_begin_run(void)
{
	g_trun = yak_realtime();
	g_wdone = g_wcur = 0;
}

/* move w from "current" to "done"; call when a new stage begins */
static void roll_stage(double w)
{
	g_wdone += g_wcur;
	g_wcur = w;
}

void ha_prog_stage(const char *name, double weight)
{
	roll_stage(weight);
	if (g_trun)
		fprintf(stderr, "[M::stage::%.1f] ==> %s\n", yak_realtime() - g_trun, name);
	else
		fprintf(stderr, "[M::stage] ==> %s\n", name);
}

void ha_prog_skip(double weight)
{
	/* account for a stage finished before this run (checkpoint resume) */
	g_wdone += weight;
}

void ha_prog_sub(const char *name)
{
	if (g_trun)
		fprintf(stderr, "[M::stage::%.1f] ==> %s\n", yak_realtime() - g_trun, name);
}

void ha_prog_init(ha_prog_t *p, const char *name, uint64_t total, double weight)
{
	if (weight >= 0) roll_stage(weight);
	memset(p, 0, sizeof(*p));
	snprintf(p->name, sizeof(p->name), "%s", name);
	p->total = total ? total : 1;
	p->next = p->total / 10 ? p->total / 10 : 1;
	p->weight = (weight >= 0) ? weight : g_wcur;
	p->t0 = yak_realtime();
	pthread_mutex_init(&p->mut, NULL);
	p->init = 1;
}

static void fmt_eta(char *b, size_t n, double sec)
{
	if (sec < 90) snprintf(b, n, "%.0fs", sec);
	else if (sec < 5400) snprintf(b, n, "%.1fm", sec/60);
	else snprintf(b, n, "%.1fh", sec/3600);
}

void ha_prog_add(ha_prog_t *p, uint64_t n)
{
	uint64_t d;
	if (!p->init) return;
	d = (uint64_t)__sync_add_and_fetch(
			(volatile long long*)&p->done, (long long)n);
	if (d < p->next) return;
	pthread_mutex_lock(&p->mut);
	if (p->done >= p->next) {
		double now = yak_realtime();
		double el = now - p->t0;
		double frac = (double)p->done / p->total;
		if (frac > 1) frac = 1;
		double seta = p->done ? el * (p->total - p->done) / p->done : 0;
		char s[32], t[32], o[64] = "";
		fmt_eta(s, sizeof(s), seta);
		/* overall estimate from stage weights */
		if (g_trun) {
			double wfrac = g_wdone + p->weight * frac;
			if (wfrac > 0.001) {
				double oeta = (now - g_trun) * (1 - wfrac) / wfrac;
				fmt_eta(t, sizeof(t), oeta > 0 ? oeta : 0);
				snprintf(o, sizeof(o), " | run ~%.0f%% eta %s", wfrac*100, t);
			}
		}
		fprintf(stderr, "[M::prog::%.1f] %s %llu/%llu %.0f%% | eta %s%s\n",
			el, p->name, (unsigned long long)p->done,
			(unsigned long long)p->total, frac*100, s, o);
		do { p->next += p->total/10 ? p->total/10 : 1; } while (p->next <= p->done);
	}
	pthread_mutex_unlock(&p->mut);
}

void ha_prog_close(ha_prog_t *p)
{
	if (!p->init) return;
	p->init = 0;
	fprintf(stderr, "[M::prog::%.1f] %s %llu/%llu 100%% | done\n",
		yak_realtime() - p->t0, p->name,
		(unsigned long long)p->done, (unsigned long long)p->total);
	pthread_mutex_destroy(&p->mut);
}

double ha_prog_round_weight(int r, int nr)
{
	/* error-correction rounds share ~60% of total runtime; each round
	   takes roughly 0.75x the previous one (measured 1Gb -t16) */
	double sum = 0;
	int i;
	for (i = 0; i < nr; ++i) sum += pow(0.75, i);
	return 0.60 * pow(0.75, r) / (sum ? sum : 1);
}
