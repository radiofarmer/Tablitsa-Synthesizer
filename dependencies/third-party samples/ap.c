// ----------------------------------------------------------------------------
// ap~ : first order allpass filter chain (PD external)
// Jari Kleimola 2009
// see http://www.acoustics.hut.fi/publications/papers/dafx09-cm/
// 
#define __attribute__(x)	;
#include "m_pd.h"

#define c_maxStages		2000

typedef struct _stage
{
	t_float x1;
	t_float y1;
} t_stage;

static t_class* ap_class;
typedef struct _ap
{
	t_object x_obj;
	t_float x_x;
	t_int numStages;
	t_stage stage[c_maxStages];
} t_ap;

// ----------------------------------------------------------------------------
//
static t_floatarg clamp(t_floatarg var, t_floatarg min, t_floatarg max)
{
	if (var < min)			var = min;
	else if (var > max)	var = max;
	return var;
}

// ----------------------------------------------------------------------------
// ctor
//
static void* ap_new()
{
	t_int s;
	t_ap* ap = (t_ap*)pd_new(ap_class);
	ap->numStages = 1;
	for (s = 0; s < c_maxStages; s++)
	{
		ap->stage[s].x1 = 0;
		ap->stage[s].y1 = 0;
	}
	inlet_new(&ap->x_obj, &ap->x_obj.ob_pd, &s_signal, &s_signal);
	inlet_new(&ap->x_obj, &ap->x_obj.ob_pd, &s_float, gensym("numStages"));
	outlet_new(&ap->x_obj, &s_signal);
	return ap;
}

// ----------------------------------------------------------------------------
//
static void ap_dtor()
{
}

// ----------------------------------------------------------------------------
// inlets
//
static void ap_numStages(t_ap* ap, t_floatarg n)
{
	ap->numStages = clamp(n, 1, 2000);
}

// ----------------------------------------------------------------------------
// the business
//
static t_int* ap_process(t_int* w)
{
	t_ap* ap = (t_ap*)w[1];
	t_sample* in  = (t_float*)w[2];
	t_sample* a1  = (t_float*)w[3];
	t_sample* out = (t_float*)w[4];
	int n = (int)w[5];
	int numStages = ap->numStages;
	int i,s;
	t_float x,a,y;

	for (i = 0; i<n; i++)
	{
		t_stage* stage = &ap->stage[0];
		x = *in++;
		a = *a1++;
		for (s = 0; s < numStages; s++)
		{
			y = stage->x1 + a * (x - stage->y1);
			stage->x1 = x;
			stage->y1 = y;
			x = y;
			stage++;
		}
		*out++ = y;
	}

	// next DSP callback object (for PD)
	return (w+6);
}

// ----------------------------------------------------------------------------
// called when DSP is turned on. Registers the actual processing method.
//
static void ap_dsp(t_ap* ap, t_signal** sp)
{
	dsp_add(ap_process, 5, ap, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[0]->s_n);
}

// ----------------------------------------------------------------------------
// entry point 
//
void ap_tilde_setup(void)
{
	ap_class = class_new(gensym("ap~"), (t_newmethod)ap_new, (t_method)ap_dtor, sizeof(t_ap), CLASS_DEFAULT, 0);
	CLASS_MAINSIGNALIN(ap_class, t_ap, x_x);
	class_addmethod(ap_class, (t_method)ap_dsp, gensym("dsp"), 0);
	class_addmethod(ap_class, (t_method)ap_numStages, gensym("numStages"), A_FLOAT, 0);
}
