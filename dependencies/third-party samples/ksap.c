// ----------------------------------------------------------------------------
// ksap~ extended karplus-strong with allpass chain inside the loop (PD external)
// Jari Kleimola 2009
// see http://www.acoustics.hut.fi/publications/papers/dafx09-cm/
//
#define __attribute__(x)	;
#include "m_pd.h"
#include <stdlib.h>
#include <time.h>

#define c_iMaxStringLength	4096
#define c_maxStages		2000

typedef struct _stage
{
	t_float x1;
	t_float y1;
} t_stage;

static t_class* ksap_class;
typedef struct _ksap
{
	t_object x_obj;
	t_float x_x;

	int m_phase;
	int m_iDelayLength;
	t_float m_f0;
	t_float m_stretch,m_damp;
	t_float m_prevOut;
	t_float m_delay[2048];

	t_int numStages;
	t_stage stage[c_maxStages];
} t_ksap;

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
static void* ksap_new()
{
	t_int s;
	t_ksap* ksap = (t_ksap*)pd_new(ksap_class);
	ksap->x_x = 0;

	// -- string
	ksap->m_phase = 0;
	ksap->m_f0 = 130;
	ksap->m_stretch = ksap->m_damp = 0;
	ksap->m_iDelayLength = 0;
	ksap->m_prevOut = 0;
	outlet_new(&ksap->x_obj, &s_signal);
	srand(time(0));

	// -- allpass chain
	ksap->numStages = 1;
	for (s = 0; s < c_maxStages; s++)
	{
		ksap->stage[s].x1 = 0;
		ksap->stage[s].y1 = 0;
	}

	inlet_new(&ksap->x_obj, &ksap->x_obj.ob_pd, &s_float, gensym("f0"));
	inlet_new(&ksap->x_obj, &ksap->x_obj.ob_pd, &s_float, gensym("damp"));
	inlet_new(&ksap->x_obj, &ksap->x_obj.ob_pd, &s_float, gensym("stretch"));
	inlet_new(&ksap->x_obj, &ksap->x_obj.ob_pd, &s_float, gensym("numStages"));

	return ksap;
}

// ----------------------------------------------------------------------------
//
static void ksap_dtor()
{
}

// ----------------------------------------------------------------------------
// inlets
//
static void ksap_damp(t_ksap* ksap, t_floatarg f)			{ ksap->m_damp = clamp(f, 0, 1);	}
static void ksap_stretch(t_ksap* ksap, t_floatarg f)		{ ksap->m_stretch = clamp(f, 0, 1);	}
static void ksap_f0(t_ksap* ksap, t_floatarg f)				{ ksap->m_f0 = f;	}
static void ksap_numStages(t_ksap* ksap, t_floatarg n)	{ ksap->numStages = clamp(n, 0, 2000); }

// ----------------------------------------------------------------------------
//
void ksap_clear(t_ksap* ksap)
{
	int s;
	for (s = 0; s < c_maxStages; s++)
	{
		ksap->stage[s].x1 = 0;
		ksap->stage[s].y1 = 0;
	}
	for (s=0; s<ksap->m_iDelayLength; s++)
		ksap->m_delay[s] = 0;
}

// ----------------------------------------------------------------------------
//
void ksap_excite(t_ksap* ksap)
{
	int i;
	t_float prev = 0;
	t_float loss = 0.8;

	ksap_clear(ksap);
	ksap->m_iDelayLength = 44100 / (ksap->m_f0) - 0.5;
	if (ksap->m_iDelayLength > (c_iMaxStringLength - 1))
		ksap->m_iDelayLength = c_iMaxStringLength - 1;

	for (i=0; i<ksap->m_iDelayLength; i++)
	{
		t_float s = ((t_float)rand()) / RAND_MAX - 0.5;
		s = (1 - loss) * s + prev * loss;
		ksap->m_delay[i] = s;
		prev = s;
	}

	ksap->m_phase = 0;
	ksap->m_prevOut = 0;
}

// ----------------------------------------------------------------------------
// the business
//
static t_int* ksap_process(t_int* w)
{
	t_ksap* ksap = (t_ksap*)w[1];
	t_sample* a1  = (t_float*)w[2];
	t_float* out = (t_float*)w[3];
	int n = (int)w[4];
	int i,is;
	int numStages = ksap->numStages;
	t_float x,a,y,ks;

	for (i = 0; i<n; i++)
	{
		t_stage* stage;
		t_float x = ksap->m_delay[ksap->m_phase];

		// -- allpass
		if (numStages > 0)
		{
			stage = &ksap->stage[0];
			a = *a1++;
			for (is = 0; is < numStages; is++)
			{
				y = stage->x1 + a * (x - stage->y1);
				stage->x1 = x;
				stage->y1 = y;
				x = y;
				stage++;
			}
		}
		else y = x;

		// -- loss
		// x = (x + ksap->m_prevOut) * 0.5;
		y = ksap->m_damp * ((1.0 - ksap->m_stretch) * y + ksap->m_stretch * ksap->m_prevOut);

		ksap->m_delay[ksap->m_phase] = y;
		ksap->m_prevOut = y;
		*out++ = y;

		ksap->m_phase++;
		if (ksap->m_phase >= ksap->m_iDelayLength)
			ksap->m_phase = 0;
	}

	// next DSP callback object (for PD)
	return (w+5);
}

// ----------------------------------------------------------------------------
// called when DSP is turned on. Registers the actual processing method.
//
static void ksap_dsp(t_ksap* ksap, t_signal** sp)
{
	dsp_add(ksap_process, 4, ksap, sp[0]->s_vec, sp[1]->s_vec, sp[0]->s_n);
}

// ----------------------------------------------------------------------------
// entry point 
//
void ksap_tilde_setup(void)
{
	ksap_class = class_new(gensym("ksap~"), (t_newmethod)ksap_new, (t_method)ksap_dtor, sizeof(t_ksap), CLASS_DEFAULT, 0);
	CLASS_MAINSIGNALIN(ksap_class, t_ksap, x_x);
	class_addmethod(ksap_class, (t_method)ksap_dsp, gensym("dsp"), 0);
	class_addmethod(ksap_class, (t_method)ksap_excite, gensym("excite"), 0);
	class_addmethod(ksap_class, (t_method)ksap_clear, gensym("clear"), 0);
	class_addmethod(ksap_class, (t_method)ksap_f0, gensym("f0"), A_FLOAT, 0);
	class_addmethod(ksap_class, (t_method)ksap_damp, gensym("damp"), A_FLOAT, 0);
	class_addmethod(ksap_class, (t_method)ksap_stretch, gensym("stretch"), A_FLOAT, 0);
	class_addmethod(ksap_class, (t_method)ksap_numStages, gensym("numStages"), A_FLOAT, 0);
}
